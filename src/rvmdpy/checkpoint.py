"""Atomic single-device and sharded distributed RVMD checkpoints."""

from __future__ import annotations

import json
import os
import shutil
from pathlib import Path
from typing import Any

import torch
import torch.distributed as dist

from .types import RestartState

PYTHON_RESTART_VERSION = 1


def _metadata(state: RestartState) -> dict[str, Any]:
    return {
        "version": state.version,
        "spatial_size": state.spatial_size,
        "sample_count": state.sample_count,
        "mode_count": state.mode_count,
        "alpha": state.alpha,
        "tolerance": state.tolerance,
        "maximum_steps": state.maximum_steps,
        "init_freq_type": state.init_freq_type,
        "init_freq_maximum": state.init_freq_maximum,
        "device": state.device,
        "fp_precision": state.fp_precision,
        "n_dc": state.n_dc,
        "is_real_input": state.is_real_input,
        "display": state.display,
        "display_interval": state.display_interval,
        "steps": state.steps,
        "data_spectrum_norm": state.data_spectrum_norm,
        "spatial_start": state.spatial_start,
        "world_size": state.world_size,
    }


def _cpu(tensor: torch.Tensor) -> torch.Tensor:
    result = tensor.detach()
    if result.device.type == "cpu":
        return result.contiguous()
    return result.to(device="cpu", copy=True).contiguous()


def _full_payload(state: RestartState) -> dict[str, Any]:
    return {
        "metadata": _metadata(state),
        "weight": _cpu(state.weight),
        "initial_frequencies": _cpu(state.initial_frequencies),
        "omega_history": _cpu(state.omega_history),
        "difference_history": _cpu(state.difference_history),
        "coefficient_spectrum": _cpu(state.coefficient_spectrum),
        "phi": _cpu(state.phi),
        "residual": _cpu(state.residual),
    }


def _state_from_payload(payload: dict[str, Any]) -> RestartState:
    metadata = payload["metadata"]
    version = int(metadata["version"])
    if version != PYTHON_RESTART_VERSION:
        raise ValueError(
            f"unsupported Python restart version {version}; "
            f"expected {PYTHON_RESTART_VERSION}"
        )
    return RestartState(
        version=version,
        spatial_size=int(metadata["spatial_size"]),
        sample_count=int(metadata["sample_count"]),
        mode_count=int(metadata["mode_count"]),
        alpha=float(metadata["alpha"]),
        weight=payload["weight"],
        tolerance=float(metadata["tolerance"]),
        maximum_steps=int(metadata["maximum_steps"]),
        init_freq_type=int(metadata["init_freq_type"]),
        init_freq_maximum=float(metadata["init_freq_maximum"]),
        initial_frequencies=payload["initial_frequencies"],
        device=str(metadata["device"]),
        fp_precision=str(metadata["fp_precision"]),
        n_dc=int(metadata["n_dc"]),
        is_real_input=bool(metadata["is_real_input"]),
        display=str(metadata["display"]),
        display_interval=int(metadata["display_interval"]),
        steps=int(metadata["steps"]),
        omega_history=payload["omega_history"],
        difference_history=payload["difference_history"],
        coefficient_spectrum=payload["coefficient_spectrum"],
        phi=payload["phi"],
        residual=payload["residual"],
        data_spectrum_norm=float(metadata["data_spectrum_norm"]),
        spatial_start=int(metadata.get("spatial_start", 0)),
        world_size=int(metadata.get("world_size", 1)),
    )


def _atomic_torch_save(payload: dict[str, Any], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(f"{path}.tmp")
    previous = Path(f"{path}.prev")
    temporary.unlink(missing_ok=True)
    torch.save(payload, temporary)
    previous.unlink(missing_ok=True)
    if path.exists():
        os.replace(path, previous)
    os.replace(temporary, path)


def save_checkpoint(path: str | os.PathLike[str], state: RestartState) -> None:
    """Atomically save a non-distributed restart state."""

    if state.world_size != 1 or state.spatial_start != 0:
        raise ValueError("use save_distributed_checkpoint for a sharded state")
    _atomic_torch_save(_full_payload(state), Path(path))


def load_checkpoint(path: str | os.PathLike[str]) -> RestartState:
    """Load a Python restart checkpoint onto CPU.

    A file is a single-device checkpoint. A directory is a distributed
    checkpoint and is assembled for one process when no process group exists.

    Args:
        path: Single-device checkpoint file or distributed checkpoint directory.

    Returns:
        A CPU-resident restart state suitable for :func:`rvmdpy.rvmd`.
    """

    checkpoint = Path(path)
    if checkpoint.is_dir():
        return load_distributed_checkpoint(checkpoint)
    payload = torch.load(checkpoint, map_location="cpu", weights_only=True)
    return _state_from_payload(payload)


def _distributed_rank_world() -> tuple[int, int]:
    if dist.is_available() and dist.is_initialized():
        return dist.get_rank(), dist.get_world_size()
    return 0, 1


def _barrier() -> None:
    if dist.is_available() and dist.is_initialized():
        dist.barrier()


def _balanced_range(size: int, world_size: int, rank: int) -> tuple[int, int]:
    base, remainder = divmod(size, world_size)
    start = rank * base + min(rank, remainder)
    stop = start + base + (1 if rank < remainder else 0)
    return start, stop


def save_distributed_checkpoint(
    path: str | os.PathLike[str], state: RestartState
) -> None:
    """Save spatial shards in parallel to a shared checkpoint directory.

    Every rank must call this collective with the same path on shared storage.
    """

    rank, world_size = _distributed_rank_world()
    if state.world_size != world_size:
        raise ValueError("restart state world_size does not match process group")

    checkpoint = Path(path)
    temporary = Path(f"{checkpoint}.tmp")
    previous = Path(f"{checkpoint}.prev")
    if rank == 0:
        shutil.rmtree(temporary, ignore_errors=True)
        temporary.mkdir(parents=True, exist_ok=False)
    _barrier()

    local_payload = {
        "spatial_start": state.spatial_start,
        "spatial_stop": state.spatial_stop,
        "weight": _cpu(state.weight),
        "phi": _cpu(state.phi),
        "residual": _cpu(state.residual),
    }
    local_path = temporary / f"rank_{rank:05d}.pt"
    torch.save(local_payload, local_path)

    if rank == 0:
        global_payload = {
            "metadata": _metadata(state),
            "initial_frequencies": _cpu(state.initial_frequencies),
            "omega_history": _cpu(state.omega_history),
            "difference_history": _cpu(state.difference_history),
            "coefficient_spectrum": _cpu(state.coefficient_spectrum),
        }
        torch.save(global_payload, temporary / "global.pt")

    range_device = (
        torch.device("cuda", torch.cuda.current_device())
        if world_size > 1 and dist.get_backend() == "nccl"
        else torch.device("cpu")
    )
    local_range = torch.tensor(
        [state.spatial_start, state.spatial_stop],
        device=range_device,
        dtype=torch.int64,
    )
    if world_size > 1:
        gathered_ranges = [torch.empty_like(local_range) for _ in range(world_size)]
        dist.all_gather(gathered_ranges, local_range)
        starts = [
            tuple(int(item) for item in bounds.cpu()) for bounds in gathered_ranges
        ]
    else:
        starts = [(state.spatial_start, state.spatial_stop)]
    if rank == 0:
        manifest = {
            "version": PYTHON_RESTART_VERSION,
            "world_size": world_size,
            "spatial_size": state.spatial_size,
            "shards": [
                {
                    "file": f"rank_{index:05d}.pt",
                    "start": int(bounds[0]),
                    "stop": int(bounds[1]),
                }
                for index, bounds in enumerate(starts)
            ],
        }
        with (temporary / "manifest.json").open("w", encoding="utf-8") as stream:
            json.dump(manifest, stream, indent=2)
            stream.write("\n")
    _barrier()

    if rank == 0:
        shutil.rmtree(previous, ignore_errors=True)
        if checkpoint.exists():
            os.replace(checkpoint, previous)
        os.replace(temporary, checkpoint)
    _barrier()


def load_distributed_checkpoint(
    path: str | os.PathLike[str],
) -> RestartState:
    """Load and reshard a directory checkpoint for the active process group.

    Each rank reads only old shard files overlapping its new balanced row range.
    Without a process group, the full spatial state is assembled on CPU.

    Args:
        path: Directory created by a distributed RVMD checkpoint write.

    Returns:
        This rank's CPU-resident restart shard, or the assembled full state when
        no multi-rank process group is active.
    """

    checkpoint = Path(path)
    with (checkpoint / "manifest.json").open(encoding="utf-8") as stream:
        manifest = json.load(stream)
    if int(manifest["version"]) != PYTHON_RESTART_VERSION:
        raise ValueError(f"unsupported Python restart version {manifest['version']}")

    global_payload = torch.load(
        checkpoint / "global.pt", map_location="cpu", weights_only=True
    )
    metadata = global_payload["metadata"]
    rank, world_size = _distributed_rank_world()
    spatial_size = int(metadata["spatial_size"])
    start, stop = _balanced_range(spatial_size, world_size, rank)
    local_size = stop - start
    mode_count = int(metadata["mode_count"])
    spectrum_length = global_payload["coefficient_spectrum"].shape[0]
    real_dtype = global_payload["coefficient_spectrum"].real.dtype
    phi_dtype = (
        real_dtype
        if bool(metadata["is_real_input"])
        else global_payload["coefficient_spectrum"].dtype
    )
    weight = torch.empty(local_size, dtype=real_dtype)
    phi = torch.empty((local_size, mode_count), dtype=phi_dtype)
    residual = torch.empty(
        (local_size, spectrum_length),
        dtype=global_payload["coefficient_spectrum"].dtype,
    )

    covered = torch.zeros(local_size, dtype=torch.bool)
    for shard_info in manifest["shards"]:
        old_start = int(shard_info["start"])
        old_stop = int(shard_info["stop"])
        overlap_start = max(start, old_start)
        overlap_stop = min(stop, old_stop)
        if overlap_start >= overlap_stop:
            continue
        shard = torch.load(
            checkpoint / shard_info["file"], map_location="cpu", weights_only=True
        )
        source = slice(overlap_start - old_start, overlap_stop - old_start)
        target = slice(overlap_start - start, overlap_stop - start)
        weight[target] = shard["weight"][source]
        phi[target] = shard["phi"][source]
        residual[target] = shard["residual"][source]
        covered[target] = True
    if not bool(torch.all(covered)):
        raise ValueError("distributed checkpoint does not cover the spatial domain")

    payload = {
        "metadata": {
            **metadata,
            "spatial_start": start,
            "world_size": world_size,
        },
        "weight": weight,
        "initial_frequencies": global_payload["initial_frequencies"],
        "omega_history": global_payload["omega_history"],
        "difference_history": global_payload["difference_history"],
        "coefficient_spectrum": global_payload["coefficient_spectrum"],
        "phi": phi,
        "residual": residual,
    }
    return _state_from_payload(payload)
