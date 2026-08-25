"""Run one spatially sharded RVMD decomposition with torchrun.

Example:
    torchrun --standalone --nproc-per-node=4 examples_python/multigpu_rvmd.py \
        snapshots.npy --output rvmd_output --k 8 --alpha 1000 \
        --maximum-steps 1000 --checkpoint rvmd_checkpoint

The input is an S-by-T NumPy .npy array. Every rank memory-maps the file and
loads only its contiguous spatial rows. The output directory contains one
spatial-mode shard per rank and one replicated temporal result written by rank 0.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np
import torch
import torch.distributed as dist

from rvmdpy import rvmd_distributed


def balanced_range(size: int, world_size: int, rank: int) -> tuple[int, int]:
    base, remainder = divmod(size, world_size)
    start = rank * base + min(rank, remainder)
    return start, start + base + (1 if rank < remainder else 0)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", nargs="?", type=Path, help="S-by-T NumPy .npy file")
    parser.add_argument("--restart", type=Path, help="distributed checkpoint directory")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path)
    parser.add_argument("--k", type=int)
    parser.add_argument("--alpha", type=float)
    parser.add_argument("--maximum-steps", type=int, default=500)
    parser.add_argument("--tolerance", type=float, default=5e-3)
    parser.add_argument("--checkpoint-interval", type=int, default=50)
    parser.add_argument("--display-interval", type=int, default=20)
    parser.add_argument("--precision", choices=("single", "double"), default="single")
    arguments = parser.parse_args()
    if arguments.restart is None and (
        arguments.input is None or arguments.k is None or arguments.alpha is None
    ):
        parser.error("a new run requires input, --k, and --alpha")
    if arguments.restart is not None and arguments.input is not None:
        parser.error("omit input when --restart is supplied")
    return arguments


def main() -> None:
    arguments = parse_arguments()
    local_rank = int(os.environ["LOCAL_RANK"])
    torch.cuda.set_device(local_rank)
    dist.init_process_group("nccl")
    rank, world_size = dist.get_rank(), dist.get_world_size()

    if arguments.restart is None:
        mapped = np.load(arguments.input, mmap_mode="r")
        if mapped.ndim != 2:
            raise ValueError("input must be an S-by-T NumPy array")
        start, stop = balanced_range(mapped.shape[0], world_size, rank)
        local_q = np.ascontiguousarray(mapped[start:stop])
        mode, info, _ = rvmd_distributed(
            local_q,
            arguments.k,
            arguments.alpha,
            device=f"cuda:{local_rank}",
            fp_precision=arguments.precision,
            tolerance=arguments.tolerance,
            maximum_steps=arguments.maximum_steps,
            display="iter",
            display_interval=arguments.display_interval,
            checkpoint_file=arguments.checkpoint,
            checkpoint_interval=arguments.checkpoint_interval,
        )
    else:
        mode, info, _ = rvmd_distributed(
            restart=arguments.restart,
            device=f"cuda:{local_rank}",
            tolerance=arguments.tolerance,
            maximum_steps=arguments.maximum_steps,
            display="iter",
            display_interval=arguments.display_interval,
            checkpoint_file=arguments.checkpoint,
            checkpoint_interval=arguments.checkpoint_interval,
        )

    if rank == 0:
        arguments.output.mkdir(parents=True, exist_ok=True)
    dist.barrier()
    torch.save(
        {
            "spatial_start": mode.spatial_start,
            "spatial_stop": mode.spatial_stop,
            "phi": mode.phi,
        },
        arguments.output / f"phi_rank_{rank:05d}.pt",
    )
    if rank == 0:
        torch.save(
            {
                "c": mode.c,
                "omega": mode.omega,
                "energy": mode.energy,
                "omega_history": info.iteration.omega,
                "difference_history": info.iteration.difference,
            },
            arguments.output / "temporal.pt",
        )
        manifest = {
            "spatial_size": info.spatial_size,
            "sample_count": info.sample_count,
            "mode_count": info.mode_count,
            "world_size": world_size,
            "stop_reason": info.stop_reason,
            "steps": info.iteration.steps,
        }
        with (arguments.output / "manifest.json").open("w", encoding="utf-8") as stream:
            json.dump(manifest, stream, indent=2)
            stream.write("\n")
    dist.barrier()
    dist.destroy_process_group()


if __name__ == "__main__":
    main()
