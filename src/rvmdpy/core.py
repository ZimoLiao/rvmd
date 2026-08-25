"""PyTorch RVMD solver for CPU, CUDA, and spatially sharded multi-GPU runs."""

from __future__ import annotations

import math
import os
import time
from collections.abc import Callable
from dataclasses import dataclass
from typing import Any, Literal

import numpy as np
import torch
import torch.distributed as dist

from .checkpoint import (
    PYTHON_RESTART_VERSION,
    load_checkpoint,
    load_distributed_checkpoint,
    save_checkpoint,
    save_distributed_checkpoint,
)
from .types import Info, IterationInfo, Mode, Progress, RestartState, StopReason

OutputFunction = Callable[[Progress, Literal["init", "iter", "done"]], Any]


@dataclass(slots=True)
class _Context:
    distributed: bool
    rank: int
    world_size: int
    device: torch.device
    spatial_start: int = 0
    spatial_size: int = 0

    @property
    def primary(self) -> bool:
        return self.rank == 0

    @property
    def spatial_stop(self) -> int:
        return self.spatial_start + self.spatial_size

    def sum_(self, tensor: torch.Tensor) -> torch.Tensor:
        if self.world_size > 1:
            dist.all_reduce(tensor, op=dist.ReduceOp.SUM)
        return tensor

    def broadcast_(self, tensor: torch.Tensor, source: int = 0) -> torch.Tensor:
        if self.world_size > 1:
            dist.broadcast(tensor, src=source)
        return tensor


@dataclass(slots=True)
class _Settings:
    tolerance: float
    maximum_steps: int
    init_freq_type: int
    init_freq_maximum: float
    device: torch.device
    fp_precision: str
    n_dc: int
    initial_frequencies: torch.Tensor | None
    display: str
    display_interval: int
    output_fcn: OutputFunction | None
    time_limit: float
    checkpoint_file: str | None
    checkpoint_interval: int


def _positive_integer(value: Any, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, (int, np.integer)) or value < 1:
        raise ValueError(f"{name} must be a positive integer")
    return int(value)


def _nonnegative_integer(value: Any, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, (int, np.integer)) or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return int(value)


def _nonnegative_scalar(value: Any, name: str, *, allow_inf: bool = False) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exception:
        raise ValueError(f"{name} must be a nonnegative real scalar") from exception
    if math.isnan(result) or result < 0 or (math.isinf(result) and not allow_inf):
        raise ValueError(f"{name} must be a nonnegative real scalar")
    return result


def _resolve_device(
    device: str | torch.device | None, distributed: bool
) -> torch.device:
    if device is None:
        device = "cuda" if distributed and torch.cuda.is_available() else "cpu"
    requested = torch.device(device)
    if requested.type == "cuda":
        if not torch.cuda.is_available():
            raise RuntimeError("CUDA was requested but is not available")
        if distributed and requested.index is None:
            requested = torch.device("cuda", int(os.environ.get("LOCAL_RANK", "0")))
        elif requested.index is None:
            requested = torch.device("cuda", 0)
        torch.cuda.set_device(requested)
    elif requested.type != "cpu":
        raise ValueError("device must select a CPU or CUDA device")
    return requested


def _initialize_context(
    device: str | torch.device | None, distributed: bool
) -> _Context:
    resolved = _resolve_device(device, distributed)
    environment_world_size = int(os.environ.get("WORLD_SIZE", "1"))
    if distributed and environment_world_size > 1 and not dist.is_initialized():
        backend = "nccl" if resolved.type == "cuda" else "gloo"
        dist.init_process_group(backend=backend, init_method="env://")
    if distributed and dist.is_initialized():
        rank = dist.get_rank()
        world_size = dist.get_world_size()
        backend = dist.get_backend()
        if resolved.type == "cuda" and backend != "nccl":
            raise RuntimeError(
                "distributed CUDA execution requires an NCCL process group"
            )
        if resolved.type == "cpu" and backend == "nccl":
            raise RuntimeError("an NCCL process group requires CUDA devices")
        if resolved.type == "cuda":
            local_rank = int(os.environ.get("LOCAL_RANK", str(resolved.index or 0)))
            resolved = torch.device("cuda", local_rank)
            torch.cuda.set_device(resolved)
        return _Context(True, rank, world_size, resolved)
    return _Context(distributed, 0, 1, resolved)


def _real_complex_dtype(precision: str) -> tuple[torch.dtype, torch.dtype]:
    if precision == "single":
        return torch.float32, torch.complex64
    if precision == "double":
        return torch.float64, torch.complex128
    raise ValueError("fp_precision must be 'single' or 'double'")


def _as_tensor(value: Any, device: torch.device) -> torch.Tensor:
    if isinstance(value, torch.Tensor):
        return value.detach().to(device=device)
    return torch.as_tensor(value, device=device)


def _validate_initial_frequencies(
    frequencies: Any,
    mode_count: int,
    n_dc: int,
    real_dtype: torch.dtype,
    device: torch.device,
) -> torch.Tensor | None:
    if frequencies is None:
        return None
    raw = _as_tensor(frequencies, device)
    if (
        torch.is_complex(raw)
        or raw.ndim > 2
        or (raw.ndim == 2 and raw.shape[0] > 1 and raw.shape[1] > 1)
    ):
        raise ValueError("initial_frequencies must be a finite real K-vector")
    result = raw.to(dtype=real_dtype).reshape(-1)
    if result.numel() != mode_count or not bool(torch.isfinite(result).all()):
        raise ValueError("initial_frequencies must be a finite K-vector")
    if bool((result < 0).any()) or bool((result > 0.5).any()):
        raise ValueError("initial_frequencies must lie in [0, 0.5]")
    if n_dc and bool((result[:n_dc] != 0).any()):
        raise ValueError("the first n_dc initial frequencies must be zero")
    return result


def _new_settings(
    *,
    mode_count: int,
    tolerance: float | None,
    maximum_steps: int | None,
    init_freq_type: int | None,
    init_freq_maximum: float | None,
    device: torch.device,
    fp_precision: str | None,
    n_dc: int | None,
    initial_frequencies: Any,
    display: str | None,
    display_interval: int | None,
    output_fcn: OutputFunction | None,
    time_limit: float,
    checkpoint_file: str | os.PathLike[str] | None,
    checkpoint_interval: int,
) -> _Settings:
    if fp_precision is not None and not isinstance(fp_precision, str):
        raise ValueError("fp_precision must be 'single' or 'double'")
    precision = "single" if fp_precision is None else fp_precision.lower()
    real_dtype, _ = _real_complex_dtype(precision)
    dc_count = 0 if n_dc is None else _nonnegative_integer(n_dc, "n_dc")
    if dc_count > mode_count:
        raise ValueError("n_dc must not exceed K")
    if init_freq_type is None:
        init_type = 1
    elif isinstance(init_freq_type, bool) or not isinstance(
        init_freq_type, (int, np.integer)
    ):
        raise ValueError("init_freq_type must be -1, 0, or 1")
    else:
        init_type = int(init_freq_type)
    if init_type not in (-1, 0, 1):
        raise ValueError("init_freq_type must be -1, 0, or 1")
    init_maximum = min(
        _nonnegative_scalar(
            0.5 if init_freq_maximum is None else init_freq_maximum,
            "init_freq_maximum",
        ),
        0.5,
    )
    initial = _validate_initial_frequencies(
        initial_frequencies, mode_count, dc_count, real_dtype, device
    )
    if display is not None and not isinstance(display, str):
        raise ValueError("display must be 'off', 'final', or 'iter'")
    display_value = "off" if display is None else display.lower()
    if display_value not in ("off", "final", "iter"):
        raise ValueError("display must be 'off', 'final', or 'iter'")
    if output_fcn is not None and not callable(output_fcn):
        raise TypeError("output_fcn must be callable or None")
    return _Settings(
        tolerance=_nonnegative_scalar(
            5e-3 if tolerance is None else tolerance, "tolerance"
        ),
        maximum_steps=_positive_integer(
            500 if maximum_steps is None else maximum_steps, "maximum_steps"
        ),
        init_freq_type=init_type,
        init_freq_maximum=init_maximum,
        device=device,
        fp_precision=precision,
        n_dc=dc_count,
        initial_frequencies=initial,
        display=display_value,
        display_interval=_positive_integer(
            20 if display_interval is None else display_interval, "display_interval"
        ),
        output_fcn=output_fcn,
        time_limit=_nonnegative_scalar(time_limit, "time_limit", allow_inf=True),
        checkpoint_file=None if checkpoint_file is None else os.fspath(checkpoint_file),
        checkpoint_interval=_positive_integer(
            checkpoint_interval, "checkpoint_interval"
        ),
    )


def _restart_settings(
    state: RestartState,
    *,
    tolerance: float | None,
    maximum_steps: int | None,
    device: torch.device,
    display: str | None,
    display_interval: int | None,
    output_fcn: OutputFunction | None,
    time_limit: float,
    checkpoint_file: str | os.PathLike[str] | None,
    checkpoint_interval: int,
) -> _Settings:
    maximum = (
        state.maximum_steps
        if maximum_steps is None
        else _positive_integer(maximum_steps, "maximum_steps")
    )
    if maximum < state.steps:
        raise ValueError("maximum_steps cannot be less than completed restart steps")
    if display is not None and not isinstance(display, str):
        raise ValueError("display must be 'off', 'final', or 'iter'")
    display_value = state.display if display is None else display.lower()
    if display_value not in ("off", "final", "iter"):
        raise ValueError("display must be 'off', 'final', or 'iter'")
    if output_fcn is not None and not callable(output_fcn):
        raise TypeError("output_fcn must be callable or None")
    return _Settings(
        tolerance=state.tolerance
        if tolerance is None
        else _nonnegative_scalar(tolerance, "tolerance"),
        maximum_steps=maximum,
        init_freq_type=state.init_freq_type,
        init_freq_maximum=state.init_freq_maximum,
        device=device,
        fp_precision=state.fp_precision,
        n_dc=state.n_dc,
        initial_frequencies=state.initial_frequencies.to(device=device),
        display=display_value,
        display_interval=state.display_interval
        if display_interval is None
        else _positive_integer(display_interval, "display_interval"),
        output_fcn=output_fcn,
        time_limit=_nonnegative_scalar(time_limit, "time_limit", allow_inf=True),
        checkpoint_file=None if checkpoint_file is None else os.fspath(checkpoint_file),
        checkpoint_interval=_positive_integer(
            checkpoint_interval, "checkpoint_interval"
        ),
    )


def _configure_spatial_partition(context: _Context, local_size: int) -> int:
    if local_size < 1:
        raise ValueError("every distributed rank must own at least one spatial row")
    if context.world_size == 1:
        context.spatial_start = 0
        context.spatial_size = local_size
        return local_size
    local = torch.tensor([local_size], device=context.device, dtype=torch.int64)
    sizes = [torch.empty_like(local) for _ in range(context.world_size)]
    dist.all_gather(sizes, local)
    counts = [int(value.item()) for value in sizes]
    context.spatial_start = sum(counts[: context.rank])
    context.spatial_size = local_size
    return sum(counts)


def _validate_distributed_problem(
    context: _Context,
    sample_count: int,
    mode_count: int,
    alpha: float,
    is_real_input: bool,
) -> None:
    if context.world_size == 1:
        return
    values = torch.tensor(
        [sample_count, mode_count, int(is_real_input)],
        device=context.device,
        dtype=torch.int64,
    )
    gathered = [torch.empty_like(values) for _ in range(context.world_size)]
    dist.all_gather(gathered, values)
    if any(not torch.equal(values, other) for other in gathered):
        raise ValueError("all ranks must use the same T, K, and real/complex type")
    alpha_tensor = torch.tensor(alpha, device=context.device, dtype=torch.float64)
    alpha_values = [torch.empty_like(alpha_tensor) for _ in range(context.world_size)]
    dist.all_gather(alpha_values, alpha_tensor)
    if any(float(value.item()) != alpha for value in alpha_values):
        raise ValueError("all ranks must use the same alpha")


def _validate_distributed_settings(
    context: _Context,
    settings: _Settings,
) -> None:
    if context.world_size == 1:
        return
    values = torch.tensor(
        [
            settings.tolerance,
            settings.maximum_steps,
            settings.n_dc,
            settings.time_limit,
            settings.checkpoint_interval,
            int(settings.checkpoint_file is not None),
        ],
        device=context.device,
        dtype=torch.float64,
    )
    minimum, maximum = values.clone(), values.clone()
    dist.all_reduce(minimum, op=dist.ReduceOp.MIN)
    dist.all_reduce(maximum, op=dist.ReduceOp.MAX)
    mismatch = not torch.equal(minimum, maximum)

    checkpoint = (settings.checkpoint_file or "").encode("utf-8")
    expected_length = torch.tensor(
        len(checkpoint) if context.primary else 0,
        device=context.device,
        dtype=torch.int64,
    )
    context.broadcast_(expected_length)
    expected = torch.zeros(
        int(expected_length.item()), device=context.device, dtype=torch.uint8
    )
    if context.primary and checkpoint:
        expected.copy_(torch.tensor(list(checkpoint), device=context.device))
    context.broadcast_(expected)
    local = torch.tensor(list(checkpoint), device=context.device, dtype=torch.uint8)
    mismatch = (
        mismatch or local.shape != expected.shape or not torch.equal(local, expected)
    )
    mismatch_status = torch.tensor(
        int(mismatch), device=context.device, dtype=torch.int32
    )
    dist.all_reduce(mismatch_status, op=dist.ReduceOp.MAX)
    if mismatch_status.item():
        raise ValueError(
            "all ranks must use the same tolerance, maximum_steps, n_dc, "
            "time_limit, checkpoint interval, and checkpoint path"
        )


def _normalize_weight(
    weight: Any,
    local_size: int,
    global_size: int,
    real_dtype: torch.dtype,
    context: _Context,
) -> torch.Tensor:
    if weight is None:
        result = torch.ones(local_size, device=context.device, dtype=real_dtype)
    else:
        raw = _as_tensor(weight, context.device)
        if (
            torch.is_complex(raw)
            or raw.ndim > 2
            or (raw.ndim == 2 and raw.shape[0] > 1 and raw.shape[1] > 1)
        ):
            raise ValueError("weight must contain positive finite real values")
        result = raw.to(dtype=real_dtype).reshape(-1)
        if result.numel() == 1:
            result = result.expand(local_size).clone()
        elif result.numel() != local_size:
            raise ValueError("weight must be scalar or match the local spatial size")
    if not bool(torch.isfinite(result).all()) or bool((result <= 0).any()):
        raise ValueError("weight must contain positive finite values")
    total = result.sum()
    context.sum_(total)
    return result / (total / global_size)


def _initialize_frequencies(
    mode_count: int,
    n_dc: int,
    init_type: int,
    maximum: float,
    real_dtype: torch.dtype,
    context: _Context,
) -> torch.Tensor:
    result = torch.zeros(mode_count, device=context.device, dtype=real_dtype)
    free_modes = mode_count - n_dc
    if free_modes:
        if init_type == -1 and context.primary:
            result[n_dc:] = (
                torch.rand(free_modes, device=context.device, dtype=real_dtype)
                * maximum
            )
        elif init_type == 1:
            result[n_dc:] = (
                torch.arange(1, free_modes + 1, device=context.device, dtype=real_dtype)
                / free_modes
                * maximum
            )
    context.broadcast_(result)
    return result


def _initial_phi(
    local_size: int,
    global_size: int,
    mode_count: int,
    weight: torch.Tensor,
    is_real_input: bool,
    real_dtype: torch.dtype,
    complex_dtype: torch.dtype,
    context: _Context,
) -> torch.Tensor:
    dtype = real_dtype if is_real_input else complex_dtype
    phi = torch.zeros((local_size, mode_count), device=context.device, dtype=dtype)
    for index in range(mode_count):
        global_index = index % global_size
        if context.spatial_start <= global_index < context.spatial_stop:
            local_index = global_index - context.spatial_start
            phi[local_index, index] = 1 / torch.sqrt(weight[local_index])
    return phi


def _global_norm(tensor: torch.Tensor, context: _Context) -> torch.Tensor:
    squared = torch.sum(torch.abs(tensor) ** 2)
    context.sum_(squared)
    return torch.sqrt(squared)


def _canonicalize_phase(
    phi: torch.Tensor, is_real_input: bool, context: _Context
) -> torch.Tensor:
    local_magnitude, local_index = torch.max(torch.abs(phi), dim=0)
    if context.world_size == 1:
        pivot_value = phi[local_index, torch.arange(phi.shape[1], device=phi.device)]
        if is_real_input:
            signs = torch.where(pivot_value < 0, -1, 1).to(phi.dtype)
            return phi * signs
        pivot_abs = torch.abs(pivot_value)
        valid = torch.isfinite(pivot_abs) & (pivot_abs > 0)
        factors = torch.ones_like(pivot_value)
        factors[valid] = torch.conj(pivot_value[valid]) / pivot_abs[valid]
        phi = phi * factors
        phi[
            local_index[valid], torch.arange(phi.shape[1], device=phi.device)[valid]
        ] = pivot_abs[valid].to(phi.dtype)
        return phi

    gathered_magnitude = [
        torch.empty_like(local_magnitude) for _ in range(context.world_size)
    ]
    dist.all_gather(gathered_magnitude, local_magnitude)
    magnitude_by_rank = torch.stack(gathered_magnitude, dim=0)
    owner = torch.argmax(magnitude_by_rank, dim=0)
    for mode_index in range(phi.shape[1]):
        owner_rank = int(owner[mode_index].item())
        if context.rank == owner_rank:
            pivot = phi[local_index[mode_index], mode_index].reshape(1)
        else:
            pivot = torch.zeros(1, device=phi.device, dtype=phi.dtype)
        context.broadcast_(pivot, owner_rank)
        if is_real_input:
            if bool((pivot < 0).item()):
                phi[:, mode_index] = -phi[:, mode_index]
        else:
            pivot_abs = torch.abs(pivot)
            if bool((torch.isfinite(pivot_abs) & (pivot_abs > 0)).item()):
                phi[:, mode_index] *= torch.conj(pivot[0]) / pivot_abs[0]
                if context.rank == owner_rank:
                    phi[local_index[mode_index], mode_index] = pivot_abs[0]
    return phi


def _fallback_phi(
    phi: torch.Tensor,
    mode_index: int,
    weight: torch.Tensor,
    is_real_input: bool,
    global_size: int,
    context: _Context,
) -> torch.Tensor:
    norm_squared = torch.sum(torch.abs(phi) ** 2 * weight)
    context.sum_(norm_squared)
    norm = torch.sqrt(norm_squared)
    if bool((torch.isfinite(norm) & (norm > 0)).item()):
        result = phi / norm
    else:
        result = torch.zeros_like(phi)
        global_index = mode_index % global_size
        if context.spatial_start <= global_index < context.spatial_stop:
            local_index = global_index - context.spatial_start
            result[local_index] = 1 / torch.sqrt(weight[local_index])
    return result.real if is_real_input else result


def _progress(
    step: int,
    difference: float,
    omega: torch.Tensor,
    elapsed: float,
    stop_reason: str = "",
) -> Progress:
    return Progress(
        step=step,
        difference=float(difference),
        omega=omega.detach().cpu().double().numpy().copy(),
        elapsed_time=float(elapsed),
        stop_reason=stop_reason,
    )


def _invoke_output_function(
    function: OutputFunction | None,
    progress: Progress,
    phase: Literal["init", "iter", "done"],
    context: _Context,
) -> bool:
    if function is None:
        return False
    callback_error: BaseException | None = None
    stop = False
    if context.primary:
        try:
            result = function(progress, phase)
            if result is not None:
                if not np.isscalar(result) or not np.isfinite(result):
                    raise ValueError("output_fcn must return None or a finite scalar")
                stop = bool(result)
        except BaseException as exception:  # synchronize ranks before propagating
            callback_error = exception
    status = torch.tensor(
        [int(stop), int(callback_error is not None)],
        device=context.device,
        dtype=torch.int32,
    )
    context.broadcast_(status)
    if int(status[1].item()):
        if callback_error is not None:
            raise callback_error
        raise RuntimeError("output_fcn failed on rank 0")
    return bool(status[0].item())


def _display_iteration(
    step: int, difference: float, omega: torch.Tensor, header: bool
) -> None:
    if header:
        print(" Iteration    Difference      Center frequencies (cycles/sample)")
    values = np.sort(omega.detach().cpu().double().numpy())
    print(
        f"{step:10d}    {difference:10.3e}      {np.array2string(values, precision=6)}"
    )


def _termination(
    completed_steps: int,
    difference: float,
    settings: _Settings,
    output_stop: bool,
    time_stop: bool,
    elapsed: float,
) -> tuple[int, StopReason, str]:
    if completed_steps > 0 and difference <= settings.tolerance:
        exit_flag, reason, description = 1, "converged", "converged"
    elif output_stop:
        exit_flag, reason, description = -1, "output_function", "stopped by output_fcn"
    elif time_stop:
        exit_flag, reason, description = -2, "time_limit", "reached time_limit"
    else:
        exit_flag, reason, description = 0, "maximum_steps", "reached maximum_steps"
    return (
        exit_flag,
        reason,
        f"RVMD {description} after {completed_steps} iterations ({elapsed:.3g} s).",
    )


def _make_state(
    settings: _Settings,
    context: _Context,
    global_size: int,
    sample_count: int,
    mode_count: int,
    alpha: float,
    weight: torch.Tensor,
    is_real_input: bool,
    initial_frequencies: torch.Tensor,
    completed_steps: int,
    frequency_history: torch.Tensor,
    difference_history: torch.Tensor,
    coefficient_spectrum: torch.Tensor,
    phi: torch.Tensor,
    residual: torch.Tensor,
    data_spectrum_norm: float,
) -> RestartState:
    return RestartState(
        version=PYTHON_RESTART_VERSION,
        spatial_size=global_size,
        sample_count=sample_count,
        mode_count=mode_count,
        alpha=float(alpha),
        weight=weight.detach().cpu().clone(),
        tolerance=settings.tolerance,
        maximum_steps=settings.maximum_steps,
        init_freq_type=settings.init_freq_type,
        init_freq_maximum=settings.init_freq_maximum,
        initial_frequencies=initial_frequencies.detach().cpu().clone(),
        device=settings.device.type,
        fp_precision=settings.fp_precision,
        n_dc=settings.n_dc,
        is_real_input=is_real_input,
        display=settings.display,
        display_interval=settings.display_interval,
        steps=completed_steps,
        omega_history=frequency_history[:, : completed_steps + 1]
        .detach()
        .cpu()
        .clone(),
        difference_history=difference_history[:completed_steps].detach().cpu().clone(),
        coefficient_spectrum=coefficient_spectrum.detach().cpu().clone(),
        phi=phi.detach().cpu().clone(),
        residual=residual.detach().cpu().clone(),
        data_spectrum_norm=float(data_spectrum_norm),
        spatial_start=context.spatial_start,
        world_size=context.world_size,
    )


def _save_state(path: str, state: RestartState, context: _Context) -> None:
    if context.distributed:
        save_distributed_checkpoint(path, state)
    elif context.primary:
        save_checkpoint(path, state)


def _load_restart(
    restart: RestartState | str | os.PathLike[str], context: _Context
) -> RestartState:
    if isinstance(restart, RestartState):
        state = restart
    elif context.distributed:
        state = load_distributed_checkpoint(restart)
    else:
        state = load_checkpoint(restart)
    if state.version != PYTHON_RESTART_VERSION:
        raise ValueError(f"unsupported Python restart version {state.version}")
    if context.distributed:
        if state.world_size != context.world_size:
            raise ValueError(
                "in-memory distributed restart world_size differs; reload its "
                "checkpoint path to reshard"
            )
        context.spatial_start = state.spatial_start
        context.spatial_size = state.phi.shape[0]
    elif state.spatial_start != 0 or state.phi.shape[0] != state.spatial_size:
        raise ValueError("load a distributed checkpoint path to assemble a full state")
    return state


@torch.no_grad()
def _solve(
    q: Any | None,
    k: int | None,
    alpha: float | None,
    *,
    restart: RestartState | str | os.PathLike[str] | None,
    weight: Any,
    tolerance: float | None,
    maximum_steps: int | None,
    init_freq_type: int | None,
    init_freq_maximum: float | None,
    device: str | torch.device | None,
    fp_precision: str | None,
    n_dc: int | None,
    initial_frequencies: Any,
    display: str | None,
    display_interval: int | None,
    output_fcn: OutputFunction | None,
    time_limit: float,
    checkpoint_file: str | os.PathLike[str] | None,
    checkpoint_interval: int,
    distributed: bool,
) -> tuple[Mode, Info, RestartState]:
    call_start = time.perf_counter()
    context = _initialize_context(device, distributed)
    restarting = restart is not None

    if restarting:
        if q is not None or k is not None or alpha is not None:
            raise ValueError("q, k, and alpha must be omitted when restart is supplied")
        if any(
            value is not None
            for value in (
                weight,
                init_freq_type,
                init_freq_maximum,
                fp_precision,
                n_dc,
                initial_frequencies,
            )
        ):
            raise ValueError(
                "initialization, precision, n_dc, and weight are immutable"
            )
        state = _load_restart(restart, context)
        if device is None and not distributed:
            desired = _resolve_device(state.device, distributed)
            if desired != context.device:
                context = _initialize_context(desired, distributed)
                if distributed:
                    context.spatial_start = state.spatial_start
                    context.spatial_size = state.phi.shape[0]
        settings = _restart_settings(
            state,
            tolerance=tolerance,
            maximum_steps=maximum_steps,
            device=context.device,
            display=display,
            display_interval=display_interval,
            output_fcn=output_fcn,
            time_limit=time_limit,
            checkpoint_file=checkpoint_file,
            checkpoint_interval=checkpoint_interval,
        )
        global_size = state.spatial_size
        sample_count = state.sample_count
        mode_count = state.mode_count
        alpha_value = state.alpha
        is_real_input = state.is_real_input
        real_dtype, complex_dtype = _real_complex_dtype(state.fp_precision)
        local_weight = state.weight.to(context.device, dtype=real_dtype)
        phi = state.phi.to(
            context.device, dtype=real_dtype if is_real_input else complex_dtype
        )
        coefficient_spectrum = state.coefficient_spectrum.to(
            context.device, dtype=complex_dtype
        )
        residual = state.residual.to(context.device, dtype=complex_dtype)
        completed_steps = state.steps
        initial = state.initial_frequencies.to(context.device, dtype=real_dtype)
        stored_omega = state.omega_history.to(context.device, dtype=real_dtype)
        stored_difference = state.difference_history.to(
            context.device, dtype=real_dtype
        )
        data_spectrum_norm = state.data_spectrum_norm
    else:
        if q is None or k is None or alpha is None:
            raise ValueError("q, k, and alpha are required for a new decomposition")
        mode_count = _positive_integer(k, "k")
        alpha_value = _nonnegative_scalar(alpha, "alpha")
        provisional_precision = (
            "single" if fp_precision is None else fp_precision.lower()
        )
        real_dtype, complex_dtype = _real_complex_dtype(provisional_precision)
        settings = _new_settings(
            mode_count=mode_count,
            tolerance=tolerance,
            maximum_steps=maximum_steps,
            init_freq_type=init_freq_type,
            init_freq_maximum=init_freq_maximum,
            device=context.device,
            fp_precision=fp_precision,
            n_dc=n_dc,
            initial_frequencies=initial_frequencies,
            display=display,
            display_interval=display_interval,
            output_fcn=output_fcn,
            time_limit=time_limit,
            checkpoint_file=checkpoint_file,
            checkpoint_interval=checkpoint_interval,
        )
        q_tensor = _as_tensor(q, context.device)
        if q_tensor.ndim != 2 or q_tensor.numel() == 0:
            raise ValueError("q must be a nonempty S-by-T matrix")
        if not bool(torch.isfinite(q_tensor).all()):
            raise ValueError("q must contain finite values")
        is_real_input = not torch.is_complex(q_tensor)
        q_tensor = q_tensor.to(
            dtype=real_dtype if is_real_input else complex_dtype
        ).contiguous()
        local_size, sample_count = q_tensor.shape
        global_size = _configure_spatial_partition(context, local_size)
        _validate_distributed_problem(
            context, sample_count, mode_count, alpha_value, is_real_input
        )
        local_weight = _normalize_weight(
            weight, local_size, global_size, real_dtype, context
        )
        if settings.initial_frequencies is None:
            initial = _initialize_frequencies(
                mode_count,
                settings.n_dc,
                settings.init_freq_type,
                settings.init_freq_maximum,
                real_dtype,
                context,
            )
        else:
            initial = settings.initial_frequencies.to(context.device)
            context.broadcast_(initial)
        completed_steps = 0
        stored_omega = torch.empty(
            (mode_count, 0), device=context.device, dtype=real_dtype
        )
        stored_difference = torch.empty(0, device=context.device, dtype=real_dtype)

        extended_length = 2 * sample_count
        half = math.ceil(sample_count / 2)
        extended = torch.cat(
            (
                torch.flip(q_tensor[:, :half], dims=(1,)),
                q_tensor,
                torch.flip(q_tensor[:, half:], dims=(1,)),
            ),
            dim=1,
        )
        q_spectrum = torch.fft.fft(extended, dim=1)
        if is_real_input:
            q_spectrum = q_spectrum[:, : sample_count + 1]
        else:
            q_spectrum = torch.fft.fftshift(q_spectrum, dim=1)
        data_spectrum_norm = float(_global_norm(q_spectrum, context).item())
        phi = _initial_phi(
            local_size,
            global_size,
            mode_count,
            local_weight,
            is_real_input,
            real_dtype,
            complex_dtype,
            context,
        )
        spectrum_length = q_spectrum.shape[1]
        if data_spectrum_norm == 0:
            coefficient_spectrum = torch.zeros(
                (spectrum_length, mode_count),
                device=context.device,
                dtype=complex_dtype,
            )
        else:
            coefficient_spectrum = torch.full(
                (spectrum_length, mode_count),
                torch.finfo(real_dtype).eps,
                device=context.device,
                dtype=complex_dtype,
            )
        residual = q_spectrum - phi.to(complex_dtype) @ coefficient_spectrum.transpose(
            0, 1
        )
        del q_spectrum, extended, q_tensor

    _validate_distributed_settings(context, settings)
    local_size = phi.shape[0]
    if context.spatial_size == 0:
        context.spatial_size = local_size
    extended_length = 2 * sample_count
    half = math.ceil(sample_count / 2)
    if is_real_input:
        spectrum_length = sample_count + 1
        frequency = (
            torch.arange(sample_count + 1, device=context.device, dtype=real_dtype)
            / extended_length
        )
        frequency_weight = torch.ones(
            spectrum_length, device=context.device, dtype=real_dtype
        )
        frequency_weight[[0, -1]] = 0.5
    else:
        spectrum_length = extended_length
        frequency = (
            torch.arange(
                -sample_count, sample_count, device=context.device, dtype=real_dtype
            )
            / extended_length
        )
        frequency_weight = torch.ones(
            spectrum_length, device=context.device, dtype=real_dtype
        )
    absolute_frequency = torch.abs(frequency)

    expected_phi_dtype = real_dtype if is_real_input else complex_dtype
    if phi.shape != (local_size, mode_count) or phi.dtype != expected_phi_dtype:
        raise ValueError("restart phi has inconsistent shape or dtype")
    if coefficient_spectrum.shape != (spectrum_length, mode_count):
        raise ValueError("restart coefficient spectrum has inconsistent dimensions")
    if residual.shape != (local_size, spectrum_length):
        raise ValueError("restart residual has inconsistent dimensions")

    history_capacity = max(settings.maximum_steps + 1, completed_steps + 1)
    frequency_history = torch.zeros(
        (mode_count, history_capacity), device=context.device, dtype=real_dtype
    )
    difference_history = torch.zeros(
        max(settings.maximum_steps, completed_steps),
        device=context.device,
        dtype=real_dtype,
    )
    if restarting:
        frequency_history[:, : completed_steps + 1] = stored_omega[
            :, : completed_steps + 1
        ]
        if completed_steps:
            difference_history[:completed_steps] = stored_difference[:completed_steps]
    else:
        frequency_history[:, 0] = initial

    scale_floor = torch.tensor(
        torch.finfo(real_dtype).eps * max(data_spectrum_norm, 1.0),
        device=context.device,
        dtype=real_dtype,
    )
    difference = (
        math.inf
        if completed_steps == 0
        else float(difference_history[completed_steps - 1].item())
    )
    iteration = completed_steps
    output_stop = _invoke_output_function(
        settings.output_fcn,
        _progress(
            iteration,
            difference,
            frequency_history[:, iteration],
            time.perf_counter() - call_start,
        ),
        "init",
        context,
    )
    time_stop = False
    last_displayed: int | None = None
    last_checkpoint: int | None = None
    mode_buffer = torch.empty(
        (local_size, spectrum_length), device=context.device, dtype=complex_dtype
    )

    while (
        iteration < settings.maximum_steps
        and difference > settings.tolerance
        and not output_stop
    ):
        elapsed = time.perf_counter() - call_start
        if context.world_size > 1:
            time_status = torch.tensor(
                int(context.primary and elapsed >= settings.time_limit),
                device=context.device,
                dtype=torch.int32,
            )
            context.broadcast_(time_status)
            time_stop = bool(time_status.item())
        else:
            time_stop = elapsed >= settings.time_limit
        if time_stop:
            break

        difference_accumulator = torch.zeros(
            (), device=context.device, dtype=real_dtype
        )
        current_frequency = frequency_history[:, iteration]
        next_frequency = frequency_history[:, iteration + 1]

        for mode_index in range(mode_count):
            old_phi = phi[:, mode_index].clone()
            old_coefficient = coefficient_spectrum[:, mode_index].clone()
            torch.mul(
                old_phi[:, None].to(complex_dtype),
                old_coefficient[None, :],
                out=mode_buffer,
            )
            residual.add_(mode_buffer)

            projection = residual @ (
                torch.conj(coefficient_spectrum[:, mode_index]) * frequency_weight
            )
            if is_real_input:
                projection = projection.real
            projection_norm = torch.sum(torch.abs(projection) ** 2 * local_weight)
            context.sum_(projection_norm)
            projection_norm = torch.sqrt(projection_norm)
            if bool((torch.isfinite(projection_norm) & (projection_norm > 0)).item()):
                phi[:, mode_index] = projection / projection_norm
            else:
                phi[:, mode_index] = _fallback_phi(
                    phi[:, mode_index],
                    mode_index,
                    local_weight,
                    is_real_input,
                    global_size,
                    context,
                )
            phi[:, mode_index : mode_index + 1] = _canonicalize_phase(
                phi[:, mode_index : mode_index + 1], is_real_input, context
            )

            denominator = (
                1
                + 2
                * alpha_value
                * (absolute_frequency - current_frequency[mode_index]) ** 2
            )
            numerator = residual.transpose(0, 1) @ (
                torch.conj(phi[:, mode_index]).to(complex_dtype) * local_weight
            )
            context.sum_(numerator)
            coefficient_spectrum[:, mode_index] = numerator / denominator

            coefficient_energy = torch.sum(
                frequency_weight * torch.abs(coefficient_spectrum[:, mode_index]) ** 2
            )
            if mode_index < settings.n_dc:
                next_frequency[mode_index] = 0
            elif bool(
                (torch.isfinite(coefficient_energy) & (coefficient_energy > 0)).item()
            ):
                next_frequency[mode_index] = (
                    torch.sum(
                        frequency_weight
                        * absolute_frequency
                        * torch.abs(coefficient_spectrum[:, mode_index]) ** 2
                    )
                    / coefficient_energy
                )
            else:
                next_frequency[mode_index] = current_frequency[mode_index]

            new_phi = phi[:, mode_index]
            new_coefficient = coefficient_spectrum[:, mode_index]
            torch.mul(
                new_phi[:, None].to(complex_dtype),
                new_coefficient[None, :],
                out=mode_buffer,
            )
            residual.sub_(mode_buffer)

            spatial_products = torch.stack(
                (
                    torch.sum(torch.abs(old_phi) ** 2).to(complex_dtype),
                    torch.sum(torch.abs(new_phi) ** 2).to(complex_dtype),
                    torch.sum(torch.conj(old_phi) * new_phi).to(complex_dtype),
                )
            )
            context.sum_(spatial_products)
            old_coefficient_norm = torch.sum(torch.abs(old_coefficient) ** 2)
            new_coefficient_norm = torch.sum(torch.abs(new_coefficient) ** 2)
            coefficient_cross = torch.sum(torch.conj(old_coefficient) * new_coefficient)
            old_norm_squared = spatial_products[0].real * old_coefficient_norm
            new_norm_squared = spatial_products[1].real * new_coefficient_norm
            cross = spatial_products[2] * coefficient_cross
            change_norm_squared = torch.clamp(
                old_norm_squared + new_norm_squared - 2 * cross.real, min=0
            )
            old_norm = torch.sqrt(old_norm_squared)
            change_norm = torch.sqrt(change_norm_squared)
            difference_accumulator.add_(
                change_norm / torch.maximum(old_norm, scale_floor)
            )

        iteration += 1
        difference_history[iteration - 1] = difference_accumulator
        difference = float(difference_accumulator.item())

        if (
            context.primary
            and settings.display == "iter"
            and (
                iteration == completed_steps + 1
                or iteration % settings.display_interval == 0
            )
        ):
            _display_iteration(
                iteration,
                difference,
                frequency_history[:, iteration],
                last_displayed is None,
            )
            last_displayed = iteration

        output_stop = _invoke_output_function(
            settings.output_fcn,
            _progress(
                iteration,
                difference,
                frequency_history[:, iteration],
                time.perf_counter() - call_start,
            ),
            "iter",
            context,
        )

        if (
            settings.checkpoint_file is not None
            and iteration % settings.checkpoint_interval == 0
        ):
            checkpoint_state = _make_state(
                settings,
                context,
                global_size,
                sample_count,
                mode_count,
                alpha_value,
                local_weight,
                is_real_input,
                initial,
                iteration,
                frequency_history,
                difference_history,
                coefficient_spectrum,
                phi,
                residual,
                data_spectrum_norm,
            )
            _save_state(settings.checkpoint_file, checkpoint_state, context)
            last_checkpoint = iteration

    completed_steps = iteration
    if (
        context.primary
        and settings.display == "iter"
        and completed_steps > 0
        and last_displayed != completed_steps
    ):
        _display_iteration(
            completed_steps,
            difference,
            frequency_history[:, completed_steps],
            last_displayed is None,
        )

    termination_elapsed = time.perf_counter() - call_start
    exit_flag, stop_reason, message = _termination(
        completed_steps,
        difference,
        settings,
        output_stop,
        time_stop,
        termination_elapsed,
    )
    restart_state = _make_state(
        settings,
        context,
        global_size,
        sample_count,
        mode_count,
        alpha_value,
        local_weight,
        is_real_input,
        initial,
        completed_steps,
        frequency_history,
        difference_history,
        coefficient_spectrum,
        phi,
        residual,
        data_spectrum_norm,
    )
    if settings.checkpoint_file is not None and last_checkpoint != completed_steps:
        _save_state(settings.checkpoint_file, restart_state, context)
        last_checkpoint = completed_steps

    if is_real_input:
        full_spectrum = torch.zeros(
            (extended_length, mode_count), device=context.device, dtype=complex_dtype
        )
        full_spectrum[:spectrum_length] = coefficient_spectrum
        full_spectrum[spectrum_length:] = torch.conj(
            torch.flip(coefficient_spectrum[1:sample_count], dims=(0,))
        )
        coefficient = torch.fft.ifft(full_spectrum, dim=0).real
        output_phi = phi.real
    else:
        coefficient = torch.fft.ifft(
            torch.fft.ifftshift(coefficient_spectrum, dim=0), dim=0
        )
        output_phi = phi
    coefficient = coefficient[half : half + sample_count]
    final_frequency = frequency_history[:, completed_steps]
    order = torch.argsort(final_frequency)
    mode = Mode(
        phi=output_phi[:, order].detach().cpu().clone(),
        c=coefficient[:, order].detach().cpu().clone(),
        omega=final_frequency[order].detach().cpu().clone(),
        energy=torch.sum(torch.abs(coefficient[:, order]) ** 2, dim=0)
        .detach()
        .cpu()
        .clone(),
        spatial_start=context.spatial_start,
        global_spatial_size=global_size,
    )
    final_elapsed = time.perf_counter() - call_start
    info = Info(
        spatial_size=global_size,
        sample_count=sample_count,
        mode_count=mode_count,
        alpha=alpha_value,
        weight=local_weight.detach().cpu().clone(),
        tolerance=settings.tolerance,
        maximum_steps=settings.maximum_steps,
        init_freq_type=settings.init_freq_type,
        init_freq_maximum=settings.init_freq_maximum,
        device=str(context.device),
        fp_precision=settings.fp_precision,
        n_dc=settings.n_dc,
        is_real_input=is_real_input,
        initial_frequencies=initial.detach().cpu().clone(),
        display_interval=settings.display_interval,
        time_limit=settings.time_limit,
        elapsed_time=final_elapsed,
        exit_flag=exit_flag,
        stop_reason=stop_reason,
        message=message,
        checkpoint_file=settings.checkpoint_file,
        last_checkpoint_step=last_checkpoint,
        iteration=IterationInfo(
            steps=completed_steps,
            omega=frequency_history[:, : completed_steps + 1].detach().cpu().clone(),
            difference=difference_history[:completed_steps].detach().cpu().clone(),
            converged=exit_flag == 1,
        ),
        spatial_start=context.spatial_start,
    )

    if context.primary and settings.display != "off":
        print(message)
    _invoke_output_function(
        settings.output_fcn,
        _progress(
            completed_steps,
            difference,
            frequency_history[:, completed_steps],
            final_elapsed,
            stop_reason,
        ),
        "done",
        context,
    )
    return mode, info, restart_state


def rvmd(
    q: Any | None = None,
    k: int | None = None,
    alpha: float | None = None,
    *,
    restart: RestartState | str | os.PathLike[str] | None = None,
    weight: Any = None,
    tolerance: float | None = None,
    maximum_steps: int | None = None,
    init_freq_type: int | None = None,
    init_freq_maximum: float | None = None,
    device: str | torch.device | None = None,
    fp_precision: str | None = None,
    n_dc: int | None = None,
    initial_frequencies: Any = None,
    display: str | None = None,
    display_interval: int | None = None,
    output_fcn: OutputFunction | None = None,
    time_limit: float = math.inf,
    checkpoint_file: str | os.PathLike[str] | None = None,
    checkpoint_interval: int = 50,
) -> tuple[Mode, Info, RestartState]:
    """Decompose one full ``S``-by-``T`` matrix on CPU or one CUDA GPU.

    Args:
        q: Finite real or complex matrix. Omit on restart.
        k: Positive mode count. Omit on restart.
        alpha: Nonnegative bandwidth penalty. Omit on restart.
        restart: A :class:`RestartState` or Python checkpoint path. Supplying it
            makes initialization, precision, weights, ``k``, and ``alpha``
            immutable.
        weight: Positive scalar or length-``S`` vector, normalized to mean one.
            The new-run default is one.
        tolerance: Nonnegative convergence tolerance. Defaults to ``5e-3`` for
            a new run and to the saved value on restart.
        maximum_steps: Positive total sweep limit, including restart steps.
            Defaults to 500 for a new run and to the saved value on restart.
        init_freq_type: ``-1`` random, ``0`` zero, or ``1`` uniform. Default 1.
        init_freq_maximum: Initial upper frequency, capped at 0.5. Default 0.5.
        device: ``"cpu"`` or a CUDA device. A new run defaults to CPU; a
            single-device restart defaults to its saved device type.
        fp_precision: ``"single"`` or ``"double"``. Default ``"single"``.
        n_dc: Number of leading internal modes fixed at zero frequency.
            Default 0.
        initial_frequencies: Optional finite length-``K`` values in ``[0,0.5]``.
        display: ``"off"``, ``"final"``, or ``"iter"``. Default ``"off"``.
        display_interval: Positive progress-print interval. Default 20.
        output_fcn: ``callback(progress, phase)`` for phases ``init``, ``iter``,
            and ``done``. A nonzero result at ``init`` or ``iter`` stops safely.
        time_limit: Per-call seconds, checked before each complete sweep.
        checkpoint_file: Atomic Python checkpoint path, or ``None`` to disable.
        checkpoint_interval: Positive completed-sweep interval. Default 50.

    Returns:
        ``(mode, info, restart_state)``. All public result tensors are on CPU.
        The reconstruction is ``mode.phi @ mode.c.T`` without conjugation.
    """

    return _solve(
        q,
        k,
        alpha,
        restart=restart,
        weight=weight,
        tolerance=tolerance,
        maximum_steps=maximum_steps,
        init_freq_type=init_freq_type,
        init_freq_maximum=init_freq_maximum,
        device=device,
        fp_precision=fp_precision,
        n_dc=n_dc,
        initial_frequencies=initial_frequencies,
        display=display,
        display_interval=display_interval,
        output_fcn=output_fcn,
        time_limit=time_limit,
        checkpoint_file=checkpoint_file,
        checkpoint_interval=checkpoint_interval,
        distributed=False,
    )


def rvmd_distributed(
    local_q: Any | None = None,
    k: int | None = None,
    alpha: float | None = None,
    *,
    restart: RestartState | str | os.PathLike[str] | None = None,
    local_weight: Any = None,
    tolerance: float | None = None,
    maximum_steps: int | None = None,
    init_freq_type: int | None = None,
    init_freq_maximum: float | None = None,
    device: str | torch.device | None = None,
    fp_precision: str | None = None,
    n_dc: int | None = None,
    initial_frequencies: Any = None,
    display: str | None = None,
    display_interval: int | None = None,
    output_fcn: OutputFunction | None = None,
    time_limit: float = math.inf,
    checkpoint_file: str | os.PathLike[str] | None = None,
    checkpoint_interval: int = 50,
) -> tuple[Mode, Info, RestartState]:
    """Decompose one matrix sharded by spatial rows across torchrun ranks.

    Every rank supplies contiguous ``local_q`` rows and an optional scalar or
    matching ``local_weight``. Rank order defines the global row order. CUDA
    runs use one process per GPU and NCCL; CPU process groups use Gloo. The
    numerical options and defaults match :func:`rvmd`. ``device=None`` selects
    the rank-local CUDA device when CUDA is available and otherwise selects
    CPU, including during restart. All ranks must use identical numerical and
    checkpoint settings; only their local input and weights differ.

    A checkpoint path names a directory. Each rank writes its local spatial
    state while rank 0 writes replicated state. Loading the path can reshard to
    a different world size. An in-memory restart requires the same world size.

    Args:
        local_q: This rank's finite, contiguous spatial rows of the common
            ``S``-by-``T`` matrix. Omit on restart.
        k: Positive shared mode count. Omit on restart.
        alpha: Nonnegative shared bandwidth penalty. Omit on restart.
        restart: A local :class:`RestartState` shard or distributed checkpoint
            directory. Supplying a directory permits resharding.
        local_weight: Positive scalar or vector matching this rank's rows.
            Global normalization is performed collectively.
        tolerance: Shared convergence tolerance; see :func:`rvmd`.
        maximum_steps: Shared total sweep limit; see :func:`rvmd`.
        init_freq_type: Shared frequency initialization type; see :func:`rvmd`.
        init_freq_maximum: Shared initial-frequency upper bound; see
            :func:`rvmd`.
        device: Rank-local CPU or CUDA device. ``None`` chooses local CUDA when
            available, otherwise CPU.
        fp_precision: Shared ``"single"`` or ``"double"`` precision.
        n_dc: Shared number of leading zero-frequency modes.
        initial_frequencies: Shared explicit initial frequency vector.
        display: Rank-0 display mode; see :func:`rvmd`.
        display_interval: Shared completed-sweep display interval.
        output_fcn: Rank-0 callback; its stop result is broadcast to all ranks.
        time_limit: Shared per-call limit checked by rank 0 at sweep boundaries.
        checkpoint_file: Shared checkpoint directory, or ``None`` to disable.
        checkpoint_interval: Shared completed-sweep checkpoint interval.

    Returns:
        ``(mode, info, restart_state)`` on every rank. ``mode.phi``,
        ``info.weight``, and the spatial restart arrays are local shards;
        ``mode.c``, ``mode.omega``, and histories are replicated CPU tensors.
    """

    return _solve(
        local_q,
        k,
        alpha,
        restart=restart,
        weight=local_weight,
        tolerance=tolerance,
        maximum_steps=maximum_steps,
        init_freq_type=init_freq_type,
        init_freq_maximum=init_freq_maximum,
        device=device,
        fp_precision=fp_precision,
        n_dc=n_dc,
        initial_frequencies=initial_frequencies,
        display=display,
        display_interval=display_interval,
        output_fcn=output_fcn,
        time_limit=time_limit,
        checkpoint_file=checkpoint_file,
        checkpoint_interval=checkpoint_interval,
        distributed=True,
    )


def gather_mode(mode: Mode, destination: int = 0) -> Mode | None:
    """Gather spatial modes to one rank without replicating them everywhere.

    All ranks must call this collective. The destination receives a full
    :class:`Mode`; other ranks receive ``None``. Do not call it when the full
    spatial mode matrix cannot fit in destination memory.

    Args:
        mode: The local :class:`Mode` returned by :func:`rvmd_distributed`.
        destination: Rank that receives the concatenated spatial rows.

    Returns:
        A full :class:`Mode` on ``destination``, ``None`` on other ranks, or
        the input unchanged when no multi-rank process group exists.
    """

    if (
        not dist.is_available()
        or not dist.is_initialized()
        or dist.get_world_size() == 1
    ):
        return mode
    rank, world_size = dist.get_rank(), dist.get_world_size()
    backend = dist.get_backend()
    if backend == "nccl":
        device = torch.device("cuda", int(os.environ.get("LOCAL_RANK", "0")))
    else:
        device = torch.device("cpu")
    destination_tensor = torch.tensor(destination, device=device, dtype=torch.int64)
    destination_minimum = destination_tensor.clone()
    destination_maximum = destination_tensor.clone()
    dist.all_reduce(destination_minimum, op=dist.ReduceOp.MIN)
    dist.all_reduce(destination_maximum, op=dist.ReduceOp.MAX)
    if (
        destination_minimum.item() != destination_maximum.item()
        or destination < 0
        or destination >= world_size
    ):
        raise ValueError("all ranks must specify the same valid destination")
    local_phi = mode.phi.to(device)
    local_rows = torch.tensor([local_phi.shape[0]], device=device, dtype=torch.int64)
    row_counts = [torch.empty_like(local_rows) for _ in range(world_size)]
    dist.all_gather(row_counts, local_rows)
    counts = [int(value.item()) for value in row_counts]
    maximum_rows = max(counts)
    padded = torch.zeros(
        (maximum_rows, local_phi.shape[1]), device=device, dtype=local_phi.dtype
    )
    padded[: local_phi.shape[0]] = local_phi
    gathered = (
        [torch.empty_like(padded) for _ in range(world_size)]
        if rank == destination
        else None
    )
    dist.gather(padded, gathered, dst=destination)
    if rank != destination:
        return None
    assert gathered is not None
    full_phi = torch.cat(
        [part[:count].cpu() for part, count in zip(gathered, counts, strict=True)],
        dim=0,
    )
    return Mode(
        phi=full_phi,
        c=mode.c,
        omega=mode.omega,
        energy=mode.energy,
        spatial_start=0,
        global_spatial_size=full_phi.shape[0],
    )
