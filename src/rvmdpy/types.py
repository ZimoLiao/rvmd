"""Public result and restart-state types for :mod:`rvmdpy`."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
import scipy.sparse
import torch

StopReason = Literal["converged", "maximum_steps", "output_function", "time_limit"]


@dataclass(slots=True)
class Mode:
    """RVMD modes.

    ``phi`` contains either all spatial rows or the local spatial shard.
    ``c`` is the full time-coefficient matrix on every rank.
    """

    phi: torch.Tensor
    c: torch.Tensor
    omega: torch.Tensor
    energy: torch.Tensor
    spatial_start: int = 0
    global_spatial_size: int | None = None

    @property
    def spatial_stop(self) -> int:
        return self.spatial_start + self.phi.shape[0]

    @property
    def is_distributed(self) -> bool:
        return self.global_spatial_size is not None and (
            self.spatial_start != 0 or self.phi.shape[0] != self.global_spatial_size
        )


@dataclass(slots=True)
class IterationInfo:
    """Unsorted internal iteration history."""

    steps: int
    omega: torch.Tensor
    difference: torch.Tensor
    converged: bool


@dataclass(slots=True)
class Info:
    """Effective settings and termination information for an RVMD call."""

    spatial_size: int
    sample_count: int
    mode_count: int
    alpha: float
    weight: torch.Tensor
    tolerance: float
    maximum_steps: int
    init_freq_type: int
    init_freq_maximum: float
    device: str
    fp_precision: str
    n_dc: int
    is_real_input: bool
    initial_frequencies: torch.Tensor
    display_interval: int
    time_limit: float
    elapsed_time: float
    exit_flag: int
    stop_reason: StopReason
    message: str
    checkpoint_file: str | None
    last_checkpoint_step: int | None
    iteration: IterationInfo
    spatial_start: int = 0


@dataclass(slots=True)
class RestartState:
    """Opaque exact continuation state returned by RVMD."""

    version: int
    spatial_size: int
    sample_count: int
    mode_count: int
    alpha: float
    weight: torch.Tensor
    tolerance: float
    maximum_steps: int
    init_freq_type: int
    init_freq_maximum: float
    initial_frequencies: torch.Tensor
    device: str
    fp_precision: str
    n_dc: int
    is_real_input: bool
    display: str
    display_interval: int
    steps: int
    omega_history: torch.Tensor
    difference_history: torch.Tensor
    coefficient_spectrum: torch.Tensor
    phi: torch.Tensor
    residual: torch.Tensor
    data_spectrum_norm: float
    spatial_start: int = 0
    world_size: int = 1

    @property
    def spatial_stop(self) -> int:
        return self.spatial_start + self.phi.shape[0]


@dataclass(slots=True)
class Progress:
    """Snapshot supplied to an RVMD output callback."""

    step: int
    difference: float
    omega: np.ndarray
    elapsed_time: float
    stop_reason: str = ""


@dataclass(slots=True)
class HilbertAnalysis:
    """Hilbert spectral analysis of selected RVMD coefficients."""

    sample_rate: float
    time: np.ndarray
    mode_indices: np.ndarray
    center_frequencies: np.ndarray
    method: str
    analytic_signal: np.ndarray
    amplitude: np.ndarray
    phase: np.ndarray
    instantaneous_frequency: np.ndarray
    instantaneous_energy: np.ndarray
    frequency_edges: np.ndarray
    frequency_bins: np.ndarray
    frequency_limits: tuple[float, float]
    hilbert_spectrum: scipy.sparse.csr_matrix
    marginal_spectrum: np.ndarray
    mode_energy: np.ndarray
    hilbert_energy: np.ndarray
    amplitude_threshold: float
    mirror_extension: bool
