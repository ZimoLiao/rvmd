"""Toolbox-independent Hilbert spectral analysis for RVMD coefficients."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np
import scipy.fft
import scipy.sparse
import torch

from .types import HilbertAnalysis, Mode


def _numpy(value: Any) -> np.ndarray:
    if isinstance(value, torch.Tensor):
        return value.detach().cpu().numpy()
    return np.asarray(value)


def _analytic_signal(signal: np.ndarray) -> np.ndarray:
    sample_count = signal.shape[0]
    multiplier = np.zeros(sample_count, dtype=signal.real.dtype)
    multiplier[0] = 1
    if sample_count % 2 == 0:
        multiplier[1 : sample_count // 2] = 2
        multiplier[sample_count // 2] = 1
    else:
        multiplier[1 : (sample_count + 1) // 2] = 2
    return scipy.fft.ifft(scipy.fft.fft(signal, axis=0) * multiplier[:, None], axis=0)


def _phase_derivative(phase: np.ndarray, sample_rate: float) -> np.ndarray:
    frequency = np.empty_like(phase)
    if phase.shape[0] == 1:
        frequency.fill(np.nan)
        return frequency
    scale = sample_rate / (2 * np.pi)
    frequency[0] = (phase[1] - phase[0]) * scale
    frequency[-1] = (phase[-1] - phase[-2]) * scale
    if phase.shape[0] > 2:
        frequency[1:-1] = (phase[2:] - phase[:-2]) * (scale / 2)
    return frequency


def rvmdhilbert(
    mode: Mode | Any,
    sample_rate: float,
    *,
    mode_indices: Sequence[int] | None = None,
    mirror_extension: bool | float = True,
    amplitude_threshold: float = 1e-6,
    frequency_bins: int = 256,
    frequency_limits: tuple[float, float] | None = None,
) -> HilbertAnalysis:
    """Analyze instantaneous amplitude, phase, frequency, and Hilbert energy.

    Args:
        mode: RVMD result with finite ``T``-by-``K`` ``c`` and ``K`` ``omega``.
        sample_rate: Positive samples per time unit. Returned frequencies use
            the corresponding inverse-time unit.
        mode_indices: Unique zero-based indices. Default: all modes.
        mirror_extension: Zero disables reflection; any finite nonzero scalar
            enables it for real coefficients with more than two samples.
        amplitude_threshold: Nonnegative fraction of each mode's maximum
            amplitude. Lower-amplitude frequencies become NaN and do not enter
            the Hilbert spectrum. Default: ``1e-6``.
        frequency_bins: Positive Hilbert-spectrum bin count. Default: 256.
        frequency_limits: Increasing frequency pair. Defaults to ``[0,Fs/2]``
            for real coefficients and ``[-Fs/2,Fs/2]`` for complex ones.

    Returns:
        A :class:`HilbertAnalysis`. Its sparse ``B``-by-``T`` spectrum sums
        selected-mode ``amplitude**2`` values in each time-frequency bin.

    Real coefficients use an FFT analytic signal. Complex coefficients are
    treated as already analytic-valued and may have signed frequencies.
    """

    try:
        sample_rate = float(sample_rate)
    except (TypeError, ValueError) as exception:
        raise ValueError("sample_rate must be a positive finite scalar") from exception
    if not np.isfinite(sample_rate) or sample_rate <= 0:
        raise ValueError("sample_rate must be a positive finite scalar")
    coefficient = _numpy(mode.c)
    omega = _numpy(mode.omega).reshape(-1)
    if (
        coefficient.ndim != 2
        or coefficient.size == 0
        or not np.isfinite(coefficient).all()
        or omega.size != coefficient.shape[1]
        or not np.isfinite(omega).all()
    ):
        raise ValueError("mode must contain finite T-by-K c and K-element omega")
    sample_count, mode_count = coefficient.shape
    if mode_indices is None:
        indices = np.arange(mode_count, dtype=np.int64)
    else:
        raw_indices = np.atleast_1d(np.asarray(mode_indices))
        if (
            raw_indices.ndim != 1
            or raw_indices.size == 0
            or not np.issubdtype(raw_indices.dtype, np.integer)
        ):
            raise ValueError("mode_indices must contain unique integer indices")
        indices = raw_indices.astype(np.int64, copy=False)
        if (
            np.any(indices < 0)
            or np.any(indices >= mode_count)
            or np.unique(indices).size != indices.size
        ):
            raise ValueError("mode_indices must contain unique valid indices")
    coefficient = coefficient[:, indices]

    mirror_value = np.asarray(mirror_extension)
    is_numeric = np.issubdtype(mirror_value.dtype, np.number) or np.issubdtype(
        mirror_value.dtype, np.bool_
    )
    if mirror_value.ndim != 0 or not is_numeric or not bool(np.isfinite(mirror_value)):
        raise ValueError("mirror_extension must be a finite scalar")
    use_mirror = bool(mirror_value != 0)
    try:
        amplitude_threshold = float(amplitude_threshold)
    except (TypeError, ValueError) as exception:
        raise ValueError(
            "amplitude_threshold must be nonnegative and finite"
        ) from exception
    if not np.isfinite(amplitude_threshold) or amplitude_threshold < 0:
        raise ValueError("amplitude_threshold must be nonnegative and finite")
    if (
        isinstance(frequency_bins, bool)
        or not isinstance(frequency_bins, (int, np.integer))
        or frequency_bins < 1
    ):
        raise ValueError("frequency_bins must be a positive integer")
    frequency_bins = int(frequency_bins)

    applied_mirror = False
    if np.isrealobj(coefficient):
        if use_mirror and sample_count > 2:
            pad_length = min(sample_count // 2, sample_count - 1)
            extended = np.concatenate(
                (
                    np.flip(coefficient[1 : pad_length + 1], axis=0),
                    coefficient,
                    np.flip(coefficient[-pad_length - 1 : -1], axis=0),
                ),
                axis=0,
            )
            analytic_extended = _analytic_signal(extended)
            analytic = analytic_extended[pad_length : pad_length + sample_count]
            method = "hilbert-mirror"
            applied_mirror = True
        else:
            analytic = _analytic_signal(coefficient)
            method = "hilbert"
    else:
        analytic = coefficient.copy()
        method = "complex-coefficient"

    amplitude = np.abs(analytic)
    phase = np.unwrap(np.angle(analytic), axis=0)
    instantaneous_frequency = _phase_derivative(phase, sample_rate)
    maximum = np.max(amplitude, axis=0)
    valid_amplitude = amplitude >= amplitude_threshold * maximum[None, :]
    valid_amplitude[:, maximum == 0] = False
    instantaneous_frequency[~valid_amplitude] = np.nan
    instantaneous_energy = amplitude**2

    if frequency_limits is None:
        limits = (
            (0.0, sample_rate / 2)
            if np.isrealobj(coefficient)
            else (-sample_rate / 2, sample_rate / 2)
        )
    else:
        values = np.asarray(frequency_limits, dtype=float).reshape(-1)
        if values.size != 2 or not np.isfinite(values).all() or values[0] >= values[1]:
            raise ValueError("frequency_limits must be an increasing finite pair")
        limits = (float(values[0]), float(values[1]))

    edges = np.linspace(limits[0], limits[1], frequency_bins + 1)
    centers = (edges[:-1] + edges[1:]) / 2
    finite_frequency = np.isfinite(instantaneous_frequency)
    bins = np.zeros(instantaneous_frequency.shape, dtype=np.int64)
    bins[finite_frequency] = np.floor(
        (instantaneous_frequency[finite_frequency] - limits[0])
        / (limits[1] - limits[0])
        * frequency_bins
    ).astype(np.int64, copy=False)
    bins[instantaneous_frequency == limits[1]] = frequency_bins - 1
    valid = (
        finite_frequency
        & (bins >= 0)
        & (bins < frequency_bins)
        & np.isfinite(instantaneous_energy)
    )
    time_columns = np.broadcast_to(
        np.arange(sample_count, dtype=np.int64)[:, None], bins.shape
    )
    spectrum = scipy.sparse.coo_matrix(
        (
            instantaneous_energy[valid].astype(float, copy=False),
            (bins[valid], time_columns[valid]),
        ),
        shape=(frequency_bins, sample_count),
    ).tocsr()

    return HilbertAnalysis(
        sample_rate=sample_rate,
        time=np.arange(sample_count, dtype=float) / sample_rate,
        mode_indices=indices,
        center_frequencies=omega[indices].astype(float, copy=False) * sample_rate,
        method=method,
        analytic_signal=analytic,
        amplitude=amplitude,
        phase=phase,
        instantaneous_frequency=instantaneous_frequency,
        instantaneous_energy=instantaneous_energy,
        frequency_edges=edges,
        frequency_bins=centers,
        frequency_limits=limits,
        hilbert_spectrum=spectrum,
        marginal_spectrum=np.asarray(spectrum.sum(axis=1)).ravel() / sample_rate,
        mode_energy=np.sum(np.abs(coefficient.astype(np.complex128)) ** 2, axis=0)
        / sample_rate,
        hilbert_energy=np.sum(instantaneous_energy.astype(float), axis=0) / sample_rate,
        amplitude_threshold=amplitude_threshold,
        mirror_extension=applied_mirror,
    )
