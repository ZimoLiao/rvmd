"""Matplotlib diagnostics for RVMD and Hilbert spectral analysis."""

from __future__ import annotations

import numpy as np

from .types import HilbertAnalysis, Info


def rvmdplot(info: Info, *, sample_rate: float = 1.0):
    """Plot convergence difference and unsorted center-frequency histories.

    Args:
        info: Termination and iteration history returned by an RVMD call.
        sample_rate: Positive multiplier converting cycles/sample to the desired
            frequency unit. Default 1 leaves normalized frequency unchanged.

    Returns:
        A Matplotlib figure with convergence and frequency-history axes.
    """

    import matplotlib.pyplot as plt

    try:
        sample_rate = float(sample_rate)
    except (TypeError, ValueError) as exception:
        raise ValueError("sample_rate must be a positive finite scalar") from exception
    if not np.isfinite(sample_rate) or sample_rate <= 0:
        raise ValueError("sample_rate must be a positive finite scalar")
    difference = info.iteration.difference.detach().cpu().double().numpy()
    omega = info.iteration.omega.detach().cpu().double().numpy() * sample_rate
    steps = info.iteration.steps
    figure, axes = plt.subplots(2, 1, constrained_layout=True)
    if difference.size:
        axes[0].semilogy(
            np.arange(1, steps + 1), np.maximum(difference, np.finfo(float).tiny)
        )
    axes[0].grid(True)
    axes[0].set(
        xlabel="Iteration", ylabel="Difference", title=f"RVMD: {info.stop_reason}"
    )
    axes[1].plot(np.arange(steps + 1), omega.T)
    axes[1].grid(True)
    axes[1].set(
        xlabel="Iteration",
        ylabel=(
            "Center frequency (cycles/sample)"
            if sample_rate == 1
            else "Center frequency"
        ),
    )
    return figure


def rvmdhilbertplot(analysis: HilbertAnalysis):
    """Plot Hilbert energy, marginal energy, and instantaneous frequencies.

    Args:
        analysis: Result returned by :func:`rvmdpy.rvmdhilbert`.

    Returns:
        A Matplotlib figure containing spectrum, marginal, and frequency axes.
    """

    import matplotlib.pyplot as plt

    energy = analysis.hilbert_spectrum.toarray()
    positive = energy[energy > 0]
    if positive.size:
        floor = max(float(np.min(positive)), float(np.max(positive)) * 1e-12)
        log_energy = np.log10(np.maximum(energy, floor))
    else:
        log_energy = np.zeros_like(energy)
    figure = plt.figure(constrained_layout=True)
    grid = figure.add_gridspec(2, 2)
    spectrum_axis = figure.add_subplot(grid[:, 0])
    marginal_axis = figure.add_subplot(grid[0, 1])
    frequency_axis = figure.add_subplot(grid[1, 1])
    image = spectrum_axis.imshow(
        log_energy,
        origin="lower",
        aspect="auto",
        extent=(
            analysis.time[0],
            analysis.time[-1] if analysis.time.size > 1 else analysis.time[0],
            analysis.frequency_edges[0],
            analysis.frequency_edges[-1],
        ),
    )
    spectrum_axis.set(
        xlabel="Time", ylabel="Frequency", title="Hilbert energy spectrum"
    )
    figure.colorbar(image, ax=spectrum_axis)
    marginal_axis.plot(analysis.marginal_spectrum, analysis.frequency_bins)
    marginal_axis.grid(True)
    marginal_axis.set(
        xlabel="Integrated energy", ylabel="Frequency", title="Marginal spectrum"
    )
    frequency_axis.plot(analysis.time, analysis.instantaneous_frequency)
    frequency_axis.grid(True)
    frequency_axis.set(
        xlabel="Time", ylabel="Frequency", title="Instantaneous frequencies"
    )
    return figure
