from __future__ import annotations

import matplotlib
import numpy as np

matplotlib.use("Agg")

import matplotlib.pyplot as plt

from rvmdpy import rvmd, rvmdhilbert, rvmdhilbertplot, rvmdplot


def test_diagnostic_plots_render_without_a_display() -> None:
    sample_count = 32
    time = np.arange(sample_count)
    q = np.vstack(
        (
            np.cos(2 * np.pi * 3 * time / sample_count),
            np.sin(2 * np.pi * 7 * time / sample_count),
        )
    )
    mode, info, _ = rvmd(
        q,
        2,
        30,
        maximum_steps=3,
        tolerance=0,
        fp_precision="double",
        initial_frequencies=[3 / sample_count, 7 / sample_count],
    )
    analysis = rvmdhilbert(mode, sample_rate=2.0)

    convergence_figure = rvmdplot(info, sample_rate=2.0)
    hilbert_figure = rvmdhilbertplot(analysis)

    assert len(convergence_figure.axes) == 2
    assert len(hilbert_figure.axes) == 4
    convergence_figure.canvas.draw()
    hilbert_figure.canvas.draw()
    plt.close("all")
