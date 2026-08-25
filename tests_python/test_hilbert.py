from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse
import torch

from rvmdpy import Mode, rvmdhilbert


def make_mode(coefficient: np.ndarray, omega: np.ndarray) -> Mode:
    coefficient = np.asarray(coefficient)
    if coefficient.ndim == 1:
        coefficient = coefficient[:, None]
    return Mode(
        phi=torch.ones((1, coefficient.shape[1])),
        c=torch.from_numpy(coefficient),
        omega=torch.as_tensor(omega).reshape(-1),
        energy=torch.from_numpy(np.sum(np.abs(coefficient) ** 2, axis=0)),
        global_spatial_size=1,
    )


def test_real_tone_instantaneous_frequency_and_sparse_spectrum() -> None:
    sample_rate = 256
    time = np.arange(256) / sample_rate
    frequency = 32
    mode = make_mode(np.cos(2 * np.pi * frequency * time), [frequency / sample_rate])
    analysis = rvmdhilbert(
        mode, sample_rate, mirror_extension=False, frequency_bins=128
    )
    np.testing.assert_allclose(analysis.amplitude[2:-2], 1, atol=1e-10)
    np.testing.assert_allclose(
        analysis.instantaneous_frequency[2:-2], frequency, atol=1e-10
    )
    assert scipy.sparse.issparse(analysis.hilbert_spectrum)
    assert analysis.hilbert_spectrum.shape == (128, 256)
    peak = np.argmax(analysis.marginal_spectrum)
    assert abs(analysis.frequency_bins[peak] - frequency) <= sample_rate / 128


def test_complex_tone_keeps_signed_frequency() -> None:
    sample_rate = 200
    time = np.arange(200) / sample_rate
    frequency = -25
    mode = make_mode(
        np.exp(2j * np.pi * frequency * time), [abs(frequency) / sample_rate]
    )
    analysis = rvmdhilbert(mode, sample_rate, frequency_bins=100)
    assert analysis.method == "complex-coefficient"
    assert analysis.frequency_limits == (-100, 100)
    np.testing.assert_allclose(
        analysis.instantaneous_frequency[2:-2], frequency, atol=1e-10
    )


def test_zero_amplitude_mask_and_mode_selection() -> None:
    zero = make_mode(np.zeros((32, 2)), [0, 0.1])
    analysis = rvmdhilbert(zero, 10)
    assert np.isnan(analysis.instantaneous_frequency).all()
    assert analysis.hilbert_spectrum.nnz == 0
    np.testing.assert_array_equal(analysis.mode_energy, [0, 0])

    time = np.arange(64) / 64
    selected = make_mode(
        np.column_stack((np.cos(2 * np.pi * 4 * time), np.cos(2 * np.pi * 12 * time))),
        [4 / 64, 12 / 64],
    )
    result = rvmdhilbert(selected, 64, mode_indices=[1], mirror_extension=False)
    assert result.analytic_signal.shape == (64, 1)
    np.testing.assert_allclose(result.center_frequencies, [12])


def test_invalid_hilbert_inputs() -> None:
    mode = make_mode(np.ones(16), [0.1])
    with pytest.raises(ValueError):
        rvmdhilbert(mode, 0)
    with pytest.raises(ValueError):
        rvmdhilbert(mode, 1, mode_indices=[0, 0])
    with pytest.raises(ValueError):
        rvmdhilbert(mode, 1, frequency_limits=(1, 0))
    with pytest.raises(ValueError):
        rvmdhilbert(mode, 1, mirror_extension="yes")
