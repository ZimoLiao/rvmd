from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest
import torch

from rvmdpy import load_checkpoint, rvmd


def two_tone(sample_count: int = 64) -> np.ndarray:
    n = np.arange(sample_count)
    return np.vstack(
        (
            np.cos(2 * np.pi * 5 * n / sample_count)
            + 0.2 * np.cos(2 * np.pi * 13 * n / sample_count),
            np.sin(2 * np.pi * 5 * n / sample_count)
            - 0.1 * np.cos(2 * np.pi * 13 * n / sample_count),
        )
    )


def test_real_decomposition_shapes_frequencies_and_reconstruction() -> None:
    sample_count = 128
    n = np.arange(sample_count)
    q = np.array([1, -0.4, 0.7])[:, None] * np.cos(
        2 * np.pi * 8 * n / sample_count
    ) + np.array([0.2, 1, -0.5])[:, None] * np.sin(2 * np.pi * 24 * n / sample_count)
    mode, info, state = rvmd(
        q,
        2,
        500,
        maximum_steps=100,
        tolerance=1e-10,
        initial_frequencies=[8 / 128, 24 / 128],
        fp_precision="double",
    )
    assert mode.phi.shape == (3, 2)
    assert mode.c.shape == (128, 2)
    assert torch.max(torch.abs(mode.omega - torch.tensor([8, 24]) / 128)) < 2 / 128
    assert (
        torch.linalg.matrix_norm(torch.from_numpy(q) - mode.phi @ mode.c.T)
        / torch.linalg.matrix_norm(torch.from_numpy(q))
        < 0.25
    )
    assert info.iteration.steps > 0
    assert state.version == 1


def test_complex_input_preserves_phase_and_signed_data() -> None:
    sample_count = 64
    n = np.arange(sample_count)
    spatial = np.array([1 + 2j, -0.3 + 0.7j, 0.2 - 0.1j])[:, None]
    coefficient = np.exp(2j * np.pi * 3 * n / sample_count) + 0.25 * np.exp(
        -2j * np.pi * 7 * n / sample_count
    )
    q = spatial * coefficient
    mode, _, _ = rvmd(
        q,
        1,
        0,
        maximum_steps=2,
        tolerance=0,
        init_freq_type=0,
        fp_precision="double",
    )
    pivot = torch.argmax(torch.abs(mode.phi[:, 0]))
    assert mode.phi[pivot, 0].real > 0
    assert abs(mode.phi[pivot, 0].imag) < 1e-12
    relative_error = torch.linalg.matrix_norm(
        torch.from_numpy(q) - mode.phi @ mode.c.T
    ) / torch.linalg.matrix_norm(torch.from_numpy(q))
    assert relative_error < 1e-12


def test_weighted_modes_and_zero_input_are_finite() -> None:
    generator = np.random.default_rng(11)
    q = generator.standard_normal((3, 24))
    mode, info, _ = rvmd(
        q,
        2,
        20,
        weight=[1, 2, 5],
        maximum_steps=2,
        tolerance=0,
        fp_precision="double",
    )
    expected_weight = torch.tensor([1, 2, 5], dtype=torch.float64) / (8 / 3)
    weighted_norm = torch.sqrt(
        torch.sum(torch.abs(mode.phi) ** 2 * expected_weight[:, None], dim=0)
    )
    assert torch.allclose(weighted_norm, torch.ones(2, dtype=torch.float64), atol=1e-12)
    assert torch.equal(info.weight, expected_weight)

    zero_mode, zero_info, zero_state = rvmd(
        np.zeros((3, 16)), 2, 10, maximum_steps=3, tolerance=0, fp_precision="double"
    )
    for tensor in (
        zero_mode.phi,
        zero_mode.c,
        zero_mode.omega,
        zero_info.iteration.difference,
        zero_state.coefficient_spectrum,
    ):
        assert torch.isfinite(tensor).all()


@pytest.mark.parametrize("complex_input", [False, True])
@pytest.mark.parametrize("sample_count", [31, 32])
def test_split_restart_is_bitwise_identical(
    complex_input: bool, sample_count: int
) -> None:
    q = two_tone(sample_count)
    if complex_input:
        n = np.arange(sample_count)
        q = q + 0.2j * np.vstack(
            (
                np.sin(2 * np.pi * 4 * n / sample_count),
                np.cos(2 * np.pi * 6 * n / sample_count),
            )
        )
    common = dict(
        k=2,
        alpha=30,
        tolerance=0,
        fp_precision="double",
        initial_frequencies=[0.1, 0.2],
    )
    full, full_info, _ = rvmd(q, maximum_steps=7, **common)
    _, _, state = rvmd(q, maximum_steps=3, **common)
    resumed, resumed_info, _ = rvmd(restart=state, maximum_steps=7, tolerance=0)
    assert torch.equal(resumed.phi, full.phi)
    assert torch.equal(resumed.c, full.c)
    assert torch.equal(resumed.omega, full.omega)
    assert torch.equal(resumed_info.iteration.omega, full_info.iteration.omega)
    assert torch.equal(
        resumed_info.iteration.difference, full_info.iteration.difference
    )


def test_output_function_and_time_limit_return_restartable_states() -> None:
    q = two_tone(48)
    calls: list[tuple[str, int]] = []

    def stop_after_three(progress, phase):
        calls.append((phase, progress.step))
        return phase == "iter" and progress.step >= 3

    full, _, _ = rvmd(
        q,
        2,
        80,
        maximum_steps=6,
        tolerance=0,
        fp_precision="double",
        initial_frequencies=[0.05, 0.2],
    )
    _, stopped_info, state = rvmd(
        q,
        2,
        80,
        maximum_steps=10,
        tolerance=0,
        fp_precision="double",
        initial_frequencies=[0.05, 0.2],
        output_fcn=stop_after_three,
    )
    resumed, _, _ = rvmd(restart=state, maximum_steps=6, tolerance=0)
    assert stopped_info.exit_flag == -1
    assert stopped_info.iteration.steps == 3
    assert calls[0][0] == "init" and calls[-1][0] == "done"
    assert torch.equal(resumed.phi, full.phi)
    assert torch.equal(resumed.c, full.c)

    _, timed_info, timed_state = rvmd(
        q,
        2,
        80,
        maximum_steps=6,
        tolerance=0,
        fp_precision="double",
        initial_frequencies=[0.05, 0.2],
        time_limit=0,
    )
    timed_resumed, _, _ = rvmd(restart=timed_state, maximum_steps=6, tolerance=0)
    assert timed_info.exit_flag == -2
    assert timed_info.iteration.steps == 0
    assert torch.equal(timed_resumed.c, full.c)


def test_atomic_checkpoint_and_previous_file(tmp_path: Path) -> None:
    path = tmp_path / "restart.pt"
    q = two_tone(40)
    full, _, _ = rvmd(
        q,
        2,
        30,
        maximum_steps=6,
        tolerance=0,
        fp_precision="double",
        initial_frequencies=[0.08, 0.2],
    )
    _, info, _ = rvmd(
        q,
        2,
        30,
        maximum_steps=3,
        tolerance=0,
        fp_precision="double",
        initial_frequencies=[0.08, 0.2],
        checkpoint_file=path,
        checkpoint_interval=2,
    )
    state = load_checkpoint(path)
    previous = load_checkpoint(Path(f"{path}.prev"))
    resumed, _, _ = rvmd(restart=state, maximum_steps=6, tolerance=0)
    assert info.last_checkpoint_step == 3
    assert state.steps == 3 and previous.steps == 2
    assert torch.equal(resumed.phi, full.phi)
    assert torch.equal(resumed.c, full.c)


def test_single_precision_dc_and_explicit_initial_frequencies() -> None:
    q = np.vstack((np.ones(32), np.cos(2 * np.pi * 4 * np.arange(32) / 32))).astype(
        np.float32
    )
    initial = [0.0, 0.125]
    mode, info, _ = rvmd(
        q,
        2,
        40,
        n_dc=1,
        initial_frequencies=initial,
        maximum_steps=5,
        tolerance=0,
        fp_precision="single",
    )
    assert mode.phi.dtype == torch.float32
    assert mode.c.dtype == torch.float32
    assert torch.all(info.iteration.omega[0] == 0)
    assert torch.equal(info.iteration.omega[:, 0], torch.tensor(initial))


@pytest.mark.skipif(not torch.cuda.is_available(), reason="CUDA is unavailable")
def test_cuda_matches_cpu_and_restarts_exactly() -> None:
    q = two_tone(128)
    common = dict(
        k=2,
        alpha=500,
        maximum_steps=20,
        tolerance=0,
        fp_precision="double",
        initial_frequencies=[0.04, 0.12],
    )
    cpu, _, _ = rvmd(q, device="cpu", **common)
    gpu, _, _ = rvmd(q, device="cuda", **common)
    _, _, partial = rvmd(q, device="cuda", **{**common, "maximum_steps": 7})
    resumed, _, _ = rvmd(restart=partial, device="cuda", maximum_steps=20, tolerance=0)
    assert torch.allclose(gpu.phi, cpu.phi, atol=1e-12, rtol=1e-12)
    assert torch.allclose(gpu.c, cpu.c, atol=1e-12, rtol=1e-12)
    assert torch.equal(gpu.phi, resumed.phi)
    assert torch.equal(gpu.c, resumed.c)


def test_invalid_inputs_are_rejected() -> None:
    q = np.ones((2, 8))
    with pytest.raises(ValueError):
        rvmd(q, 0, 10)
    with pytest.raises(ValueError):
        rvmd(q, 2, -1)
    with pytest.raises(ValueError):
        rvmd(q, 2, 10, weight=[1, -1])
    with pytest.raises(ValueError):
        rvmd(q, 2, 10, n_dc=3)
    with pytest.raises(ValueError):
        rvmd(q, 2, 10, initial_frequencies=[0.1, 0.2, 0.3])
    with pytest.raises(ValueError):
        rvmd(q, 2, 10, init_freq_type=1.2)
    with pytest.raises(ValueError):
        rvmd(q, 2, 10, weight=[[1, 1], [1, 1]])
    with pytest.raises(ValueError):
        rvmd(q, 2, 10, display=1)


@pytest.mark.skipif(shutil.which("octave-cli") is None, reason="Octave is unavailable")
def test_matches_octave_reference(tmp_path: Path) -> None:
    output = tmp_path / "reference.mat"
    root = Path(__file__).resolve().parents[1]
    expression = (
        f"addpath('{root.as_posix()}'); "
        "T=32; n=0:T-1; q=[cos(2*pi*3*n/T);sin(2*pi*7*n/T)]; "
        "[mode,info]=rvmd(q,2,30,'MaximumSteps',5,'Tolerance',0,"
        "'FPPrecision','double','InitialFrequencies',[0.1;0.2]); "
        "phi=mode.phi; c=mode.c; omega=mode.omega; "
        "omega_hist=info.Iteration.omega; difference=info.Iteration.difference; "
        f"save('-mat7-binary','{output.as_posix()}',"
        "'phi','c','omega','omega_hist','difference');"
    )
    subprocess.run(
        ["octave-cli", "--quiet", "--no-gui", "--eval", expression],
        check=True,
        cwd=root,
    )
    from scipy.io import loadmat

    expected = loadmat(output)
    n = np.arange(32)
    q = np.vstack((np.cos(2 * np.pi * 3 * n / 32), np.sin(2 * np.pi * 7 * n / 32)))
    mode, info, _ = rvmd(
        q,
        2,
        30,
        maximum_steps=5,
        tolerance=0,
        fp_precision="double",
        initial_frequencies=[0.1, 0.2],
    )
    np.testing.assert_allclose(mode.phi, expected["phi"], rtol=2e-8, atol=2e-8)
    np.testing.assert_allclose(mode.c, expected["c"], rtol=2e-8, atol=2e-8)
    np.testing.assert_allclose(
        mode.omega[:, None], expected["omega"], rtol=2e-9, atol=2e-11
    )
    np.testing.assert_allclose(
        info.iteration.omega, expected["omega_hist"], rtol=2e-8, atol=2e-9
    )
    np.testing.assert_allclose(
        info.iteration.difference[None], expected["difference"], rtol=2e-8
    )
