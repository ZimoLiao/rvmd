from __future__ import annotations

import importlib.resources
from pathlib import Path

import numpy as np
import pytest
import torch
from scipy.io import savemat

from examples_python.cylinder_jfm2023 import (
    PAPER_INITIAL_FREQUENCIES,
    PAPER_MODE_COUNT,
    load_cylinder,
    resolve_device,
)


def test_jfm_example_loads_matlab_order_and_removes_spatial_means(
    tmp_path: Path,
) -> None:
    snapshots = np.arange(24, dtype=float).reshape((4, 6), order="F")
    path = tmp_path / "cylinder.mat"
    savemat(path, {"S": 4, "T": 6, "velocity": snapshots})

    loaded = load_cylinder(path)

    assert loaded.shape == (4, 6)
    assert np.allclose(loaded, snapshots - snapshots.mean(axis=1, keepdims=True))
    assert np.allclose(loaded.mean(axis=1), 0)


def test_jfm_example_defaults_and_device_resolution() -> None:
    assert PAPER_INITIAL_FREQUENCIES.shape == (PAPER_MODE_COUNT,)
    assert PAPER_INITIAL_FREQUENCIES[0] == 0
    assert PAPER_INITIAL_FREQUENCIES[-1] == pytest.approx(0.15)
    assert resolve_device("cpu") == "cpu"
    assert resolve_device("auto") == ("cuda" if torch.cuda.is_available() else "cpu")


def test_installed_package_declares_typing_support() -> None:
    marker = importlib.resources.files("rvmdpy").joinpath("py.typed")
    assert marker.is_file()
