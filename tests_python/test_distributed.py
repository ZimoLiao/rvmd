from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pytest

from rvmdpy import load_distributed_checkpoint


@pytest.mark.skipif(shutil.which("torchrun") is None, reason="torchrun is unavailable")
def test_two_rank_spatial_sharding_restart_checkpoint_and_gather(
    tmp_path: Path,
) -> None:
    checkpoint = tmp_path / "distributed_checkpoint"
    worker = Path(__file__).with_name("distributed_worker.py")
    environment = {
        **os.environ,
        "PYTHONPATH": str(Path(__file__).resolve().parents[1] / "src"),
    }
    subprocess.run(
        [
            "torchrun",
            "--standalone",
            "--nproc-per-node=2",
            str(worker),
            "--checkpoint",
            str(checkpoint),
        ],
        check=True,
        env=environment,
        timeout=120,
    )
    state = load_distributed_checkpoint(checkpoint)
    previous = load_distributed_checkpoint(Path(f"{checkpoint}.prev"))
    assert state.steps == 3
    assert previous.steps == 2
    assert state.spatial_start == 0
    assert state.phi.shape[0] == state.spatial_size == 5
