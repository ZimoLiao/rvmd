from __future__ import annotations

import argparse
import os
from pathlib import Path

import numpy as np
import torch
import torch.distributed as dist

from rvmdpy import gather_mode, rvmd, rvmd_distributed


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint", type=Path, required=True)
    arguments = parser.parse_args()
    device_type = os.environ.get("RVMD_TEST_DEVICE", "cpu")
    if device_type not in {"cpu", "cuda"}:
        raise ValueError("RVMD_TEST_DEVICE must be 'cpu' or 'cuda'")
    if device_type == "cuda":
        local_rank = int(os.environ["LOCAL_RANK"])
        torch.cuda.set_device(local_rank)
        device = f"cuda:{local_rank}"
        backend = "nccl"
    else:
        device = "cpu"
        backend = "gloo"
    dist.init_process_group(backend)
    rank = dist.get_rank()
    assert dist.get_world_size() == 2

    sample_count = 40
    n = np.arange(sample_count)
    full_q = np.vstack(
        [
            np.cos(2 * np.pi * 3 * n / sample_count),
            np.sin(2 * np.pi * 7 * n / sample_count),
            0.4 * np.cos(2 * np.pi * 11 * n / sample_count),
            np.sin(2 * np.pi * 5 * n / sample_count),
            0.2 * np.cos(2 * np.pi * 13 * n / sample_count),
        ]
    )
    bounds = ((0, 2), (2, 5))[rank]
    local_q = full_q[bounds[0] : bounds[1]]
    common = dict(
        k=3,
        alpha=40,
        tolerance=0,
        fp_precision="double",
        initial_frequencies=[0.05, 0.15, 0.28],
        device=device,
    )
    try:
        rvmd_distributed(local_q, maximum_steps=2 + rank, **common)
    except ValueError as exception:
        assert "all ranks must use the same" in str(exception)
    else:
        raise AssertionError("inconsistent distributed settings were accepted")

    distributed_full, full_info, _ = rvmd_distributed(
        local_q, maximum_steps=7, **common
    )
    _, _, partial = rvmd_distributed(local_q, maximum_steps=3, **common)
    resumed, resumed_info, _ = rvmd_distributed(
        restart=partial, maximum_steps=7, tolerance=0, device=device
    )
    assert torch.equal(distributed_full.phi, resumed.phi)
    assert torch.equal(distributed_full.c, resumed.c)
    assert torch.equal(full_info.iteration.omega, resumed_info.iteration.omega)

    checkpointed, checkpoint_info, _ = rvmd_distributed(
        local_q,
        maximum_steps=3,
        checkpoint_file=arguments.checkpoint,
        checkpoint_interval=2,
        **common,
    )
    assert checkpoint_info.last_checkpoint_step == 3
    assert checkpointed.phi.shape[0] == local_q.shape[0]

    gathered = gather_mode(distributed_full, destination=0)
    if rank == 0:
        assert gathered is not None
        single, _, _ = rvmd(full_q, maximum_steps=7, **common)
        assert torch.allclose(gathered.phi, single.phi, atol=2e-12, rtol=2e-12)
        assert torch.allclose(gathered.c, single.c, atol=2e-12, rtol=2e-12)
        assert torch.allclose(gathered.omega, single.omega, atol=2e-12, rtol=2e-12)
    dist.barrier()
    dist.destroy_process_group()


if __name__ == "__main__":
    main()
