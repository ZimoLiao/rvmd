"""Small repeatable CPU/CUDA throughput benchmark for the Python backend."""

from __future__ import annotations

import argparse
import time

import numpy as np
import torch

from rvmdpy import rvmd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--spatial-size", type=int, default=2048)
    parser.add_argument("--sample-count", type=int, default=1024)
    parser.add_argument("--modes", type=int, default=4)
    parser.add_argument("--steps", type=int, default=5)
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--precision", choices=("single", "double"), default="single")
    arguments = parser.parse_args()

    generator = np.random.default_rng(42)
    dtype = np.float32 if arguments.precision == "single" else np.float64
    q = generator.standard_normal(
        (arguments.spatial_size, arguments.sample_count), dtype=dtype
    )
    initial = np.linspace(0.5 / arguments.modes, 0.5, arguments.modes)
    if arguments.device.startswith("cuda"):
        rvmd(
            q[:32, :64],
            arguments.modes,
            1000,
            maximum_steps=1,
            tolerance=0,
            fp_precision=arguments.precision,
            initial_frequencies=initial,
            device=arguments.device,
        )
        torch.cuda.synchronize()
    start = time.perf_counter()
    _, info, _ = rvmd(
        q,
        arguments.modes,
        1000,
        maximum_steps=arguments.steps,
        tolerance=0,
        fp_precision=arguments.precision,
        initial_frequencies=initial,
        device=arguments.device,
    )
    if arguments.device.startswith("cuda"):
        torch.cuda.synchronize()
    elapsed = time.perf_counter() - start
    print(
        f"device={arguments.device} shape={q.shape} K={arguments.modes} "
        f"steps={info.iteration.steps} elapsed={elapsed:.3f}s "
        f"per_step={elapsed / max(info.iteration.steps, 1):.3f}s"
    )


if __name__ == "__main__":
    main()
