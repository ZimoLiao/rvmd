"""Run the JFM 2023 transient-cylinder RVMD and Hilbert-analysis workflow.

The numerical defaults reproduce the paper setup in normalized sample units:
K=10, alpha=1000, tolerance=0.002, and ten initial center frequencies
uniformly distributed over [0, 0.15].  With sample_rate=4, that interval is
[0, 0.6] in the paper's physical frequency units.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import scipy.sparse
import torch
from scipy.io import loadmat

from rvmdpy import (
    HilbertAnalysis,
    Info,
    Mode,
    Progress,
    rvmd,
    rvmdhilbert,
    rvmdhilbertplot,
    rvmdplot,
)

PAPER_MODE_COUNT = 10
PAPER_ALPHA = 1000.0
PAPER_TOLERANCE = 0.002
PAPER_SAMPLE_RATE = 4.0
PAPER_INITIAL_FREQUENCIES = np.linspace(0.0, 0.15, PAPER_MODE_COUNT)
# Result from the original 2023 implementation with the settings above.
REFERENCE_PHYSICAL_FREQUENCIES = np.array(
    [
        0.0039296,
        0.1310116,
        0.1313817,
        0.1553421,
        0.1579847,
        0.3059522,
        0.3077050,
        0.4666592,
        0.4692709,
        0.6265233,
    ]
)


def default_data_path() -> Path:
    return (
        Path(__file__).resolve().parents[1] / "case3_CylinderWake" / "data_cylinder.mat"
    )


def load_cylinder(path: Path) -> np.ndarray:
    """Load and mean-subtract the paper's S-by-T velocity snapshot matrix."""

    data = loadmat(path)
    required = {"S", "T", "velocity"}
    missing = sorted(required.difference(data))
    if missing:
        raise ValueError(f"{path} is missing variables: {', '.join(missing)}")
    spatial_size = int(np.asarray(data["S"]).squeeze())
    sample_count = int(np.asarray(data["T"]).squeeze())
    velocity = np.asarray(data["velocity"])
    if velocity.size != spatial_size * sample_count:
        raise ValueError(
            "velocity size does not match S*T: "
            f"{velocity.size} != {spatial_size}*{sample_count}"
        )
    snapshots = velocity.reshape((spatial_size, sample_count), order="F")
    if not np.isfinite(snapshots).all():
        raise ValueError("velocity contains non-finite values")
    return snapshots - snapshots.mean(axis=1, keepdims=True)


def resolve_device(requested: str) -> str:
    """Resolve ``auto`` and reject an unavailable requested CUDA backend."""

    if requested == "auto":
        return "cuda" if torch.cuda.is_available() else "cpu"
    if requested == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("--device cuda was requested, but PyTorch cannot use CUDA")
    return requested


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=default_data_path())
    parser.add_argument("--output", type=Path, default=Path("rvmd_jfm2023_output"))
    parser.add_argument("--restart", type=Path, help="resume a Python .pt checkpoint")
    parser.add_argument("--device", choices=("auto", "cpu", "cuda"), default="auto")
    parser.add_argument("--precision", choices=("single", "double"), default="single")
    parser.add_argument("--maximum-steps", type=int, default=500)
    parser.add_argument(
        "--time-limit",
        type=float,
        default=math.inf,
        help="seconds for this call; checked at complete-sweep boundaries",
    )
    parser.add_argument("--checkpoint-interval", type=int, default=50)
    parser.add_argument("--display-interval", type=int, default=20)
    parser.add_argument(
        "--stop-file",
        type=Path,
        help="stop safely when this file exists (default: OUTPUT/STOP_RVMD)",
    )
    parser.add_argument("--no-plots", action="store_true")
    return parser.parse_args()


def save_results(
    output: Path,
    mode: Mode,
    info: Info,
    analysis: HilbertAnalysis,
    *,
    make_plots: bool,
) -> None:
    output.mkdir(parents=True, exist_ok=True)
    physical_frequencies = PAPER_SAMPLE_RATE * mode.omega.numpy()
    maximum_frequency_error = (
        float(np.max(np.abs(physical_frequencies - REFERENCE_PHYSICAL_FREQUENCIES)))
        if info.iteration.converged
        else None
    )
    torch.save(
        {
            "phi": mode.phi,
            "c": mode.c,
            "omega": mode.omega,
            "energy": mode.energy,
        },
        output / "mode.pt",
    )
    np.savez_compressed(
        output / "hilbert_analysis.npz",
        time=analysis.time,
        center_frequencies=analysis.center_frequencies,
        amplitude=analysis.amplitude,
        phase=analysis.phase,
        instantaneous_frequency=analysis.instantaneous_frequency,
        instantaneous_energy=analysis.instantaneous_energy,
        frequency_edges=analysis.frequency_edges,
        frequency_bins=analysis.frequency_bins,
        marginal_spectrum=analysis.marginal_spectrum,
    )
    scipy.sparse.save_npz(output / "hilbert_spectrum.npz", analysis.hilbert_spectrum)
    summary = {
        "device": info.device,
        "precision": info.fp_precision,
        "stop_reason": info.stop_reason,
        "converged": info.iteration.converged,
        "steps": info.iteration.steps,
        "difference": float(info.iteration.difference[-1])
        if info.iteration.difference.numel()
        else None,
        "elapsed_seconds": info.elapsed_time,
        "sample_rate": PAPER_SAMPLE_RATE,
        "center_frequencies_cycles_per_sample": mode.omega.numpy().tolist(),
        "center_frequencies_physical": physical_frequencies.tolist(),
        "reference_maximum_frequency_error": maximum_frequency_error,
        "reconstruction": "phi @ c.T (nonconjugating transpose)",
    }
    with (output / "summary.json").open("w", encoding="utf-8") as stream:
        json.dump(summary, stream, indent=2)
        stream.write("\n")

    if make_plots:
        convergence = rvmdplot(info, sample_rate=PAPER_SAMPLE_RATE)
        convergence.savefig(output / "convergence.png", dpi=180)
        hilbert = rvmdhilbertplot(analysis)
        hilbert.savefig(output / "hilbert_spectrum.png", dpi=180)
        import matplotlib.pyplot as plt

        plt.close(convergence)
        plt.close(hilbert)


def main() -> None:
    arguments = parse_arguments()
    device = resolve_device(arguments.device)
    arguments.output.mkdir(parents=True, exist_ok=True)
    stop_file = arguments.stop_file or arguments.output / "STOP_RVMD"
    checkpoint = arguments.output / "restart.pt"

    def monitor(progress: Progress, phase: str) -> bool:
        if phase == "iter" and progress.step % arguments.display_interval == 0:
            physical = np.sort(progress.omega) * PAPER_SAMPLE_RATE
            print(
                f"step={progress.step} difference={progress.difference:.4e} "
                f"frequency={np.array2string(physical, precision=5)}",
                flush=True,
            )
        return phase in {"init", "iter"} and stop_file.is_file()

    if arguments.restart is None:
        snapshots = load_cylinder(arguments.data)
        mode, info, _ = rvmd(
            snapshots,
            PAPER_MODE_COUNT,
            PAPER_ALPHA,
            tolerance=PAPER_TOLERANCE,
            maximum_steps=arguments.maximum_steps,
            initial_frequencies=PAPER_INITIAL_FREQUENCIES,
            device=device,
            fp_precision=arguments.precision,
            display="off",
            output_fcn=monitor,
            time_limit=arguments.time_limit,
            checkpoint_file=checkpoint,
            checkpoint_interval=arguments.checkpoint_interval,
        )
    else:
        mode, info, _ = rvmd(
            restart=arguments.restart,
            maximum_steps=arguments.maximum_steps,
            device=device,
            display="off",
            output_fcn=monitor,
            time_limit=arguments.time_limit,
            checkpoint_file=checkpoint,
            checkpoint_interval=arguments.checkpoint_interval,
        )

    analysis = rvmdhilbert(mode, PAPER_SAMPLE_RATE, mirror_extension=True)
    save_results(
        arguments.output,
        mode,
        info,
        analysis,
        make_plots=not arguments.no_plots,
    )
    physical = PAPER_SAMPLE_RATE * mode.omega.numpy()
    print(
        f"finished: reason={info.stop_reason} steps={info.iteration.steps} "
        f"elapsed={info.elapsed_time:.3f}s",
        flush=True,
    )
    print(f"physical center frequencies: {np.array2string(physical, precision=7)}")
    if info.iteration.converged:
        reference_error = np.max(np.abs(physical - REFERENCE_PHYSICAL_FREQUENCIES))
        print(f"difference from the 2023 reference run: {reference_error:.3e}")
    print(f"results: {arguments.output.resolve()}")


if __name__ == "__main__":
    main()
