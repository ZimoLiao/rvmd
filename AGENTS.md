# RVMD repository guide for coding agents

## Repository scope

The root-level interfaces documented by this repository are:

- `rvmd.m`: decomposition, progress reporting, checkpointing, and restart.
- `rvmdplot.m`: convergence and center-frequency history plots.
- `rvmdhilbert.m`: Hilbert spectral analysis of `mode.c`.
- `rvmdhilbertplot.m`: plots an `rvmdhilbert` result.
- `src/rvmdpy`: matching typed Python API implemented with PyTorch.
- `examples_python/multigpu_rvmd.py`: torchrun driver for one spatially
  sharded decomposition.

The `case*` and `tutorial_CylinderWake` directories are executable examples.
Their local helpers, including the BVH utilities and the tutorial frequency-
spectrum and mirror-extension functions, are example-internal rather than
supported public interfaces. Do not expand or redesign them unless a task
explicitly targets the corresponding example.

User-facing behavior is documented in each root function's help block, in
`README.md`, and in `functionSignatures.json`. Keep those three sources aligned.

## Contracts that changes must preserve

- A new decomposition uses `rvmd(Q,K,Alpha,Name,Value,...)`.
- A continuation uses `rvmd('Restart',state,Name,Value,...)`.
- The reconstruction convention is `mode.phi * mode.c.'`; the transpose is
  nonconjugate.
- `mode.phi`, `mode.c`, and `mode.omega` are sorted by final center frequency.
  `info.Iteration.omega` and restart arrays retain internal mode order.
- Center frequencies are stored in cycles per sample in `[0,0.5]`.
- Real input uses the nonnegative half-spectrum. Complex input uses the full
  shifted spectrum and the distance `abs(f)-omega`.
- `MaximumSteps` is the total sweep limit, including steps already present in
  a restart state. It is not an additional-step count.
- Output callbacks, time-limit checks, and checkpoint writes occur only at
  complete-sweep boundaries. A callback stop must return a restartable state.
- Restart state is opaque to callers. Version 4 stores `residual_n` and
  `dataSpectrumNorm` and does not store the original `Q`. Do not change restart
  fields or the version without an upgrade path and restart-equivalence tests.
- A checkpoint MAT-file contains the variable `restart`. Its previous valid
  version is retained with the `.prev` suffix.
- `rvmdhilbert` requires no Signal Processing Toolbox. It constructs the
  analytic signal in the Fourier domain for real coefficients and treats
  complex coefficients as already analytic-valued.
- `rvmdhilbert` reports frequencies in the units of its `sampleRate` input and
  stores the combined Hilbert energy spectrum as a sparse frequency-by-time
  matrix.
- `rvmdhilbert` keeps the existing permissive `MirrorExtension` input contract:
  logical values and finite numeric scalars are accepted; zero is false and a
  nonzero value is true.
- Python uses the same filter, update order, center-frequency definition,
  nonconjugate reconstruction, defaults, and total-step restart semantics.
  Python option and field names use snake_case and mode indices are zero-based.
- `rvmd_distributed` shards only the spatial dimension. Local `phi`, `residual`,
  and `weight` are sharded; `c`, `omega`, and iteration histories are replicated.
  The function accelerates one decomposition and is not batch parallelism.
- Distributed CUDA execution uses one torchrun process per GPU and NCCL. Every
  rank must enter collective operations in the same order.
- A Python single-device checkpoint is one atomic `.pt` file. A distributed
  checkpoint is an atomic directory containing spatial shards and global state.
  Python checkpoints are not MATLAB restart files.

## Change guidelines

- Keep the numerical iteration readable and procedural. Add abstractions only
  when they remove repeated logic or isolate a testable responsibility.
- Do not change the RVMD filter, update ordering, normalization, phase
  convention, convergence metric, or mirror extension as incidental cleanup.
- Preserve real and complex behavior, odd and even sample counts, spatial
  weights, single and double precision, and zero-input handling.
- Do not add a toolbox dependency to the CPU or Hilbert-analysis path. The GPU
  path may use MATLAB `gpuArray` support.
- When a public parameter, output field, default, or error contract changes,
  update the function help, both README language sections,
  `functionSignatures.json`, and tests in the same change.
- Do not claim compatibility with a specific MATLAB, Octave, GPU, or toolbox
  version unless that combination has been executed.
- Keep CPU, single-CUDA, and distributed Python execution in the shared solver
  implementation. Do not fork the numerical algorithm into backend copies.
- Do not gather a distributed residual or full spatial mode merely to save a
  checkpoint. Preserve sharded checkpoint I/O and restart resharding.
- Collective changes require a two-rank Gloo regression locally and a real NCCL
  run before claiming multi-GPU hardware validation.

## Verification

From the repository root, run the GNU Octave compatibility suite:

```bash
octave-cli --no-gui --quiet --eval "addpath('tests'); run_octave_tests;"
```

In MATLAB, run both function-based suites:

```matlab
results = runtests({'tests/test_rvmd.m', 'tests/test_rvmdhilbert.m'});
assertSuccess(results)
```

Also validate repository metadata and whitespace:

```bash
jq empty functionSignatures.json
git diff --check
```

Run the Python suite from the repository root:

```bash
python -m ruff check src tests_python examples_python benchmarks
python -m ruff format --check src tests_python examples_python benchmarks
python -m pytest tests_python
```

The Python suite includes CPU, checkpoint/restart, Hilbert analysis, an Octave
numerical comparison when Octave is available, CUDA when a GPU is available,
and a two-process Gloo simulation of the spatially sharded solver.

On a machine with at least two CUDA GPUs, run the same distributed correctness
worker through NCCL before claiming multi-GPU hardware validation:

```bash
RVMD_TEST_DEVICE=cuda torchrun --standalone --nproc-per-node=2 \
  tests_python/distributed_worker.py --checkpoint /tmp/rvmd-nccl-checkpoint
```

Restart-related changes require a split-run comparison: stopping and restarting
must reproduce the corresponding uninterrupted run for the same precision and
device. Hilbert-analysis changes require at least a known real tone and a known
signed complex tone.

## Known boundary

`Ctrl+C` is an immediate host interruption and is not a graceful-stop API.
Use `OutputFcn` or `TimeLimit` for a returned consistent state. Periodic
checkpoints provide recovery from the most recently completed checkpoint.

The Python tests can validate distributed algorithm and checkpoint behavior
with two CPU processes, but that is not evidence of NCCL multi-GPU performance.
State the actual GPU hardware used whenever reporting a performance result.
