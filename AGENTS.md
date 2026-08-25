# RVMD repository guide for coding agents

## Repository scope

The root-level interfaces documented by this repository are:

- `rvmd.m`: decomposition, progress reporting, checkpointing, and restart.
- `rvmdplot.m`: convergence and center-frequency history plots.
- `rvmdhilbert.m`: Hilbert spectral analysis of `mode.c`.
- `rvmdhilbertplot.m`: plots an `rvmdhilbert` result.

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

Restart-related changes require a split-run comparison: stopping and restarting
must reproduce the corresponding uninterrupted run for the same precision and
device. Hilbert-analysis changes require at least a known real tone and a known
signed complex tone.

## Known boundary

`Ctrl+C` is an immediate host interruption and is not a graceful-stop API.
Use `OutputFcn` or `TimeLimit` for a returned consistent state. Periodic
checkpoints provide recovery from the most recently completed checkpoint.
