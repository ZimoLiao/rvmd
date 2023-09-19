# Reduced-order Variational Mode Decomposition (RVMD) in MATLAB

`rvmd()` is a MATLAB implementation of reduced-order variational mode decomposition (RVMD). 

RVMD is a multi-component extension of variational mode decomposition (VMD), which combines the low-order representation in modal analysis (for fluid dynamics). RVMD can adaptively extract **low-order dynamics featured by transient or non-stationary properties** in space-time data.

## Files

| File/Folder                  | Description                               |
| ---------------------------- | ----------------------------------------- |
| `rvmd.m`                     | RVMD function in MATLAB                            |
| `functionSignatures.json`    | MATLAB function signatures for `rvmd()`   |
| `case1_NonstationarySignal/` | Test case for a 1D non-stationary signal |
| `case2_LorenzAttractor/`     | Test case for Lorenz Attractor            |
| `case3_CylinderWake/`        | Test case for transient cylinder wake     |
| `case4_MotionCapture/`       | Test case for motion capture data         |

## Reference

[1] Liao, Z.-M., Zhao, Z., Chen, L.-B., Wan, Z.-H., Liu, N.-S. & Lu, X.-Y. 2023 Reduced-order variational mode decomposition to reveal transient and non-stationary dynamics in fluid flows. J. Fluid Mech., 966, A7. [doi:10.1017/jfm.2023.435](https://doi.org/10.1017/jfm.2023.435)
