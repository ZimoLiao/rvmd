# Reduced-order Variational Mode Decomposition (RVMD) in MATLAB

`rvmd()` is a MATLAB implementation of reduced-order variational mode decomposition (RVMD). RVMD is a multi-component extension of variational mode decomposition (VMD), which combines the low-order representation in modal analysis (for fluid dynamics). RVMD can adaptively extract **low-order dynamics featured by transient or nonstationary properties** in space-time data.

## Files

| File              | Description                                |
| ----------------- | ------------------------------------------ |
| `rvmd.m`          | RVMD function                              |
| `test_1d.m`       | 1D nonstationary signal decomposition test |
| `test_cylinder.m` | transient cylinder wakes test              |
| `test_lorenz.m`   | lorenz system                              |

## Reference

[1] Liao, Z.-M., Zhao, Z., Chen, L.-B., Wan, Z.-H., Liu, N.-S. & Lu, X.-Y. 2023 Reduced-order variational mode decomposition to reveal transient and non-stationary dynamics in fluid flows. J. Fluid Mech., 966, A7. doi:10.1017/jfm.2023.435
