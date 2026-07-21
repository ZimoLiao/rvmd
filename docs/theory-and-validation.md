# RVMD theory and implementation contract

This document fixes the mathematical convention implemented by `rvmd.m`.
It is intended to prevent the real- and complex-valued code paths from
silently drifting apart.

## Inner product and model

The spatial inner product is linear in its first argument:

\[
\langle f,g\rangle_{x,W}
=\int_\Omega f(x)\overline{g(x)}W(x)\,\mathrm dx,
\qquad W(x)>0.
\]

RVMD approximates

\[
q(x,t)\approx\sum_{k=1}^K\phi_k(x)c_k(t),
\qquad \|\phi_k\|_{x,W}=1.
\]

The MATLAB input `Weight` is converted to a column vector and divided by
its mean. This leaves relative spatial weighting unchanged and makes the
meaning of `Alpha` independent of a global rescaling of the supplied
weights.

## Real-valued problem

For real `Q`, `phi` and `c`, Hermitian symmetry allows the problem to be
written on nonnegative frequencies:

\[
\min \int_0^\infty
\left\|\hat q-\sum_k\phi_k\hat c_k\right\|_{x,W}^2\,\mathrm d\omega
+\sum_k\int_0^\infty
2\alpha(\omega-\omega_k)^2|\hat c_k|^2\,\mathrm d\omega.
\]

Holding all other modes fixed and writing their exclusion residual as
`r_k`, one Gauss--Seidel block update is

\[
b_k=\operatorname{Re}\int_0^\infty
\overline{\hat c_k}\hat r_k\,\mathrm d\omega,
\qquad
\phi_k=\frac{b_k}{\|b_k\|_{x,W}},
\]

\[
\hat c_k(\omega)=
\frac{\langle\hat r_k(\cdot,\omega),\phi_k\rangle_{x,W}}
{1+2\alpha(\omega-\omega_k)^2},
\]

\[
\omega_k=
\frac{\int_0^\infty\omega|\hat c_k(\omega)|^2\,\mathrm d\omega}
{\int_0^\infty|\hat c_k(\omega)|^2\,\mathrm d\omega}.
\]

The denominator is therefore `1 + 2*Alpha*(...)^2`. The former
`develop` implementation used `4*Alpha`; that is inconsistent with both
the objective and equation (3.14) of the published paper. For this filter,
the distance between its two half-power points is

\[
\Delta=\sqrt{(2\sqrt{2}-2)/\alpha}.
\]

## Complex-valued extension

For complex `Q`, `phi`, and `c`, the repository implements the explicitly
defined full-spectrum objective

\[
\min \int_{-\infty}^{\infty}
\left\|\hat q-\sum_k\phi_k\hat c_k\right\|_{x,W}^2\,\mathrm d\omega
+\sum_k\int_{-\infty}^{\infty}
2\alpha(|\omega|-\omega_k)^2|\hat c_k|^2\,\mathrm d\omega.
\]

The updates are

\[
b_k=\int_{-\infty}^{\infty}
\overline{\hat c_k}\hat r_k\,\mathrm d\omega,
\qquad
\phi_k=\frac{b_k}{\|b_k\|_{x,W}},
\]

\[
\hat c_k(\omega)=
\frac{\langle\hat r_k(\cdot,\omega),\phi_k\rangle_{x,W}}
{1+2\alpha(|\omega|-\omega_k)^2},
\]

\[
\omega_k=
\frac{\int_{-\infty}^{\infty}|\omega||\hat c_k(\omega)|^2\,\mathrm d\omega}
{\int_{-\infty}^{\infty}|\hat c_k(\omega)|^2\,\mathrm d\omega}.
\]

No real-part projection is used in the complex spatial update. With the
inner-product convention above, the numerator of the coefficient update is
`r.' * conj(phi)`, not a conjugate transpose of `r`.

The use of `abs(omega)` means that one nonnegative center describes energy
near the pair `+omega_k` and `-omega_k`. This is the intended extension in
this repository. It does not assign a signed center frequency.

## Discrete FFT convention

Frequencies are normalized cycles per sample. Mirror extension doubles the
time length before the FFT.

- Real data use bins `0:1/(2T):0.5` and reconstruct the negative half by
  Hermitian symmetry.
- Complex data use all shifted bins `-0.5:1/(2T):0.5-1/(2T)`.
- In the real one-sided sums, DC and Nyquist receive half weight. Interior
  bins receive unit weight. These are the exact trapezoidal multiplicities
  obtained by reducing the full Hermitian spectrum; the weights cancel out
  of the pointwise `c` update but matter in the `phi` and `omega` updates.

## Normalization, phase, and uniqueness

Unit weighted norm removes the arbitrary amplitude exchange between
`phi_k` and `c_k`. A complex mode still has the invariance

\[
(\phi_k,c_k)\mapsto
(\mathrm e^{\mathrm i\theta}\phi_k,
 \mathrm e^{-\mathrm i\theta}c_k).
\]

`rvmd.m` removes this reporting ambiguity by making the largest-magnitude
entry of every complex spatial mode real and positive. The analogous sign
rule is used for real modes. Output modes are then sorted by increasing
center frequency.

This gauge convention does **not** make the non-convex RVMD optimization
globally unique. Different initial frequencies, random initialization, or
degenerate/equal-frequency data can lead to different stationary points or
mode bases. The repository promises one canonical implementation file and
a deterministic representation of each returned mode, not a theorem of
global-solution uniqueness.

## Degenerate blocks and restart

If a spatial projection is zero, every unit-norm spatial mode solves that
block subproblem. The implementation retains a normalized previous mode or
uses a deterministic basis vector. If a coefficient spectrum has zero
energy, its previous center frequency is retained. These rules keep zero
data and over-specified decompositions finite.

Restart state stores the unsorted internal Gauss--Seidel order and the live
exclusion residual, so current-format checkpoints preserve the floating-
point trajectory exactly when resumed on the same device and precision.
Legacy states without a saved residual remain loadable and reconstruct it
from the saved modes. On restart, `MaximumSteps` is the total iteration cap,
not a number of additional steps. A staged run to steps 3 and 6 must match a
single run to step 6.

## Evidence ledger

| Claim | Status | Evidence | Failure guarded against |
|---|---|---|---|
| Real block updates above minimize their coordinate subproblems | proved here and in the paper | expansion, Cauchy--Schwarz, first variation | misplaced conjugate or real projection |
| Complex block updates minimize the stated full-spectrum subproblems | proved here | same first-variation calculation over the full spectrum | conjugating the wrong factor |
| Filter denominator is `1 + 2*Alpha*d^2` | known and independently derived | JFM (2023) eq. 3.14; derivative of the stated objective | former factor-of-two bug |
| Weighted modes have unit `W/mean(W)` norm | tested | `tests/test_rvmd.m` | row-weight broadcasting and unweighted normalization |
| Zero input and zero-energy modes remain finite | tested | `tests/test_rvmd.m` | division by zero and NaN centers |
| Restart is trajectory-equivalent to a single run | tested | `tests/test_rvmd.m` | lost precision/options/history |
| Real and complex synthetic tones recover frequency magnitude | experiment-only | `tests/test_rvmd.m` | wrong spectrum branch or signed-frequency handling |

Primary source: [Liao et al., *Journal of Fluid Mechanics* 966 (2023),
A7](https://doi.org/10.1017/jfm.2023.435).
