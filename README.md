# Reduced-order Variational Mode Decomposition (RVMD)

[English](#english) · [中文](#中文)

## English

`rvmd.m` is the single MATLAB implementation in this repository for both
real- and complex-valued reduced-order variational mode decomposition.
Given an `S × T` space-time matrix `Q`, it returns spatial modes `phi`,
time-evolution coefficients `c`, and nonnegative center frequencies
`omega` such that

```matlab
Q_reconstructed = mode.phi * mode.c.';
```

The transpose is intentionally nonconjugating (`.'`).

### Installation

Add the repository root to the MATLAB path:

```matlab
addpath('/path/to/rvmd')
```

CPU computation needs no optional MATLAB toolbox. GPU computation requires
Parallel Computing Toolbox.

### Quick start

Real input automatically uses real spatial modes and a nonnegative
half-spectrum:

```matlab
[mode, info] = rvmd(Qreal, 4, 1000, ...
    'FPPrecision', 'double', 'Tolerance', 1e-6);
```

Complex input automatically uses complex modes and the complete shifted
spectrum:

```matlab
[mode, info] = rvmd(Qcomplex, 4, 1000, ...
    'FPPrecision', 'double', 'Tolerance', 1e-6);
```

For complex data, the bandwidth distance is `abs(f) - omega(k)`. Each
reported center is therefore a frequency magnitude in `[0, 0.5]`, not a
signed frequency.

### API

```matlab
[mode, info] = rvmd(Q, K, Alpha)
[mode, info, restart] = rvmd(Q, K, Alpha, Name, Value, ...)
[mode, info, restart] = rvmd('Restart', restart, Name, Value, ...)
```

Required inputs:

| Input | Contract |
|---|---|
| `Q` | Finite, nonempty numeric matrix, `S × T`; may be real or complex. |
| `K` | Positive integer number of modes. |
| `Alpha` | Nonnegative scalar bandwidth penalty. |

Options:

| Name | Default | Contract |
|---|---:|---|
| `'Weight'` | `1` | Positive scalar or vector of length `S`. It is reshaped to a column and normalized by its mean. |
| `'Tolerance'` | `5e-3` | Nonnegative stopping tolerance for the sum of relative spectral-mode changes. |
| `'MaximumSteps'` | `500` | Positive total iteration cap. On restart this is not an additional-step count. |
| `'InitFreqType'` | `1` | `-1` random, `0` all zero, `1` uniformly distributed. |
| `'InitFreqMaximum'` | `0.5` | Upper initialization frequency, capped at Nyquist `0.5`. |
| `'Device'` | `'cpu'` | `'cpu'` or `'gpu'`. |
| `'FPPrecision'` | `'single'` | `'single'` or `'double'`. |
| `'nDC'` | `0` | Number of leading internal modes whose center remains fixed at zero; must not exceed `K`. |
| `'Display'` | `'off'` | `'off'` or `'iter'`. |
| `'Restart'` | `0` | State returned by a previous call; the restart-only form above is preferred. |

Outputs:

| Field | Meaning |
|---|---|
| `mode.phi` | `S × K` weighted-unit-norm spatial modes. |
| `mode.c` | `T × K` time coefficients. |
| `mode.omega` | `K × 1` final center frequencies, sorted low to high. |
| `mode.energy` | `1 × K`, `sum(abs(mode.c).^2,1)`. |
| `info.Iteration.steps` | Number of completed Gauss--Seidel sweeps. |
| `info.Iteration.omega` | Unsorted internal center-frequency history, including initialization. |
| `info.Iteration.difference` | Convergence metric after each completed sweep. |
| `info.Iteration.converged` | Whether the final metric satisfies `Tolerance`. |
| `restart` | Unsorted internal state used to continue the identical trajectory. |

### The exact filter convention

Both code paths implement

```text
g(f) = 1 / (1 + 2*Alpha*(distance from center)^2)
```

where the distance is `f - omega(k)` for the real nonnegative spectrum and
`abs(f) - omega(k)` for the complex full spectrum. The half-power bandwidth
between the two cutoff points is

```text
bandwidth = sqrt((2*sqrt(2) - 2) / Alpha).
```

This factor of `2` is part of the published RVMD definition. See
[the derivation and validation contract](docs/theory-and-validation.md).

### Weighted RVMD

```matlab
w = cellVolume(:);                 % positive, length S
[mode, info] = rvmd(Q, K, Alpha, 'Weight', w);

% Every column is normalized in the implemented inner product:
wNormalized = w / mean(w);
sqrt(sum(abs(mode.phi).^2 .* wNormalized, 1))
```

### Restart

`MaximumSteps` always means a total cap:

```matlab
[~, ~, state] = rvmd(Q, K, Alpha, 'MaximumSteps', 200);
[mode, info] = rvmd('Restart', state, 'MaximumSteps', 1000);
```

The second call continues steps 201 through at most 1000. Immutable problem
settings (data, `K`, `Alpha`, weights, precision, real/complex branch,
initialization, and `nDC`) come from `state`.
The saved computation device is also reused by default; explicitly pass
`'Device','cpu'` or `'Device','gpu'` to move the continuation. Exact
floating-point trajectory equivalence is guaranteed when the device and
precision are unchanged.

### Canonicalization and non-uniqueness

The implementation removes amplitude ambiguity with weighted unit norm. It
removes the remaining per-mode phase/sign ambiguity by making the largest
entry of each spatial mode real and positive, then sorts outputs by center
frequency.

RVMD is nevertheless a non-convex optimization. Different initializations
or degenerate modes may converge to different stationary decompositions.
The repository provides one canonical implementation and representation;
it does not claim a globally unique optimizer.

### Examples and tests

The four `case*` directories cover a non-stationary signal, the Lorenz
attractor, a transient cylinder wake, and motion-capture data. The
`tutorial_CylinderWake` directory contains additional flow-analysis tools.

Run the MATLAB regression suite from the repository root:

```matlab
results = runtests('tests/test_rvmd.m');
assertSuccess(results)
```

GNU Octave users can run the executable compatibility suite with:

```bash
octave-cli --no-gui --quiet --eval "addpath('tests'); run_octave_tests;"
```

The suite checks the filter coefficient, real and complex synthetic
frequencies, conjugation/phase convention, weighted normalization, zero
data, input validation, and restart equivalence.

## 中文

本仓库只保留一个实现文件 `rvmd.m`，同时正确处理实值与复值时空数据。
输入 `Q` 的尺寸为 `S × T`（空间点 × 时间快照），重构方式为：

```matlab
Q_rec = mode.phi * mode.c.';  % 注意这里是非共轭转置 .'
```

### 实值与复值分支

| 输入 | 空间模态 | 频谱 | 中心频率含义 |
|---|---|---|---|
| 实值 `Q` | 实值 | 非负单边谱，并按 Hermitian 对称重构 | `[0,0.5]` 内的频率 |
| 复值 `Q` | 复值 | 完整双边移位频谱 | `abs(f)` 的能量加权平均，即频率绝对值 |

复值目标函数使用 `(|f|-omega(k))^2`，因此一个中心频率同时描述
`+omega(k)` 与 `-omega(k)` 附近的能量；它不是带符号中心频率。

### 基本调用

```matlab
[mode, info] = rvmd(Q, K, Alpha);
[mode, info, restart] = rvmd(Q, K, Alpha, ...
    'Weight', w, 'Tolerance', 1e-6, ...
    'MaximumSteps', 1000, 'FPPrecision', 'double');
```

严格采用的滤波器为

```text
1 / (1 + 2*Alpha*(距中心频率的距离)^2)
```

不是 `1 + 4*Alpha*(...)^2`。实值分支的距离为 `f-omega(k)`，复值
分支为 `abs(f)-omega(k)`。完整推导、离散 FFT 端点权重和共轭约定见
[理论与验证说明](docs/theory-and-validation.md)。

### 名称-值参数

| 参数 | 默认值 | 说明 |
|---|---:|---|
| `'Weight'` | `1` | 正标量或长度为 `S` 的正向量；内部转为列向量并除以均值。 |
| `'Tolerance'` | `5e-3` | 迭代停止阈值。 |
| `'MaximumSteps'` | `500` | 总迭代步上限；断点续算时也不是“新增步数”。 |
| `'InitFreqType'` | `1` | `-1` 随机、`0` 全零、`1` 均匀分布。 |
| `'InitFreqMaximum'` | `0.5` | 初始频率上界，最高截断到 Nyquist 频率 `0.5`。 |
| `'Device'` | `'cpu'` | `'cpu'` 或 `'gpu'`。GPU 需要 Parallel Computing Toolbox。 |
| `'FPPrecision'` | `'single'` | `'single'` 或 `'double'`。 |
| `'nDC'` | `0` | 固定为零中心频率的前置内部模态数，不得超过 `K`。 |
| `'Display'` | `'off'` | `'off'` 或 `'iter'`。 |

输出模态按中心频率从低到高排序。`mode.energy` 定义为每个时间系数
的离散平方范数。`info.Iteration.omega` 保存的是排序前的内部迭代轨迹。

### 断点续算

```matlab
[~, ~, state] = rvmd(Q, K, Alpha, 'MaximumSteps', 200);
[mode, info] = rvmd('Restart', state, 'MaximumSteps', 1000);
```

第二次调用从第 201 步继续，最多算到总计 1000 步。数据、`K`、
`Alpha`、权重、精度、实/复值分支、初始化方式和 `nDC` 均从保存状态
恢复，避免续算轨迹发生漂移。默认也会沿用保存时的计算设备；只有显式
传入新的 `'Device'` 才会迁移。保持设备和精度不变时，续算保持同一条
浮点轨迹。

### “唯一”的准确含义

仓库中只有 `rvmd.m` 这一份算法实现。代码还固定了每个模态的符号/
常相位：空间模态绝对值最大的元素被规范为正实数。这样消除了输出表达
中的任意相位，但 RVMD 目标函数仍然非凸；不同初始化或简并数据可能
得到不同驻点，不能宣称全局优化解在数学上唯一。

## Citation

If you use RVMD, please cite:

> Liao, Z.-M., Zhao, Z., Chen, L.-B., Wan, Z.-H., Liu, N.-S. & Lu,
> X.-Y. (2023). Reduced-order variational mode decomposition to reveal
> transient and non-stationary dynamics in fluid flows. *Journal of Fluid
> Mechanics*, **966**, A7. https://doi.org/10.1017/jfm.2023.435

```bibtex
@article{liao2023rvmd,
  title   = {Reduced-order variational mode decomposition to reveal transient and non-stationary dynamics in fluid flows},
  author  = {Liao, Z.-M. and Zhao, Z. and Chen, L.-B. and Wan, Z.-H. and Liu, N.-S. and Lu, X.-Y.},
  journal = {Journal of Fluid Mechanics},
  volume  = {966},
  pages   = {A7},
  year    = {2023},
  doi     = {10.1017/jfm.2023.435}
}
```

## License

MIT; see [LICENSE](LICENSE).
