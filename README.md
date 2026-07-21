# Reduced-order Variational Mode Decomposition (RVMD)

[English](#english) · [中文](#中文)

## English

`rvmd.m` provides real- and complex-valued reduced-order variational mode
decomposition in MATLAB.
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

For complex data, each reported center is a frequency magnitude in
`[0, 0.5]`.

### API

```matlab
[mode, info] = rvmd(Q, K, Alpha)
[mode, info, restart] = rvmd(Q, K, Alpha, Name, Value, ...)
[mode, info, restart] = rvmd('Restart', restart, Name, Value, ...)
```

Required inputs:

| Input | Description |
|---|---|
| `Q` | Finite, nonempty numeric matrix, `S × T`; may be real or complex. |
| `K` | Positive integer number of modes. |
| `Alpha` | Nonnegative scalar bandwidth penalty. |

Options:

| Name | Default | Description |
|---|---:|---|
| `'Weight'` | `1` | Positive scalar or vector of length `S`. It is reshaped to a column and normalized by its mean. |
| `'Tolerance'` | `5e-3` | Nonnegative stopping tolerance for the sum of relative spectral-mode changes. |
| `'MaximumSteps'` | `500` | Positive total iteration cap for initial and restarted runs. |
| `'InitFreqType'` | `1` | `-1` random, `0` all zero, `1` uniformly distributed. |
| `'InitFreqMaximum'` | `0.5` | Upper initialization frequency, capped at Nyquist `0.5`. |
| `'Device'` | `'cpu'` | `'cpu'` or `'gpu'`. |
| `'FPPrecision'` | `'single'` | `'single'` or `'double'`. |
| `'nDC'` | `0` | Number of leading internal modes whose center remains fixed at zero; must not exceed `K`. |
| `'Display'` | `'off'` | `'off'` or `'iter'`. |
| `'Restart'` | `0` | State returned by a previous call. |

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

The second call continues from step 200 to a total limit of 1000 steps. The
saved data and settings are reused. Pass `'Device','cpu'` or
`'Device','gpu'` to continue on a different device.

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

The suite covers real and complex inputs, spatial weighting, input
validation, and restart.

## 中文

本仓库提供 MATLAB 实现 `rvmd.m`，支持实值与复值时空数据。
输入 `Q` 的尺寸为 `S × T`（空间点 × 时间快照），重构方式为：

```matlab
Q_rec = mode.phi * mode.c.';  % 注意这里是非共轭转置 .'
```

### 实值与复值分支

| 输入 | 空间模态 | 频谱 | 中心频率含义 |
|---|---|---|---|
| 实值 `Q` | 实值 | 非负单边谱，并按 Hermitian 对称重构 | `[0,0.5]` 内的频率 |
| 复值 `Q` | 复值 | 完整双边移位频谱 | `abs(f)` 的能量加权平均，即频率绝对值 |

复值输入的中心频率表示 `[0,0.5]` 内的频率绝对值。

### 基本调用

```matlab
[mode, info] = rvmd(Q, K, Alpha);
[mode, info, restart] = rvmd(Q, K, Alpha, ...
    'Weight', w, 'Tolerance', 1e-6, ...
    'MaximumSteps', 1000, 'FPPrecision', 'double');
```

### 名称-值参数

| 参数 | 默认值 | 说明 |
|---|---:|---|
| `'Weight'` | `1` | 正标量或长度为 `S` 的正向量；内部转为列向量并除以均值。 |
| `'Tolerance'` | `5e-3` | 迭代停止阈值。 |
| `'MaximumSteps'` | `500` | 初始计算和断点续算的总迭代步上限。 |
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

第二次调用从第 200 步继续，总迭代步上限为 1000。数据和参数从保存
状态恢复；传入新的 `'Device'` 可以更换计算设备。

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
