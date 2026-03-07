# Reduced-order Variational Mode Decomposition (RVMD)

[English](#english) | [中文](#中文)

---

<a name="english"></a>

## Overview

**RVMD** is a MATLAB implementation of Reduced-order Variational Mode Decomposition, a method that combines variational mode decomposition (VMD) with low-order modal representation. RVMD can adaptively extract **low-order dynamics featured by transient or non-stationary properties** from space-time data.

Unlike standard VMD which operates on 1D signals, RVMD decomposes a space-time data matrix **Q** (size S × T) into a set of modes, each consisting of a spatial mode **φ**, a time-evolution coefficient **c(t)**, and a characteristic central frequency **ω**.

## Installation

1. Clone or download this repository.
2. Add the repository root to your MATLAB path:
   ```matlab
   addpath('/path/to/rvmd')
   ```

No additional toolboxes are required for CPU computation. GPU computation requires the Parallel Computing Toolbox.

## Quick Start

```matlab
% Q: data matrix of size [S, T] (S = spatial points, T = time steps)
% K: number of modes to extract
% Alpha: filtering parameter (bandwidth control)

[mode, info] = rvmd(Q, K, Alpha);

% Access results
mode.phi    % spatial modes    [S x K]
mode.c      % time-evolution coefficients [T x K]
mode.omega  % central frequencies [K x 1]
mode.energy % mode energy [1 x K]

% Reconstruct the data
Q_rec = mode.phi * mode.c.';
```

## API Reference

```matlab
[mode, info] = rvmd(Q, K, Alpha)
[mode, info] = rvmd(Q, K, Alpha, Name, Value)
[mode, info, restart] = rvmd(Q, K, Alpha, ...)
[mode, info, restart] = rvmd('Restart', restart_data, ...)
```

### Required Arguments

| Argument | Type | Description |
|----------|------|-------------|
| `Q` | matrix [S × T] | Input data matrix. Rows are spatial points, columns are time steps. |
| `K` | positive integer | Number of modes to extract. |
| `Alpha` | positive scalar | Filtering parameter controlling mode bandwidth. Larger values produce narrower bandwidth (more frequency-selective). |

### Optional Name-Value Parameters

| Name | Default | Description |
|------|---------|-------------|
| `'Tolerance'` | `5e-3` | Convergence tolerance for stopping criterion. |
| `'MaximumSteps'` | `500` | Maximum number of iteration steps. |
| `'InitFreqType'` | `1` | Central frequency initialization strategy: `-1` = random, `0` = all zero, `1` = uniformly distributed. |
| `'InitFreqMaximum'` | `0.5` | Upper bound for initial central frequencies (capped at 0.5, the Nyquist frequency). |
| `'Device'` | `'cpu'` | Computation device: `'cpu'` or `'gpu'`. |
| `'FPPrecision'` | `'single'` | Floating-point precision: `'single'` or `'double'`. |
| `'Weight'` | `1` | Weight vector of length S for spatially weighted RVMD. Default is uniform weighting. |
| `'Restart'` | `0` | Restart struct from a previous call for continuing computation. |

### Outputs

| Output | Description |
|--------|-------------|
| `mode.phi` | Spatial modes [S × K], each column is a unit-norm spatial mode. |
| `mode.c` | Time-evolution coefficients [T × K]. |
| `mode.omega` | Central frequencies [K × 1], sorted from low to high. |
| `mode.energy` | Mode energy [1 × K], defined as ‖c_k‖². |
| `info` | Struct containing parameter settings and iteration history (`info.Iteration.difference`, `info.Iteration.omega`, `info.Iteration.steps`). |
| `restart` | (Optional) Struct storing internal state for resuming computation later. |

### Restart Example

For large or difficult problems, you can split the computation into stages:

```matlab
% Run for 200 steps first
[~, ~, restart] = rvmd(Q, K, Alpha, 'MaximumSteps', 200);

% Continue for 8000 more steps from where we left off
[mode, info, ~] = rvmd('Restart', restart, 'MaximumSteps', 8000);
```

### Weighted RVMD Example

Apply spatial weighting to emphasize specific regions:

```matlab
w = ones(S, 1);
w(region_of_interest) = 10;  % higher weight for important regions
[mode, info] = rvmd(Q, K, Alpha, 'Weight', w);
```

## Parameter Tuning Guide

| Parameter | Effect | Guidance |
|-----------|--------|----------|
| **K** (number of modes) | Determines how many components are extracted. | Start small and increase. Too few modes miss dynamics; too many may produce spurious modes. |
| **Alpha** (filtering parameter) | Controls mode bandwidth. Larger Alpha → narrower bandwidth. | Typical range: 500–5000. Use larger values for data with well-separated frequencies. |
| **Tolerance** | Convergence threshold. | Use `1e-3` to `1e-5` for most cases. Tighter tolerance gives more accurate results but requires more iterations. |
| **InitFreqMaximum** | Constrains the initial frequency search range. | Set to slightly above the highest frequency of interest (normalized, max 0.5). |

## Examples

### Case 1: 1D Non-stationary Signal

RVMD reduces to VMD for 1D signals. This example decomposes a composite signal with three frequency components and noise.

```matlab
cd case1_NonstationarySignal
NonstationarySignal
```

### Case 2: Lorenz Attractor

Demonstrates RVMD on chaotic Lorenz system data (3D state space), with restart functionality for long iterative computations.

```matlab
cd case2_LorenzAttractor
LorenzAttractor
```

### Case 3: Transient Cylinder Wake

Applies RVMD to 2D velocity field data of a transient cylinder wake flow, extracting shift modes and vortex shedding modes at different stages. Supports GPU acceleration.

```matlab
cd case3_CylinderWake
CylinderWake
```

### Case 4: Motion Capture

Decomposes human motion capture data into frequency-based components (e.g., low-frequency leg-lifting vs. high-frequency jumping).

```matlab
cd case4_MotionCapture
MotionCapture
```

## File Structure

```
rvmd/
├── rvmd.m                          # RVMD function
├── functionSignatures.json         # MATLAB IDE auto-completion support
├── LICENSE                         # MIT License
├── case1_NonstationarySignal/
│   └── NonstationarySignal.m       # 1D signal decomposition example
├── case2_LorenzAttractor/
│   ├── LorenzAttractor.m           # Lorenz system example (with restart)
│   └── data_lorenz.mat             # Lorenz attractor data
├── case3_CylinderWake/
│   ├── CylinderWake.m              # Cylinder wake example
│   └── data_cylinder.mat           # Cylinder wake velocity field data
└── case4_MotionCapture/
    ├── MotionCapture.m             # Motion capture example
    ├── 49_03.bvh                   # BVH motion capture file
    ├── loadbvh.m                   # BVH file loader
    ├── bvh2snapshot.m              # BVH to snapshot converter
    └── MotionCapture.avi           # Result animation
```

## Citation

If you use RVMD in your research, please cite:

> Liao, Z.-M., Zhao, Z., Chen, L.-B., Wan, Z.-H., Liu, N.-S. & Lu, X.-Y. 2023 Reduced-order variational mode decomposition to reveal transient and non-stationary dynamics in fluid flows. *J. Fluid Mech.*, **966**, A7. [doi:10.1017/jfm.2023.435](https://doi.org/10.1017/jfm.2023.435)

```bibtex
@article{liao2023rvmd,
  title={Reduced-order variational mode decomposition to reveal transient and non-stationary dynamics in fluid flows},
  author={Liao, Z.-M. and Zhao, Z. and Chen, L.-B. and Wan, Z.-H. and Liu, N.-S. and Lu, X.-Y.},
  journal={Journal of Fluid Mechanics},
  volume={966},
  pages={A7},
  year={2023},
  doi={10.1017/jfm.2023.435}
}
```

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

---

<a name="中文"></a>

# 降阶变分模态分解 (RVMD)

## 概述

**RVMD** 是降阶变分模态分解 (Reduced-order Variational Mode Decomposition) 的 MATLAB 实现。该方法将变分模态分解 (VMD) 与低阶模态表示相结合，能够从时空数据中自适应地提取**以瞬态或非定常特性为特征的低阶动力学**。

与仅处理一维信号的标准 VMD 不同，RVMD 将时空数据矩阵 **Q**（尺寸 S × T）分解为一组模态，每个模态由空间模态 **φ**、时间演化系数 **c(t)** 和特征中心频率 **ω** 组成。

## 安装

1. 克隆或下载本仓库。
2. 将仓库根目录添加到 MATLAB 路径：
   ```matlab
   addpath('/path/to/rvmd')
   ```

CPU 计算无需额外工具箱。GPU 计算需要 Parallel Computing Toolbox。

## 快速上手

```matlab
% Q: 数据矩阵，尺寸 [S, T]（S = 空间采样点数，T = 时间步数）
% K: 提取的模态数
% Alpha: 滤波参数（带宽控制）

[mode, info] = rvmd(Q, K, Alpha);

% 获取结果
mode.phi    % 空间模态       [S x K]
mode.c      % 时间演化系数   [T x K]
mode.omega  % 中心频率       [K x 1]
mode.energy % 模态能量       [1 x K]

% 重构数据
Q_rec = mode.phi * mode.c.';
```

## 接口说明

```matlab
[mode, info] = rvmd(Q, K, Alpha)
[mode, info] = rvmd(Q, K, Alpha, Name, Value)
[mode, info, restart] = rvmd(Q, K, Alpha, ...)
[mode, info, restart] = rvmd('Restart', restart_data, ...)
```

### 必需参数

| 参数 | 类型 | 说明 |
|------|------|------|
| `Q` | 矩阵 [S × T] | 输入数据矩阵。行为空间采样点，列为时间步。 |
| `K` | 正整数 | 提取的模态数。 |
| `Alpha` | 正标量 | 滤波参数，控制模态带宽。值越大，带宽越窄（频率选择性越强）。 |

### 可选名称-值参数

| 名称 | 默认值 | 说明 |
|------|--------|------|
| `'Tolerance'` | `5e-3` | 迭代收敛容差。 |
| `'MaximumSteps'` | `500` | 最大迭代步数。 |
| `'InitFreqType'` | `1` | 中心频率初始化策略：`-1` = 随机分布，`0` = 全零，`1` = 均匀分布。 |
| `'InitFreqMaximum'` | `0.5` | 初始中心频率上界（不超过 Nyquist 频率 0.5）。 |
| `'Device'` | `'cpu'` | 计算设备：`'cpu'` 或 `'gpu'`。 |
| `'FPPrecision'` | `'single'` | 浮点精度：`'single'`（单精度）或 `'double'`（双精度）。 |
| `'Weight'` | `1` | 长度为 S 的权重向量，用于空间加权 RVMD。默认为均匀加权。 |
| `'Restart'` | `0` | 上次调用返回的重启结构体，用于继续计算。 |

### 输出

| 输出 | 说明 |
|------|------|
| `mode.phi` | 空间模态 [S × K]，每列为单位范数的空间模态。 |
| `mode.c` | 时间演化系数 [T × K]。 |
| `mode.omega` | 中心频率 [K × 1]，按从低到高排序。 |
| `mode.energy` | 模态能量 [1 × K]，定义为 ‖c_k‖²。 |
| `info` | 包含参数设置和迭代历史的结构体（`info.Iteration.difference`，`info.Iteration.omega`，`info.Iteration.steps`）。 |
| `restart` | （可选）保存内部状态的结构体，用于后续恢复计算。 |

### 重启示例

对于大规模或困难问题，可以分阶段计算：

```matlab
% 先运行 200 步
[~, ~, restart] = rvmd(Q, K, Alpha, 'MaximumSteps', 200);

% 从断点继续运行 8000 步
[mode, info, ~] = rvmd('Restart', restart, 'MaximumSteps', 8000);
```

### 加权 RVMD 示例

通过空间加权突出特定区域：

```matlab
w = ones(S, 1);
w(region_of_interest) = 10;  % 对重要区域赋予更高权重
[mode, info] = rvmd(Q, K, Alpha, 'Weight', w);
```

## 参数调节指南

| 参数 | 作用 | 建议 |
|------|------|------|
| **K**（模态数） | 决定提取多少个分量。 | 从小值开始逐步增加。模态数过少会遗漏动力学特征，过多可能产生伪模态。 |
| **Alpha**（滤波参数） | 控制模态带宽。Alpha 越大，带宽越窄。 | 常用范围：500–5000。数据频率分离度高时使用较大值。 |
| **Tolerance**（容差） | 收敛阈值。 | 一般使用 `1e-3` 到 `1e-5`。容差越小结果越精确，但需要更多迭代。 |
| **InitFreqMaximum** | 限制初始频率搜索范围。 | 设为略高于关注的最高频率（归一化值，最大 0.5）。 |

## 算例

### 算例 1：一维非定常信号

对于一维信号，RVMD 退化为 VMD。本例分解包含三个频率分量和噪声的合成信号。

```matlab
cd case1_NonstationarySignal
NonstationarySignal
```

### 算例 2：Lorenz 吸引子

在混沌 Lorenz 系统数据（三维状态空间）上演示 RVMD，并展示重启功能用于长时间迭代计算。

```matlab
cd case2_LorenzAttractor
LorenzAttractor
```

### 算例 3：瞬态圆柱绕流

将 RVMD 应用于瞬态圆柱绕流的二维速度场数据，提取位移模态和不同阶段的涡脱落模态。支持 GPU 加速。

```matlab
cd case3_CylinderWake
CylinderWake
```

### 算例 4：运动捕捉

将人体运动捕捉数据分解为基于频率的分量（如低频抬腿运动 vs. 高频跳跃运动）。

```matlab
cd case4_MotionCapture
MotionCapture
```

## 文件结构

```
rvmd/
├── rvmd.m                          # RVMD 函数
├── functionSignatures.json         # MATLAB IDE 自动补全支持
├── LICENSE                         # MIT 许可证
├── case1_NonstationarySignal/
│   └── NonstationarySignal.m       # 一维信号分解算例
├── case2_LorenzAttractor/
│   ├── LorenzAttractor.m           # Lorenz 系统算例（含重启）
│   └── data_lorenz.mat             # Lorenz 吸引子数据
├── case3_CylinderWake/
│   ├── CylinderWake.m              # 圆柱绕流算例
│   └── data_cylinder.mat           # 圆柱绕流速度场数据
└── case4_MotionCapture/
    ├── MotionCapture.m             # 运动捕捉算例
    ├── 49_03.bvh                   # BVH 运动捕捉文件
    ├── loadbvh.m                   # BVH 文件加载器
    ├── bvh2snapshot.m              # BVH 转快照工具
    └── MotionCapture.avi           # 结果动画
```

## 引用

如果您在研究中使用了 RVMD，请引用：

> Liao, Z.-M., Zhao, Z., Chen, L.-B., Wan, Z.-H., Liu, N.-S. & Lu, X.-Y. 2023 Reduced-order variational mode decomposition to reveal transient and non-stationary dynamics in fluid flows. *J. Fluid Mech.*, **966**, A7. [doi:10.1017/jfm.2023.435](https://doi.org/10.1017/jfm.2023.435)

```bibtex
@article{liao2023rvmd,
  title={Reduced-order variational mode decomposition to reveal transient and non-stationary dynamics in fluid flows},
  author={Liao, Z.-M. and Zhao, Z. and Chen, L.-B. and Wan, Z.-H. and Liu, N.-S. and Lu, X.-Y.},
  journal={Journal of Fluid Mechanics},
  volume={966},
  pages={A7},
  year={2023},
  doi={10.1017/jfm.2023.435}
}
```

## 许可证

本项目采用 MIT 许可证，详见 [LICENSE](LICENSE)。
