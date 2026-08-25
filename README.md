# Reduced-order Variational Mode Decomposition (RVMD)

[English](#english) · [中文](#中文)

## English

`rvmd.m` provides real- and complex-valued reduced-order variational mode
decomposition in MATLAB. The `rvmdpy` package provides the matching PyTorch
implementation for CPU, one CUDA GPU, and one spatially sharded multi-GPU run.
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
| `'InitialFrequencies'` | `[]` | Optional `K × 1` initial center frequencies in `[0,0.5]`; the first `nDC` entries must be zero. |
| `'Display'` | `'off'` | `'off'`, `'final'`, or `'iter'`. |
| `'DisplayInterval'` | `20` | Positive command-window display interval. |
| `'OutputFcn'` | `[]` | Callback `stop = fcn(progress,phase)` called at safe sweep boundaries. |
| `'TimeLimit'` | `Inf` | Nonnegative per-call wall-time limit, checked at safe sweep boundaries. |
| `'CheckpointFile'` | `''` | MAT-file for periodic and final checkpoints; empty disables checkpointing. |
| `'CheckpointInterval'` | `50` | Positive checkpoint interval in completed sweeps. |

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
| `info.ExitFlag` | `1` converged, `0` maximum steps, `-1` callback stop, `-2` time limit. |
| `info.StopReason` | Text termination reason. |
| `info.ElapsedTime` | Per-call elapsed wall time in seconds. |
| `restart` | Unsorted internal state used to continue the identical trajectory. |

### Weighted RVMD

```matlab
w = cellVolume(:);                 % positive, length S
[mode, info] = rvmd(Q, K, Alpha, 'Weight', w);

% Every column is normalized in the implemented inner product:
wNormalized = w / mean(w);
sqrt(sum(abs(mode.phi).^2 .* wNormalized, 1))
```

### Monitor and stop safely

`OutputFcn` runs only after a complete Gauss--Seidel sweep, so requesting a
stop always returns a consistent restart state. The callback phases are
`'init'`, `'iter'`, and `'done'`. Progress contains `step`, `difference`,
`omega`, `elapsedTime`, and (during `'done'`) `stopReason`.

```matlab
function stop = monitor(progress, phase)
    stop = false;
    if strcmp(phase, 'iter')
        fprintf('%4d  %.3e  %s\n', progress.step, progress.difference, ...
            mat2str(sort(progress.omega), 5));
        stop = isfile('STOP_RVMD');
    end
end

[mode, info, state] = rvmd(Q, K, Alpha, ...
    'OutputFcn', @monitor, 'TimeLimit', 3600);
```

`Ctrl+C` still aborts MATLAB immediately. For a graceful interactive stop,
use an output callback, a GUI stop flag, or a sentinel file as above.

### Checkpoint and restart

`MaximumSteps` always means a total cap:

```matlab
[~, ~, state] = rvmd(Q, K, Alpha, 'MaximumSteps', 200);
[mode, info] = rvmd('Restart', state, 'MaximumSteps', 1000);
```

The second call continues from step 200 to a total limit of 1000 steps. The
saved data and settings are reused. Pass `'Device','cpu'` or
`'Device','gpu'` to continue on a different device.

For long runs, RVMD can write checkpoints automatically:

```matlab
[mode, info] = rvmd(Q, K, Alpha, ...
    'MaximumSteps', 1000, ...
    'CheckpointFile', 'rvmd_checkpoint.mat', ...
    'CheckpointInterval', 50);

saved = load('rvmd_checkpoint.mat', 'restart');
[mode, info] = rvmd('Restart', saved.restart, 'MaximumSteps', 1500);
```

The previous valid file is retained as `rvmd_checkpoint.mat.prev`. Restart
version 4 stores the exact iterative state and spectrum norm without copying
the original `Q`, and skips recomputing its FFT. Older repository restart
states remain readable. The legacy name-value form
`rvmd(Q,K,Alpha,'Restart',state)` is rejected to prevent silently ignoring
new problem inputs.

### Convergence and Hilbert spectral analysis

```matlab
fs = 4;                         % samples per physical-time unit
[mode, info] = rvmd(Q, 4, 1000);

rvmdplot(info, 'SampleRate', fs);
hs = rvmdhilbert(mode, fs, ...
    'MirrorExtension', true, 'FrequencyBins', 256);
rvmdhilbertplot(hs);
```

`rvmdhilbert` requires no Signal Processing Toolbox. For real coefficients it
constructs Fourier-domain analytic signals, optionally with mirror extension;
for complex coefficients it uses them directly. It returns instantaneous
amplitudes, phases, signed instantaneous frequencies, instantaneous energies,
a memory-efficient sparse Hilbert energy spectrum, and the marginal spectrum.
Low-amplitude samples are masked because their instantaneous frequency is not
meaningful. For a complex coefficient, instantaneous frequency is obtained
from its phase directly; interpretation requires the mode to remain locally
monocomponent.

### Python: CPU and one GPU

The `rvmdpy` package uses one PyTorch tensor implementation for CPU and CUDA:

```bash
python -m pip install -e '.[plot]'
```

For CUDA, install the PyTorch build that matches your system. Check the active
device before a long run:

```bash
python -c "import torch; print(torch.__version__, torch.cuda.is_available()); print(torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'CPU')"
```

For most jobs, use one CUDA GPU. Use CPU when CUDA is unavailable. Multi-GPU is
for one decomposition that is too large for a single GPU.

```python
from rvmdpy import rvmd, rvmdhilbert, rvmdhilbertplot, rvmdplot

mode, info, state = rvmd(
    q, 4, 1000,
    device="cuda",                 # or "cpu"
    fp_precision="double",
    initial_frequencies=[0.02, 0.08, 0.16, 0.30],
    display="iter",
    display_interval=20,
    checkpoint_file="rvmd_restart.pt",
    checkpoint_interval=50,
)

mode, info, state = rvmd(
    restart="rvmd_restart.pt", maximum_steps=1500, device="cuda"
)

figure = rvmdplot(info, sample_rate=4)
analysis = rvmdhilbert(mode, 4, mirror_extension=True)
hilbert_figure = rvmdhilbertplot(analysis)
```

Python options use the MATLAB defaults and snake-case names. `Mode`, `Info`,
`Progress`, and `RestartState` are typed dataclasses containing CPU PyTorch
tensors. Python mode indices are zero-based. Python checkpoints are native
PyTorch checkpoints and are not MATLAB `.mat` restart files.

`q` can be a NumPy array, a PyTorch tensor, or another array-like object. Its
shape is `S × T`: rows are spatial points and columns are snapshots. Frequencies
are stored in cycles per sample; multiply by the sample rate for physical units:

```python
physical_frequency = mode.omega.numpy() * sample_rate
q_reconstructed = mode.phi @ mode.c.T  # T is deliberately nonconjugating
```

The principal Python options are:

| Option | Default | Meaning |
|---|---:|---|
| `weight` | `1` | Positive spatial weights. |
| `tolerance` | `5e-3` | Convergence tolerance. |
| `maximum_steps` | `500` | Total iteration limit, including completed restart steps. |
| `init_freq_type` | `1` | `-1` random, `0` all zero, `1` uniform. |
| `init_freq_maximum` | `0.5` | Upper initial frequency in cycles per sample. |
| `initial_frequencies` | `None` | Explicit length-`K` initial frequencies. |
| `device` | `"cpu"` | `"cpu"`, `"cuda"`, or a specific CUDA device such as `"cuda:1"`. |
| `fp_precision` | `"single"` | `"single"` or `"double"`. |
| `n_dc` | `0` | Number of center frequencies fixed at zero. |
| `display` | `"off"` | `"off"`, `"final"`, or `"iter"`. |
| `output_fcn` | `None` | Called after each iteration; return `True` to stop. |
| `time_limit` | infinity | Time limit for one call. |
| `checkpoint_file` | `None` | `.pt` checkpoint path. |
| `checkpoint_interval` | `50` | Iterations between checkpoints. |

`mode` contains `phi`, `c`, `omega`, and `energy`. `info` contains the settings,
convergence history, elapsed time, and stop reason. The third return value is
used to resume the run.

#### Run the JFM 2023 cylinder case

This example runs RVMD, Hilbert spectral analysis, plots, and checkpointing:

```bash
python examples_python/cylinder_jfm2023.py --device auto \
  --output rvmd_jfm2023_output
```

Results are written to that directory: `mode.pt`, `summary.json`, Hilbert data,
plots, and restart files. Add `--no-plots` if Matplotlib is not installed.

It uses `K=10`, `alpha=1000`, `tolerance=0.002`, sample rate `4`, and initial
frequencies `linspace(0,0.15,10)`. To stop safely, create
`rvmd_jfm2023_output/STOP_RVMD`. Delete that file, then resume with:

```bash
python examples_python/cylinder_jfm2023.py \
  --restart rvmd_jfm2023_output/restart.pt \
  --output rvmd_jfm2023_output --maximum-steps 1000
```

#### Monitor and stop safely from Python

```python
from pathlib import Path

def monitor(progress, phase):
    if phase == "iter":
        print(progress.step, progress.difference, sorted(progress.omega))
    return phase in {"init", "iter"} and Path("STOP_RVMD").is_file()

mode, info, state = rvmd(
    q, K, alpha,
    output_fcn=monitor,
    time_limit=3600,
    checkpoint_file="rvmd_restart.pt",
    checkpoint_interval=50,
)
```

RVMD checks stop conditions and writes checkpoints after each full iteration,
so saved states can be resumed. `Ctrl+C` stops immediately; use checkpoints for
long jobs.

### Python: one decomposition across multiple GPUs

`rvmd_distributed` shards the spatial dimension `S` across one process per GPU.
Each rank stores its local `Q` rows, spatial modes, residual, and weights. The
time coefficients and frequency history are replicated. The required spatial
reductions use NCCL; the Gauss--Seidel mode order is unchanged.
All ranks must receive the same numerical options and checkpoint path.

The ready-to-run `.npy` driver loads only each rank's rows:

```bash
torchrun --standalone --nproc-per-node=4 \
  examples_python/multigpu_rvmd.py snapshots.npy \
  --output rvmd_output --k 8 --alpha 1000 \
  --maximum-steps 1000 \
  --checkpoint rvmd_checkpoint
```

Restart, including with a different number of ranks:

```bash
torchrun --standalone --nproc-per-node=8 \
  examples_python/multigpu_rvmd.py \
  --restart rvmd_checkpoint --output rvmd_output_resumed \
  --maximum-steps 1500 --checkpoint rvmd_checkpoint_8gpu
```

Distributed checkpoints are directories containing independent spatial shard
files and one global state file. They avoid gathering the full residual onto
rank 0 and retain the previous directory with a `.prev` suffix. The checkpoint
path must be visible to every rank; use a shared filesystem for multi-node
runs. Returned `mode.phi` is the local spatial shard; `mode.c` and
`mode.omega` are available on every rank. Use `gather_mode` only when the full
spatial result fits on the destination rank.

Multi-GPU scaling is intended for a single large decomposition whose spatial
matrix operations dominate collective communication. Small problems generally
run faster on one device.

The driver writes one spatial shard per rank and one shared temporal file.
Combine them with:

```python
from pathlib import Path
import torch

output = Path("rvmd_output")
phi = torch.cat(
    [torch.load(path, weights_only=True)["phi"]
     for path in sorted(output.glob("phi_rank_*.pt"))],
    dim=0,
)
temporal = torch.load(output / "temporal.pt", weights_only=True)
q_reconstructed = phi @ temporal["c"].T
```

### Examples and tests

The four `case*` directories cover a non-stationary signal, the Lorenz
attractor, a transient cylinder wake, and motion-capture data. The
`tutorial_CylinderWake` directory contains additional flow-analysis tools.

Run the MATLAB regression suite from the repository root:

```matlab
results = runtests({'tests/test_rvmd.m', 'tests/test_rvmdhilbert.m'});
assertSuccess(results)
```

GNU Octave users can run the executable compatibility suite with:

```bash
octave-cli --no-gui --quiet --eval "addpath('tests'); run_octave_tests;"
```

Run the Python suite, including the two-process spatial-sharding test, with:

```bash
python -m pytest tests_python
```

The suite covers real and complex inputs, spatial weighting, validation,
graceful stopping, checkpoint/restart, and Hilbert spectral analysis.

## 中文

本仓库提供 MATLAB 实现 `rvmd.m` 和 PyTorch 实现 `rvmdpy`，支持实值与
复值时空数据；Python 版可运行于 CPU、单张 CUDA GPU 或空间分片的多 GPU。
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
| `'InitialFrequencies'` | `[]` | 可选的 `K × 1` 初始中心频率，范围 `[0,0.5]`；前 `nDC` 项必须为零。 |
| `'Display'` | `'off'` | `'off'`、`'final'` 或 `'iter'`。 |
| `'DisplayInterval'` | `20` | 命令行迭代信息的显示间隔。 |
| `'OutputFcn'` | `[]` | 在完整 sweep 边界调用 `stop = fcn(progress,phase)`。 |
| `'TimeLimit'` | `Inf` | 本次调用的运行时限，在安全的 sweep 边界检查。 |
| `'CheckpointFile'` | `''` | 自动 checkpoint 文件；空字符串表示关闭。 |
| `'CheckpointInterval'` | `50` | checkpoint 的 sweep 间隔。 |

输出模态按中心频率从低到高排序。`mode.energy` 定义为每个时间系数
的离散平方范数。`info.Iteration.omega` 保存的是排序前的内部迭代轨迹。

### 安全停止与运行中监看

`OutputFcn` 只在完整 Gauss--Seidel sweep 后调用，因此回调请求停止时，
函数会正常返回一致的 `mode`、`info` 和 `restart`。`progress` 包含当前步数、
收敛差、未排序中心频率和运行时间；`phase` 为 `'init'`、`'iter'` 或
`'done'`。可通过回调读取 GUI 状态或 `STOP_RVMD` 标志文件，实现安全停止。

终止原因通过 `info.ExitFlag`、`info.StopReason`、`info.Message` 和
`info.ElapsedTime` 返回。`Display='iter'` 会按 `DisplayInterval` 输出收敛差
和排序后的中心频率。

### 自动 checkpoint 与断点续算

```matlab
[~, ~, state] = rvmd(Q, K, Alpha, 'MaximumSteps', 200);
[mode, info] = rvmd('Restart', state, 'MaximumSteps', 1000);
```

第二次调用从第 200 步继续，总迭代步上限为 1000。数据和参数从保存
状态恢复；传入新的 `'Device'` 可以更换计算设备。

长任务可直接设置：

```matlab
[mode, info] = rvmd(Q, K, Alpha, ...
    'CheckpointFile', 'rvmd_checkpoint.mat', ...
    'CheckpointInterval', 50, 'TimeLimit', 3600);
saved = load('rvmd_checkpoint.mat', 'restart');
[mode, info] = rvmd('Restart', saved.restart, 'MaximumSteps', 1500);
```

写入时先生成临时文件，并保留上一份有效状态为 `.prev`。v4 checkpoint
不再复制原始 `Q`，续算时也不再重复执行输入数据 FFT。续算统一使用
`rvmd('Restart',state,...)`，避免新的 `Q/K/Alpha` 被误认为有效输入。

### Hilbert 谱后处理

```matlab
fs = 4;
hs = rvmdhilbert(mode, fs, 'MirrorExtension', true);
rvmdplot(info, 'SampleRate', fs);
rvmdhilbertplot(hs);
```

`rvmdhilbert` 输出解析时间系数、瞬时幅值、相位、瞬时频率、瞬时能量、
稀疏 Hilbert 能量谱和边际谱。实值系数采用无工具箱依赖的 FFT 解析信号；
复值系数直接按复解析信号处理。相对幅值过低的位置自动屏蔽瞬时频率，
避免把相位噪声解释成物理频率。

### Python CPU 与单 GPU

```bash
python -m pip install -e '.[plot]'
```

需要 CUDA 时，请先安装与本机环境匹配的 PyTorch。长任务开始前可以检查一下
实际使用的设备：

```bash
python -c "import torch; print(torch.__version__, torch.cuda.is_available()); print(torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'CPU')"
```

一般直接用单张 CUDA GPU；没有 CUDA 时使用 CPU。只有一个分解大到单卡放不下
或算不动时，才需要多 GPU。

```python
from rvmdpy import rvmd, rvmdhilbert, rvmdhilbertplot, rvmdplot

mode, info, state = rvmd(
    q, 4, 1000,
    device="cuda",                 # CPU 使用 "cpu"
    fp_precision="double",
    initial_frequencies=[0.02, 0.08, 0.16, 0.30],
    display="iter",
    checkpoint_file="rvmd_restart.pt",
    checkpoint_interval=50,
)

mode, info, state = rvmd(
    restart="rvmd_restart.pt", maximum_steps=1500, device="cuda"
)
analysis = rvmdhilbert(mode, 4, mirror_extension=True)
rvmdplot(info, sample_rate=4)
rvmdhilbertplot(analysis)
```

Python 版与 MATLAB 版采用相同默认值、停止边界、中心频率定义和累计
`maximum_steps` 语义。参数名使用 snake_case；结果为带类型的 dataclass，
其中张量返回到 CPU。Python 模态索引从零开始。Python checkpoint 是原生
PyTorch 格式，不与 MATLAB `.mat` restart 文件互读。

`q` 可以是 NumPy 数组、PyTorch 张量或其他 array-like 对象，形状为
`S × T`：行是空间点，列是时间快照。中心频率的单位是 cycles/sample，乘以
采样率就是物理频率：

```python
physical_frequency = mode.omega.numpy() * sample_rate
q_reconstructed = mode.phi @ mode.c.T  # T 是非共轭转置
```

主要 Python 参数如下：

| 参数 | 默认值 | 说明 |
|---|---:|---|
| `weight` | `1` | 正的空间权重。 |
| `tolerance` | `5e-3` | 收敛阈值。 |
| `maximum_steps` | `500` | 总迭代上限，续算时包含已经完成的步数。 |
| `init_freq_type` | `1` | `-1` 随机、`0` 全零、`1` 均匀分布。 |
| `init_freq_maximum` | `0.5` | 初始频率上界，单位为 cycles/sample。 |
| `initial_frequencies` | `None` | 指定长度为 `K` 的初始频率。 |
| `device` | `"cpu"` | `"cpu"`、`"cuda"` 或 `"cuda:1"` 等具体设备。 |
| `fp_precision` | `"single"` | `"single"` 或 `"double"`。 |
| `n_dc` | `0` | 固定为零的中心频率数量。 |
| `display` | `"off"` | `"off"`、`"final"` 或 `"iter"`。 |
| `output_fcn` | `None` | 每轮迭代后调用；返回 `True` 时停止。 |
| `time_limit` | 无穷 | 单次运行的时间上限。 |
| `checkpoint_file` | `None` | `.pt` 断点文件。 |
| `checkpoint_interval` | `50` | 每隔多少轮保存一次。 |

`mode` 包含 `phi`、`c`、`omega` 和 `energy`。`info` 记录参数、收敛过程、
耗时和停止原因。第三个返回值用于续算。

#### 运行 JFM 2023 圆柱尾流算例

下面的命令会完成 RVMD、Hilbert 谱分析、绘图和断点保存：

```bash
python examples_python/cylinder_jfm2023.py --device auto \
  --output rvmd_jfm2023_output
```

结果会保存在该目录，包括 `mode.pt`、`summary.json`、Hilbert 数据、图片和
断点文件。没有安装 Matplotlib 时可以加 `--no-plots`。

参数是 `K=10`、`alpha=1000`、`tolerance=0.002`、采样率 `4`，初始频率为
`linspace(0,0.15,10)`。需要停止时，创建
`rvmd_jfm2023_output/STOP_RVMD`；程序会在当前这轮结束后保存并退出。删除该文件
后，再运行续算命令：

```bash
python examples_python/cylinder_jfm2023.py \
  --restart rvmd_jfm2023_output/restart.pt \
  --output rvmd_jfm2023_output --maximum-steps 1000
```

#### Python 运行中监看与安全停止

```python
from pathlib import Path

def monitor(progress, phase):
    if phase == "iter":
        print(progress.step, progress.difference, sorted(progress.omega))
    return phase in {"init", "iter"} and Path("STOP_RVMD").is_file()

mode, info, state = rvmd(
    q, K, alpha,
    output_fcn=monitor,
    time_limit=3600,
    checkpoint_file="rvmd_restart.pt",
    checkpoint_interval=50,
)
```

RVMD 每轮结束后检查停止条件并写入 checkpoint，因此保存的状态可以直接续算。
`Ctrl+C` 会立即中断，长任务建议开启 checkpoint。

### Python 单任务多 GPU

多 GPU 版本沿空间维 `S` 分片同一个分解，每个进程独占一张 GPU。局部保存
输入行、空间模态、残差和权重；时间系数与频率历史在各 rank 复制；空间归约
通过 NCCL 完成，Gauss--Seidel 模态更新顺序保持不变。所有 rank 必须使用
相同的数值选项和 checkpoint 路径。

```bash
torchrun --standalone --nproc-per-node=4 \
  examples_python/multigpu_rvmd.py snapshots.npy \
  --output rvmd_output --k 8 --alpha 1000 \
  --maximum-steps 1000 --checkpoint rvmd_checkpoint
```

分布式 checkpoint 是分片目录，不会把完整残差聚集到 rank 0，并保留上一份
`.prev` 目录。checkpoint 路径必须对所有 rank 可见；多节点运行应使用共享
文件系统。恢复时可以改变 GPU 数量，加载器会按新 world size 重新划分
空间分片。返回的 `mode.phi` 是本 rank 的局部行；`mode.c` 和 `mode.omega`
在各 rank 都可用。只有确认完整空间模态能放入目标 rank 内存时才应调用
`gather_mode`。

多 GPU 面向空间矩阵运算占主导的超大任务。小任务的 collective 通信开销
可能高于并行收益。

程序会为每个 rank 保存一份空间模态，并单独保存公共的时间结果。组合方法：

```python
from pathlib import Path
import torch

output = Path("rvmd_output")
phi = torch.cat(
    [torch.load(path, weights_only=True)["phi"]
     for path in sorted(output.glob("phi_rank_*.pt"))],
    dim=0,
)
temporal = torch.load(output / "temporal.pt", weights_only=True)
q_reconstructed = phi @ temporal["c"].T
```

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
