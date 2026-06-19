# EOSCalculator — 状态方程与熔化计算管线

`EOSCalculator` 是 LAMMPS + DeePMD 的自动化熔化/热力学计算管线，设计模式与 `DPBuilder` 一致（`extrempy/campaign/single_element_dp.py`）。

## 目录

- [设计思路](#设计思路)
- [三阶段工作流](#三阶段工作流)
- [快速开始](#快速开始)
- [API 参考](#api-参考)
- [配置参数](#配置参数)
- [目录结构](#目录结构)
- [后处理](#后处理)
- [批量运行](#批量运行)

---

## 设计思路

```
DPBuilder (DP势构建)          EOSCalculator (EOS/熔化计算)
─────────────────             ─────────────────────────────
generate_init_aimd(segs)      generate_two_phase()
submit_init_aimd()            submit_two_phase()
generate_dpgen(segs)          generate_npt(phases)
submit_dpgen()                submit_npt()
inspect()                     analyze_two_phase() / analyze_npt()
_get_tm()                     _get_tm()
build_all_elements()          run_eos_all()
```

核心原则：

1. **所有配置集中在构造函数**，一目了然
2. **Generate / Submit 分离** — 先生成输入文件检查，再提交
3. **逐元素子类** `ElementEOSCalculator(name, ...)` — 自动查找势函数、POSCAR、熔点
4. **模板驱动** — LAMMPS 输入由 Jinja2 模板渲染，内置 5 个模板（`campaign/templates/`）

---

## 三阶段工作流

### Phase 1: Two-phase — 双相法测熔点

在不同温度下运行固液共存模拟，通过 Q4/Q6 序参量判断最终状态（固/液），确定熔点区间。

模板: `templates/two-phase.j2`

关键步骤：
1. 读入固液界面 POSCAR，沿 z 方向复制（`nz=10`）
2. z 方向切分为上下两组
3. 上半边加热至 `T_superheat` 熔化
4. 冷却至 `Tm_estimate` 平衡
5. 整体 NPT 长时间演化
6. 每 100 步输出 thermo.dat，含上下半区 RDF + Q4/Q6

### Phase 2: NPT — 状态方程计算

在 Tm 以下多个温度点跑 NPT（固体/液体），提取体积、密度、能量随温度的变化。

模板: `templates/npt-solid.j2`, `templates/npt-liquid.j2`

输出: `thermo.dat`（温度、压力、体积、密度、能量、MSD）

### Phase 3: NVT 轨迹 — 动力学计算

在熔点附近跑 NVT 产出轨迹 dump，供后续 VDOS/VACF 分析。

模板: `templates/nvt-solid-traj.j2`, `templates/nvt-liquid-traj.j2`

---

## 快速开始

### 模式 A: 从 DPGEN 项目一键联动（推荐）

如果 DP 模型和 POSCAR 都来自 DPGEN 项目，只需给一个 `dpgen_dir`：

```python
from extrempy import ElementEOSCalculator

calc = ElementEOSCalculator('Al',
    work_root    = '/share/zeng/metals/dpmd',
    dpgen_dir    = '/share/zeng/metals/sample',       # 自动搜 model + POSCAR
    machine_template = '~/template/dpgen-machine.json',
    partition    = 'gpu_share',                       # slurm 分区
    nodes        = 1,                                 # 覆盖 machine_template
    ntasks_per_node = 8,
    wall_time    = '48:00:00',
    gres         = 'gpu:1',
)

# 分步执行
calc.generate_two_phase()           # 生成双相法输入
calc.submit_two_phase()             # 提交到集群
calc.analyze_two_phase()            # 分析结果 → Tm 区间

calc.generate_npt()                 # 生成 NPT 输入（固体+液体）
calc.submit_npt()

calc.generate_nvt_traj()            # 生成 NVT 轨迹输入
calc.submit_nvt_traj()

# 或一步到位
calc.run_all(submit=True)
```

### 模式 B: 显式指定所有路径

如果模型是额外训练的、POSCAR 是额外准备的：

```python
calc = ElementEOSCalculator('Al',
    work_root    = '/share/zeng/metals/dpmd',
    dp_model_path = '/extra/train/frozen_model.pb',   # 额外训练的模型
    poscar_path   = '/extra/confs/Al-fcc.POSCAR',      # 额外准备的 POSCAR
    machine_template = '~/template/dpgen-machine.json',
    partition = 'cpu', nodes=2, ntasks_per_node=32,
)
```

### 模式 C: 混合（DP 模型走 dpgen_dir，POSCAR 走外部目录）

```python
calc = ElementEOSCalculator('Al',
    work_root = '/share/zeng/metals/dpmd',
    dpgen_dir = '/share/zeng/metals/sample',            # 从这里找 DP 模型
    poscar_dir = '/share/zeng/metals/poscar',            # 从这里找 POSCAR
    machine_template = '~/template/dpgen-machine.json',
)
```

### 仅生成不提交

```python
calc.run_all(submit=False)   # 只生成所有输入文件，不提交
```

### 单独分析 NPT 结果

```python
summary = calc.analyze_npt()
print(summary)

from extrempy.md.thermo import plot_thermo_summary
plot_thermo_summary(summary, 'Al')
```

---

## API 参考

### `EOSCalculator`

| 方法 | 说明 |
|---|---|
| `generate_two_phase()` | 在 3 个候选温度下生成双相法输入 |
| `submit_two_phase(submit=True)` | 提交双相法任务 |
| `analyze_two_phase()` → dict | 读取 dump → Q4/Q6 → 返回每个温度点的固/液判定 |
| `generate_npt(phases=('solid','liquid'))` | 在 N 个温度点生成 NPT 输入（支持指定相） |
| `submit_npt(submit=True)` | 提交 NPT 任务 |
| `analyze_npt()` → DataFrame | 遍历 NPT 目录 → 读取 thermo.dat → MSD 相判定 → 汇总 |
| `generate_nvt_traj(phases=('solid','liquid'))` | 生成 NVT 轨迹输入 |
| `submit_nvt_traj(submit=True)` | 提交 NVT 轨迹任务 |
| `run_all(submit=True)` | 按顺序执行上述所有 generate + submit |

### Hooks（可被子类覆盖）

| 方法 | 默认行为 |
|---|---|
| `_get_tm()` → float | 从 `ELEMENT_PHASE_DATA` 查熔点 |
| `_find_pot()` → str | `dp_model_path > dpgen_dir/{el}_sample/iter.*/00.train/000/` |
| `_find_poscar(idx=0)` → str | `poscar_path > poscar_dir/{el}-*POSCAR > dpgen_dir/{el}/confs/*.POSCAR` |
| `_get_two_phase_temps()` → list | `Tm + [-ΔT+shift, 0+shift, ΔT+shift]` |
| `_get_npt_temps()` → list | `Tm + [-(n//2)..(n//2)]*dT + shift` |

### `ElementEOSCalculator(element, **kwargs)`

绑定到单个元素的子类，自动实现上述 hooks。

### `run_eos_all(elements, work_root, **kwargs)`

```python
results = run_eos_all(['Al', 'Cu', 'Au'],
    work_root='/share/zeng/metals/dpmd',
    dpgen_dir='/share/zeng/metals/sample',
    machine_template='~/template/dpgen-machine.json')
# → {'Al': 'generated', 'Cu': 'generated', 'Au': 'generated'}
```

---

## 配置参数

### 模型与结构文件

| 参数 | 默认值 | 说明 |
|---|---|---|
| `dpgen_dir` | None | DPGEN 项目根目录。自动搜 model + POSCAR |
| `dp_model_path` | None | 显式 DP frozen_model 路径（跳过搜索） |
| `poscar_path` | None | 显式 POSCAR 路径（跳过搜索） |
| `poscar_dir` | None | POSCAR 目录，glob `{element}-*POSCAR` |

搜索优先级：`dp_model_path > dpgen_dir`，`poscar_path > poscar_dir > dpgen_dir/confs`

### Slurm 提交

| 参数 | 默认值 | 说明 |
|---|---|---|
| `machine_template` | None | dpgen `machine.json` 路径，用作默认资源配置 |
| `partition` | None | Slurm 分区 |
| `nodes` | 1 | 节点数 |
| `ntasks_per_node` | 32 | 每节点任务数 |
| `wall_time` | '24:00:00' | 最长运行时间 |
| `gres` | None | GPU 资源 (e.g. `'gpu:1'`) |
| `lmp_command` | `'lmp -in run.in > log.run'` | LAMMPS 运行命令 |

优先级：构造函数参数 > `machine_template(model_devi/fp)` 段

### 温度与计算

| 参数 | 默认值 | 说明 |
|---|---|---|
| `two_phase_frac` | 0.10 | 双相法温度偏移 = Tm × `two_phase_frac` |
| `two_phase_shift` | 300 | 双相法整体偏移 (K) |
| `two_phase_nz` | 10 | 双相法 z 方向复制数 |
| `npt_n` | 5 | NPT 温度点数 |
| `npt_dT` | 100 | NPT 温度间隔 (K) |
| `npt_shift` | -600 | NPT 相对 Tm 的整体偏移 (K) |
| `supercell` | (5,5,5) | 晶胞复制数 |
| `liquid_superheat` | 1.9 | 液相过热倍数 (`T_high = Tm × liquid_superheat`) |
| `equil_steps` | 100000 | 平衡步数 |
| `heat_steps` | 10000 | 加热步数 |
| `dump_freq` | 10 | NVT 轨迹 dump 频率 |
| `dt` | 0.001 | LAMMPS 时间步长 (ps) |
| `Q_cutoff` | 3.0 | Q 序参量截断半径 (Å) |
| `pressure` | 0.0001 | 压强 (万 bar) |

> **关于温度计算示例**（以 Al 为例，Tm = 933 K）：
> - Two-phase 温度: `[933-100+300, 933+0+300, 933+100+300]` = [1133, 1233, 1333] K
> - NPT 温度: `933 + [-200, -100, 0, 100, 200] - 600` = 过滤 ≥250 → [333, 433, 533] K

---

## 目录结构

```
{work_root}/{element}/
├── melt/
│   ├── {Tm1}k/
│   │   ├── run.in          # LAMMPS 输入 (two-phase.j2 渲染)
│   │   ├── confs.data      # LAMMPS 结构文件
│   │   ├── cp.pb           # DeePMD 势函数
│   │   ├── job.sbatch      # Slurm 脚本
│   │   └── traj/           # dump 输出
│   ├── {Tm2}k/
│   └── {Tm3}k/
├── npt/
│   ├── {T1}k_solid/
│   │   ├── run.in          # npt-solid.j2
│   │   ├── thermo.dat      # 热力学输出
│   │   └── rdf.txt         # RDF
│   ├── {T1}k_liquid/
│   │   ├── run.in          # npt-liquid.j2
│   │   └── ...
│   ├── {T2}k_solid/
│   └── ...
└── traj/
    ├── {Tm}k_solid/
    │   ├── run.in          # nvt-solid-traj.j2
    │   └── traj/           # NVT dump 轨迹
    └── {Tm}k_liquid/
```

---

## 后处理

### 双相法分析

```python
result = calc.analyze_two_phase()
# {'results': {1133: 'solid', 1233: 'partial', 1333: 'liquid'},
#  'Tm_interval': (1233, 1333)}
```

内部使用 Q4/Q6 序参量判断每个 dump 的相（阈值：Q4 > 0.1 & Q6 > 0.3 → 固体）。

### NPT 热力学分析

```python
df = calc.analyze_npt()      # pandas DataFrame
plot_thermo_summary(df, 'Al')  # 温度 vs 体积/密度/能量 三面板图
```

### RDF 文件读取

```python
from extrempy.md.traj import read_rdf_file
r, g_r = read_rdf_file('rdf.txt')
```

### 批量双相法分析（已有 dump 结果时）

```python
from extrempy.md.traj import batch_analyze_two_phase
results = batch_analyze_two_phase('/path/to/dump/parent/dir')
```

---

## 内置模板说明

5 个 Jinja2 模板位于 `extrempy/campaign/templates/`，与包一同分发：

| 模板 | 用途 | 关键参数 |
|---|---|---|
| `two-phase.j2` | 固液双相法测熔点 | `Tm_estimate`, `T_superheat`, `nx/ny/nz`, `Q_cutoff` |
| `npt-solid.j2` | 固体 NPT 等温等压 | `temperature`, `nx/ny/nz` |
| `npt-liquid.j2` | 液体 NPT（先熔化后降温） | `temperature`, `high_temperature` |
| `nvt-solid-traj.j2` | 固体 NVT 轨迹产出 | `temperature`, `dump_freq` |
| `nvt-liquid-traj.j2` | 液体 NVT 轨迹产出 | `temperature`, `high_temperature`, `dump_freq` |

若需自定义模板，传入自定义路径：

```python
calc = EOSCalculator(
    work_root=...,
    template_dir='/path/to/your/templates',  # 不使用内置模板
    ...)
```

注意：npt 模板中使用 `fparam ${T}` 传递温度给 DeePMD 势函数（温度相关的势函数需要此参数）。

---

## 与 notebook 工作流的对应

EOSCalculator 封装的是 `run_melt.ipynb` 的三个主要提交循环：

| Notebook Cell | EOSCalculator 方法 |
|---|---|
| Cell 4 (two-phase submit) | `generate_two_phase()` + `submit_two_phase()` |
| Cell 6 (NPT submit) | `generate_npt()` + `submit_npt()` |
| Cell 8 (NVT traj submit) | `generate_nvt_traj()` + `submit_nvt_traj()` |
| Cell 11-12 (dump 分析) | `analyze_two_phase()` |
| Cell 19 (thermo 分析) | `analyze_npt()` |
