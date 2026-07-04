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
    dpgen_dir    = '/share/zeng/metals/sample/Al/dpgen',  # element-internal
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
    dpgen_dir = '/share/zeng/metals/sample/Al/dpgen',   # element-internal
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
| `run_two_phase(submit=True)` | 仅 Phase 1：generate + submit two-phase |
| `run_property_scans(submit=True)` | 仅 Phase 2+3：generate + submit NPT + NVT（用 `Tm_refined` 如有） |
| `run_all(submit=True)` | 便捷：two-phase + property scans（**用估算 Tm**，不含 analyze） |

> **两步工作流**（推荐，用精修 Tm）：
> ```python
> calc.run_two_phase(submit=True)
> # ... 等 job 跑完 ...
> calc.analyze_two_phase()            # 设置 calc.Tm_refined
> calc.run_property_scans(submit=True)  # NPT/NVT 围绕精修 Tm
> ```
> `run_all()` 仍可用，但用估算 Tm，不包含两相法精修。

### Hooks（可被子类覆盖）

| 方法 | 默认行为 |
|---|---|
| `_get_tm()` → float | 从 `ELEMENT_PHASE_DATA` 查熔点 |
| `_find_pot()` → str | `dp_model_path > dpgen_dir/frozen_model*.pb symlink > dpgen_dir/iter.*/00.train/000/` |
| `_find_poscar(role)` → str | solid_rt: `poscar_path > poscar_dir/{el}-{rt}.POSCAR > dpgen_dir/../confs/ > ASE auto-gen`; liquid: 同前按 `-LIQ` label，找不到静默回退 solid_rt |
| `_get_two_phase_temps()` → list | `Tm + [-ΔT+shift, 0+shift, ΔT+shift]` |
| `_get_npt_temps()` → list | `Tm + [-(n//2)..(n//2)]*dT + shift`（默认 shift=0，以 Tm 为中心） |

### `ElementEOSCalculator(element, structure=None, phase_label=None, legacy=False, **kwargs)`

绑定到单个元素（及可选晶体相）的子类，自动实现上述 hooks。

- `structure`：晶体相键名（如 `'hcp'`/`'bcc'`/`'fcc'`），用于 ASE 自动生成
  POSCAR 兜底；缺省时取元素 `rt_structure`。传 `None` 显式关闭自动生成
  （仅用于外部 `poscar_path` 场景）。
- `phase_label`：路径/POSCAR 标签/job_name 用的相标签，默认取 `structure`。
  外部 POSCAR 场景可指定自定义标签如 `'exp_defect'`。
- `legacy`：`True` 时不附加 `_{phase_label}` 段，复现旧版路径布局。

### `run_eos_all(specs, work_root, **kwargs)`

```python
# 元素符号默认按 rt_structure 跑（适用于单相元素）
results = run_eos_all(['Al', 'Cu', 'Au'],
    work_root='/share/zeng/metals/dpmd',
    dpgen_dir='/share/zeng/metals/sample',  # outer; each element → {work}/{el}/dpgen
    machine_template='~/template/dpgen-machine.json')
# → {'Al': 'generated', 'Cu': 'generated', 'Au': 'generated'}

# 多相元素用 'Element-structure' 语法分别跑
results = run_eos_all(['Al', 'Ti-hcp', 'Ti-bcc'],
    work_root='/share/zeng/metals/dpmd', ...)
# → {'Al': 'generated', 'Ti-hcp': 'generated', 'Ti-bcc': 'generated'}
```

> 复杂场景（外部 POSCAR、`legacy=True` 读旧数据、自定义 `phase_label`）请直接
> 构造 `ElementEOSCalculator`，不走 `run_eos_all`。

---

## 配置参数

### 模型与结构文件

| 参数 | 默认值 | 说明 |
|---|---|---|
| `dpgen_dir` | None | Element-internal DPGEN 目录，即 `{work_root}/{element}/dpgen`（与 `DPBuilder.dpgen_dir` 一致） |
| `dp_model_path` | None | 显式 DP frozen_model 路径（跳过搜索） |
| `poscar_path` | None | 显式 POSCAR 路径（跳过搜索） |
| `poscar_dir` | None | POSCAR 目录，按 `{element}-HCP.POSCAR` / `{element}-LIQ.POSCAR` 等 label 精确匹配 |
| `structure` | (rt_structure) | 晶体相键名（`'hcp'`/`'bcc'`/...）；缺省按元素 rt_structure |
| `phase_label` | (= structure) | 路径/job_name 标签；外部 POSCAR 场景可自定义 |
| `legacy` | False | True → 路径不附加 `_{phase_label}` 段（旧版兼容） |

搜索优先级：
- **model**: `dp_model_path` > `dpgen_dir/{frozen_model.pb, frozen_model_compressed.pb}` symlink > `dpgen_dir/iter.*/00.train/000/`
- **POSCAR** (`role='solid_rt'`): `poscar_path` > `poscar_dir/{element}-{rt_structure}.POSCAR` > `dpgen_dir/../confs/` > ASE 自动生成
- **POSCAR** (`role='liquid'`): `poscar_dir/{element}-LIQ.POSCAR` > `dpgen_dir/../confs/{element}-LIQ.POSCAR` > 静默回退到 `solid_rt`

### Slurm 提交

| 参数 | 默认值 | 说明 |
|---|---|---|
| `machine_template` | None | dpgen `machine.json` 路径，**最高优先来源** |
| `partition` | None | Slurm 分区（覆盖 machine_template） |
| `nodes` | None | 节点数（覆盖 machine_template） |
| `ntasks_per_node` | None | 每节点任务数（覆盖 machine_template） |
| `wall_time` | None | 最长运行时间（覆盖 machine_template） |
| `gres` | None | GPU 资源 (e.g. `'gpu:1'`)（覆盖 machine_template） |
| `lmp_command` | None | LAMMPS 运行命令（覆盖 machine_template） |

优先级：**machine_template** → 构造函数非 `None` 参数 → 内部默认值

示例：`machine_template` 的 `model_devi` 指定 `cpu_per_node=40`、`queue_name=nudt40c`、`source_list=[...activate]`，则生成的 sbatch 自动用这些值，除非构造函数显式覆盖。

### 温度与计算

| 参数 | 默认值 | 说明 |
|---|---|---|
| `two_phase_frac` | 0.10 | 双相法温度偏移 = Tm × `two_phase_frac` |
| `two_phase_shift` | 300 | 双相法整体偏移 (K) |
| `two_phase_nz` | 10 | 双相法 z 方向复制数 |
| `npt_n` | 5 | NPT 温度点数 |
| `npt_dT` | 100 | NPT 温度间隔 (K) |
| `npt_shift` | 0 | NPT 相对 Tm 的整体偏移 (K)；默认 0 即以 Tm 为中心 |
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
> - NPT 温度: `933 + [-200, -100, 0, 100, 200]` = [733, 833, 933, 1033, 1133] K（跨越 Tm）

---

## 目录结构

```
{work_root}/{element}/
├── melt/
│   ├── {Tm1}k_{phase_label}/        # 两相熔化（每个候选温度一个目录）
│   │   ├── run.in          # LAMMPS 输入 (two-phase.j2 渲染)
│   │   ├── confs.data      # LAMMPS 结构文件
│   │   ├── cp.pb           # DeePMD 势函数
│   │   ├── job.sbatch      # Slurm 脚本
│   │   └── traj/           # dump 输出
│   ├── {Tm2}k_{phase_label}/
│   └── {Tm3}k_{phase_label}/
├── npt/
│   ├── {T1}k_{phase_label}_solid/
│   │   ├── run.in          # npt-solid.j2
│   │   ├── thermo.dat      # 热力学输出
│   │   └── rdf.txt         # RDF
│   ├── {T1}k_{phase_label}_liquid/
│   │   ├── run.in          # npt-liquid.j2
│   │   └── ...
│   ├── {T2}k_{phase_label}_solid/
│   └── ...
└── traj/
    ├── {Tm}k_{phase_label}_solid/
    │   ├── run.in          # nvt-solid-traj.j2
    │   └── traj/           # NVT dump 轨迹
    └── {Tm}k_{phase_label}_liquid/
```

`{phase_label}` 是晶体相标签，默认取 `structure`（即 `'hcp'`/`'bcc'`/`'fcc'`...）。
对单相元素（如 Al, FCC）该标签也自动附加（路径形如 `Al/melt/933k_fcc/`），
以便多元素批处理时保持一致的目录结构。要复现旧式无标签布局
（`melt/{T}k/`、`npt/{T}k_solid/`），构造 calculator 时传 `legacy=True`。

对于多相元素（如 Ti 同时有 HCP 和 BCC 两相），需分别为每个相各构造一个
calculator 实例，保证两条熔化曲线互不覆盖：

```python
for struct in ('hcp', 'bcc'):
    calc = ElementEOSCalculator('Ti', structure=struct,
                                work_root='/share/zeng/metals/dpmd',
                                dpgen_dir='/share/zeng/metals/dpmd/Ti/dpgen')
    calc.run_two_phase(submit=False)
```

外部 POSCAR 场景：用 `structure=None` 关闭 ASE 自动生成，用 `phase_label=`
显式指定路径/标签：

```python
calc = ElementEOSCalculator('Ti', structure=None, phase_label='exp_defect',
                            poscar_path='/path/to/my.POSCAR',
                            work_root='/share/zeng/metals/dpmd')
```

---

## 后处理

### 双相法分析

```python
result = calc.analyze_two_phase()
# {'results': {1133: 'solid', 1233: 'coexist', 1333: 'liquid'},
#  'coexist_temps': [1233],
#  'Tm_interval': (1233, 1233),
#  'Tm_refined': 1233,
#  'details': {1133: {...}, 1233: {...}, 1333: {...}}}
```

读取 `chunk.profile`（按 z 分层的 Q4/Q6/密度）比较上下半盒判定 `solid` /
`liquid` / `coexist`；可选地用 `rdf_top.txt` 交叉验证并标注 `confidence`。
`Tm_refined`（coexist 温度的中位数）会写回 `calc.Tm_refined`，供后续
`run_property_scans()` 自动使用。

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

### 批量双相法分析（旧 dump 路径，已弃用）

```python
from extrempy.md.traj import batch_analyze_two_phase  # deprecated
results = batch_analyze_two_phase('/path/to/dump/parent/dir')
```

新代码请用 `extrempy.campaign.chunk.diagnose_case` 基于 `chunk.profile` + RDF。

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
