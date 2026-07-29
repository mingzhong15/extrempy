# 深地大模型宽区间 (P,T) 自动化采样 — 实施计划

> **目标**：构建覆盖地球内部主流元素和矿物的第一性原理数据集，用于训练深地大模型。
>
> **考核指标**：
> - KPI 1.1（热力学区间）：T = 300 K → 10000 K，P = 0.1 MPa → 400 GPa
> - KPI 1.2（化学组分）：至少包含 16 种元素 — Mg, Si, O, Fe, Ca, Ni, H, C, N, Al, He, Na, P, S, K, Ti
>
> **本计划状态**：已完成方案审查与决策收敛，待启动实施。

---

## 一、设计原则（审查后确立）

1. **极简启动**：3 个新文件 + 4 个文件修改即可跑通核心流程，砍掉一切"以后可能有用"的字段和模块。
2. **派生量不存储**：JSON 只存热力学原始数据（`P_range`/`T_range`/`refs`），`P_explore`/`T_explore`/`earth_layers` 等全部运行时派生。
3. **运行时状态不进 JSON**：模型路径、数据路径由构造函数参数传入，不写入静态数据库。
4. **薄继承**：`CompoundDPBuilder` 只 override 必要 hook，其余全继承 `DPBuilder`，不做过度封装。
5. **向后兼容**：现有 `ElementDPBuilder` 及其测试不破坏；2D (P,T) seg 字段为可选，缺失时回退到 1D 流程。
6. **数据边界明确**：纯元素常压相留 `ELEMENT_PHASE_DATA`；16 深地元素（含高压相）+ 化合物进 `mineral_phases.json`。
7. **合并工作流解耦**：各化合物独立跑完整 DPGEN 流程；统一大模型的合并是独立后续工作流，本轮不实现。

---

## 二、架构总览

```
extrempy/
├── data/                                  ★ 新建
│   ├── mineral_phases.json                ★ 深地热力学数据库 (极简, git-tracked)
│   └── literature_mining/                 ★ DS v4 Flash 调研工作区
│       ├── prompts/
│       │   ├── phase_boundary.md          ★ 单相 P-T 边界调研提示词
│       │   └── compound_survey.md         ★ 单化合物全部相清单提示词
│       └── targets/
│           └── mg_sio_system.yaml         ★ Mg-Si-O 调研目标列表
├── lazy/
│   ├── minerals.py                        ★ 新建: 加载/查询/派生
│   ├── mc3d.py                            改: +parse_formula, +get_phases_compound
│   ├── dpgen.py                           改: +_generate_press_grid_for_phase, 升级 _set_model_devi_jobs_from_segments
│   ├── potcar_map.py                      改: POTCAR_MAP 补 9 元素
│   └── lib.py                             改: 删 Mg2SiO4-fo1/fo2/fo3 (529-679)
├── campaign/
│   ├── compound_dp.py                     ★ 新建: CompoundDPBuilder (薄继承)
│   └── single_element_dp.py               不动 (向后兼容)
└── test/
    ├── test_minerals.py                   ★ 新建
    ├── test_mc3d_compound.py              ★ 新建
    └── test_compound_dp.py                ★ 新建
```

**砍掉的模块**（审查后删除）：
- ~~`extreme_dp.py`~~ → 用 `compound_dp.py` 替代
- ~~`earth_layers.json`~~ → 并入 `minerals.py` 作为 Python 常量
- ~~`cli/minerals.py`~~ → 延后，初期手动编辑 JSON + `validate_db()`
- ~~`run_ds_survey.py` / `aggregate_to_db.py`~~ → 延后，初期只写 prompt 模板

---

## 三、`mineral_phases.json` 数据库

### 3.1 极简 schema（6 字段/相）

```json
{
  "version": "0.1",
  "compounds": {
    "<formula>": {
      "Tm_at_1bar_K": <int>,
      "phases": [
        {
          "label": "<unique-label>",
          "mc3d_uuid": "<uuid>" | null,
          "structure_source": "mc3d" | "manual" | "from_aimd_contcar",
          "P_range_GPa": [<lo>, <hi>],
          "T_range_K": [<lo>, <hi>],
          "refs": ["<citation>", ...]
        }
      ]
    }
  }
}
```

### 3.2 不存储的派生字段（运行时计算）

| 派生字段 | 计算方式 |
|---|---|
| `P_explore_GPa` | `P_range_GPa` ± 15% margin（margin 可配置） |
| `T_explore_K` | `T_range_K` ± 300 K margin |
| `earth_layers` | 由 `P_range_GPa` 查 `EARTH_LAYERS` 常量 |
| `confidence` / `manual_refined_by` | 等真需要分级时再加 |

### 3.3 初始内容

**纯元素（16 深地元素，含高压相）**：
- Mg, Si, O, Fe, Ca, Ni, H, C, N, Al, He, Na, P, S, K, Ti
- 每个元素列出常压相 + 已知高压相（如 Fe-bcc / Fe-hcp / Fe-LIQ）

**Mg-Si-O 体系化合物（试点，~15 相）**：

| Compound | Phases | P (GPa) | T (K) | Layer | Ref |
|---|---|---|---|---|---|
| MgO | B1 (periclase) | 0-400 | 300-10k | 全部 | Duffy 1995 |
| MgO | B2 (CsCl) | >400 | 300-10k | 内核 | Oganov 2003 |
| SiO2 | stishovite | 9-50 | 300-3000 | 上/过渡带 | Stishov 1961 |
| SiO2 | CaCl2 | 50-100 | 300-3000 | 下地幔 | Tsuchida 1962 |
| SiO2 | α-PbO2 (seifertite) | 100-400 | 300-3000 | 下/外核 | Murakami 2003 |
| MgSiO3 | enstatite | 0-10 | 300-2000 | 上地幔 | Presnall |
| MgSiO3 | majorite | 10-25 | 1000-2500 | 过渡带 | Gasparik |
| MgSiO3 | akimotoite | 18-25 | 1500-2200 | 过渡带 | Kato |
| MgSiO3 | bridgmanite | 23-125 | 1500-4000 | 下地幔 | Murakami 2004 |
| MgSiO3 | post-perovskite | 125-400 | 2500-5000 | D''/外核 | Oganov 2004 |
| Mg2SiO4 | forsterite (α) | 0-13 | 300-2000 | 上地幔 | Boyd 1964 |
| Mg2SiO4 | wadsleyite (β) | 13-18 | 500-2000 | 过渡带 | Ringwood 1975 |
| Mg2SiO4 | ringwoodite (γ) | 18-23 | 800-2200 | 过渡带 | Akaogi |
| 每个 compound | LIQ | 0-400 | Tm-10k | — | — |

### 3.4 手动精修工作流

```bash
# 1. 查看
python -c "from extrempy.lazy.minerals import list_compounds; list_compounds()"

# 2. 编辑 mineral_phases.json (vscode/jupyter)
#    修改 P_range_GPa, T_range_K, refs

# 3. 校验
python -c "from extrempy.lazy.minerals import validate_db; validate_db()"

# 4. git diff 检查 + 提交
git diff extrempy/data/mineral_phases.json
git commit -m "data: refine MgSiO3 bridgmanite P-T range per Murakami 2004"
```

**精修原则**：
- JSONDiff-friendly：每次只改几个字段，git diff 一目了然
- refs 强制：`validate_db()` 对无 refs 的条目警告
- 气体元素（H/He/N/P/S/O）只在化合物中出现，不作为单独 compound

---

## 四、`lazy/minerals.py` 模块

### 4.1 常量

```python
EARTH_LAYERS = {
    'upper_mantle':    {'P_GPa': (0, 13),   'T_K': (300, 2000)},
    'transition_zone': {'P_GPa': (13, 23),  'T_K': (500, 2200)},
    'lower_mantle':    {'P_GPa': (23, 125), 'T_K': (1500, 4000)},
    'd_prime':         {'P_GPa': (125, 135),'T_K': (2500, 4000)},
    'outer_core':      {'P_GPa': (135, 330),'T_K': (4000, 6000)},
    'inner_core':      {'P_GPa': (330, 400),'T_K': (5000, 10000)},
}

DEEP_EARTH_ELEMENTS = ['Mg','Si','O','Fe','Ca','Ni','H','C',
                       'N','Al','He','Na','P','S','K','Ti']

GLOBAL_P_RANGE_GPa = (0, 400)
GLOBAL_T_RANGE_K = (300, 10000)
```

### 4.2 函数

```python
def parse_formula(formula: str) -> dict[str, int]
    """'MgSiO3' → {'Mg':1, 'Si':1, 'O':3}"""

def load_mineral_db(json_path=None) -> dict
    """importlib.resources 加载, 默认 extrempy/data/mineral_phases.json"""

def list_compounds(db=None) -> list[str]

def get_compound(formula: str, db=None) -> dict
    """返回 {Tm_at_1bar_K, phases: [...]}"""

def get_compound_phases(formula: str, db=None) -> list[dict]
    """返回 phases 列表"""

def make_extreme_segments(formula: str, db=None, *,
                          p_margin=0.15, t_margin=300) -> list[dict]
    """构造 2D (P,T) segs.
    每个 seg 含:
      label, mc3d_uuid, structure_source,
      T_core, T_explore, P_core, P_explore,
      refs
    P_explore = P_range ± p_margin; T_explore = T_range ± t_margin
    """

def lookup_earth_layers(P_range_GPa) -> list[str]
    """由 P_range 查 EARTH_LAYERS, 返回匹配的层名列表"""

def validate_db(db=None) -> list[str]
    """返回 warnings 列表:
      - 无 refs 的相
      - P_range/T_range 缺失
      - mc3d_uuid=null 且 structure_source!='manual' (提示)
    """
```

### 4.3 入口统一

seg 构造**只在 `minerals.py` 的 `make_extreme_segments()` 里**。`mc3d.py` 不再提供独立的 `make_compound_segments_from_db()`，避免函数命名混乱。

---

## 五、`lazy/mc3d.py` 扩展

### 5.1 新增函数

```python
def parse_formula(formula: str) -> tuple[dict[str, int], set[str]]
    """'MgSiO3' → ({'Mg':1,'Si':1,'O':3}, {'Mg','Si','O'})"""

def get_phases_compound(formula: str, *,
                        method='pbesol-v2',
                        mode='all') -> list[dict]
    """按 formula 精确匹配 MC3D, 含 HP/HT/theoretical 相.
    返回字段与 get_phases() 一致 + phase_type."""

def backfill_mc3d_uuid(db: dict, method='pbesol-v2') -> dict
    """对 db 中 mc3d_uuid=null 且 structure_source='mc3d' 的相,
    按 (formula, spacegroup, spg_intl) 在 MC3D 查找匹配.
    返回 {'matched': int, 'unmatched': [labels]}.
    修改 db in-place."""
```

### 5.2 保留不动

- `get_phases()` 单元素查询（`ElementDPBuilder` 仍用）
- `make_phase_segments()` 单元素 seg 构造（`ElementDPBuilder` 仍用）
- `download_atoms()` / `list_phases()`

---

## 六、`lazy/dpgen.py` 扩展

### 6.1 新增函数

```python
def _generate_press_grid_for_phase(seg: dict, *, n_points=None) -> list[int] | None
    """seg 有 P_explore_GPa → log-spaced grid (bar); 无 → None.
    n_points 自适应: max(5, int(log10(P_hi/P_lo) * 3)).
    返回值单位 bar (DPGEN convention)."""
```

### 6.2 升级 `_set_model_devi_jobs_from_segments`

```python
def _set_model_devi_jobs_from_segments(self, segs, *,
                                       nsteps_per_phase=5,
                                       init_steps=None,
                                       press_grid=None,       # 全局 fallback
                                       n_press_per_phase=None, # 新增, 自适应
                                       trj_freq=20,
                                       numb_frame_per_iter_per_PT=5,
                                       ensemble='npt',
                                       sub_indices=None):
    ...
    for sys_idx, seg in enumerate(segs):
        ...
        T_list = _generate_temp_list(T_explore[0], T_explore[1])
        # NEW: per-phase press grid
        p_grid = _generate_press_grid_for_phase(seg, n_points=n_press_per_phase)
        if p_grid is None:
            p_grid = press_grid or [1, 1e1, 1e2, 1e3, 1e4]  # 旧 fallback
        ...
```

### 6.3 向后兼容

- 旧 segs（`ElementDPBuilder` 产生的，无 `P_explore_GPa`）→ `_generate_press_grid_for_phase` 返回 None → 回退全局 `press_grid`
- 现有 `test_dpgen.py` 测试不破坏

---

## 七、`campaign/compound_dp.py` — CompoundDPBuilder

### 7.1 薄继承设计

```python
from extrempy.campaign.single_element_dp import DPBuilder

class CompoundDPBuilder(DPBuilder):
    """DP construction across wide (P,T) for deep-Earth compounds/elements.

    与 ElementDPBuilder 的区别:
      - formula 支持化合物 ('MgO', 'MgSiO3') 和纯元素 ('Fe')
      - 相段来源 mineral_phases.json (含高压相)
      - 2D (P,T) seg: 每个 seg 带 P_core/P_explore
      - per-phase press_grid (log-spaced within P_explore)
      - 不走 ASE bulk fallback (化合物结构复杂)

    KPI coverage:
      P: 0.1 MPa → 400 GPa
      T: 300 K → 10000 K
      Elements: 16 (Mg,Si,O,Fe,Ca,Ni,H,C,N,Al,He,Na,P,S,K,Ti)
    """

    def __init__(self, formula, *,
                 mineral_db_path=None,
                 mc3d_method='pbesol-v2',
                 target_atoms=100,
                 supercell=None,
                 n_press_per_phase=None,   # None → 自适应
                 **kwargs):
        super().__init__(**kwargs)
        self.formula = formula
        self.elements = sorted(parse_formula(formula).keys())
        self.work_dir = os.path.join(self.work_root, formula)
        self.mc3d_method = mc3d_method
        self.target_atoms = target_atoms
        self.supercell = supercell
        self.n_press_per_phase = n_press_per_phase
        self.mineral_db = load_mineral_db(mineral_db_path)
        self._validate_formula_in_db()

    # ── override 的 hook (仅 3 个) ──────────────────────
    def _get_tm(self):
        return get_compound(self.formula, self.mineral_db).get('Tm_at_1bar_K', 2000)

    def get_phase_segments(self):
        segs = make_extreme_segments(self.formula, self.mineral_db)
        self._print_segment_table(segs)
        return segs

    def generate_poscars(self, segs):
        """Per seg:
           - structure_source='mc3d' → mc3d_source(uuid, target_atoms)
           - structure_source='from_aimd_contcar' → LIQ placeholder
           - structure_source='manual' → resolve_poscar with manual path
        不走 ASE bulk fallback."""
        ...

    # ── 全继承的方法 ──────────────────────────────────
    # generate_init_aimd / submit_init_aimd / collect_init_data
    # generate_dpgen / submit_dpgen / collect_dpgen / inspect
```

### 7.2 不包含（审查后砍掉）

- ~~`reuse_existing_dp` 参数~~ → 用户直接传 `extra_init_root`
- ~~`is_magnetic` 参数~~ → 默认 False (VASP ISPIN=1)，需要时通过 `extra_params` 传
- ~~`existing_dp_path` 字段~~ → 运行时状态不进 JSON/类属性
- ~~合并逻辑~~ → 独立后续工作流

### 7.3 批量入口

```python
def build_deep_earth_dataset(work_root, formulas=None, **kwargs):
    """批量构建, 默认覆盖 Mg-Si-O 4 个端元."""
    if formulas is None:
        formulas = ['MgO', 'SiO2', 'MgSiO3', 'Mg2SiO4']
    ...
```

---

## 八、`lazy/potcar_map.py` 扩展

在 `POTCAR_MAP['PBE54']`（第 6-62 行）补 9 个元素：

```python
'O':  {'variant': '',    'ZVAL': 0},   # 高圧下可改 _pv
'Fe': {'variant': '_pv', 'ZVAL': 0},
'Ni': {'variant': '_pv', 'ZVAL': 0},
'H':  {'variant': '',    'ZVAL': 0},
'C':  {'variant': '',    'ZVAL': 0},
'N':  {'variant': '',    'ZVAL': 0},
'He': {'variant': '',    'ZVAL': 0},
'P':  {'variant': '',    'ZVAL': 0},   # 高圧可改 _pv
'S':  {'variant': '',    'ZVAL': 0},
```

**待验证**：本地 `~/potpaw_PBE.54/` 是否有这些 POTCAR 文件（特别是 `Fe_pv`, `Ni_pv`）。缺失则需补 VASP 官方 pp。此项可在阶段 E 之前验证。

---

## 九、`lazy/lib.py` 清理

删除 `Mg2SiO4-fo1/fo2/fo3` 硬编码（第 529-679 行），共 3 个 `elif` 分支。

**核查**：grep 确认仅 `lib.py` 内部引用，无其他文件依赖（已确认，仅 `lib.py:529/582/629` 三处）。

---

## 十、数据边界（明确）

| 数据源 | 覆盖范围 | 使用者 | 操作 |
|---|---|---|---|
| `ELEMENT_PHASE_DATA` (lib.py dict) | 96 元素常压相 | `ElementDPBuilder` (老流程) | **不动** |
| `mineral_phases.json` | 16 深地元素（含高压相）+ 化合物 | `CompoundDPBuilder` (新流程) | **新建** |

- 纯 Fe 的高压相（hcp-Fe）进 JSON
- `ELEMENT_PHASE_DATA['Fe']` 的常压相（bcc/fcc）不动
- 两套数据有少量重叠（常压相），但各自独立，不冲突
- `get_viable_elements()` 不改，现有 `test_lib.py` 测试不破坏

---

## 十一、合并工作流（独立后续，本轮不实现）

```
各化合物独立跑 CompoundDPBuilder
  → 各自 frozen_model.pb + collected/
  → 后续设计独立的合并工作流:
     collect all collected/ → 统一 DPGEN (DPA-2, 16 元素 type_map)
```

**决策**：独立 DPGEN 后合并。合并机制等各化合物数据就绪后再设计，不耦合到 `CompoundDPBuilder`。

---

## 十二、DS v4 Flash 文献调研

### 12.1 交互方式（已确认）

我准备 prompt 模板 + targets 列表，用户跑 DS v4 Flash 后返回 JSON，我聚合到 `mineral_phases.json` 并标记。

### 12.2 初期交付物（本轮）

- `data/literature_mining/prompts/phase_boundary.md` — 单相 P-T 边界调研提示词
- `data/literature_mining/prompts/compound_survey.md` — 单化合物全部相清单提示词
- `data/literature_mining/targets/mg_sio_system.yaml` — Mg-Si-O 调研目标

### 12.3 延后交付物

- `run_ds_survey.py` 批量调用脚本
- `aggregate_to_db.py` 聚合脚本
- 等调研流程稳定后再做

---

## 十三、阶段任务分解

### 阶段 A：数据基础设施（5-6 天，先行）

| # | 任务 | 文件 | 工作量 |
|---|---|---|---|
| A1 | 创建 `extrempy/data/` 目录 + `.gitignore` | 新建 | 0.5h |
| A2 | 写 `data/mineral_phases.json` Mg-Si-O 4 化合物骨架 + 16 元素骨架 | 新建 | 2-3 天 |
| A3 | 写 `lazy/minerals.py` (常量 + 加载/查询/派生/validate) | 新建 | 1 天 |
| A4 | 写 DS v4 Flash prompts + targets | `data/literature_mining/` | 0.5 天 |
| A5 | 删除 `lib.py:529-679` 的 Mg2SiO4 硬编码 | 改 | 0.5h |
| A6 | 写 `test/test_minerals.py` | 新建 | 0.5 天 |

### 阶段 B：MC3D 化合物查询扩展（2-2.5 天，与 A 并行）

| # | 任务 | 文件 | 工作量 |
|---|---|---|---|
| B1 | `mc3d.py` 加 `parse_formula`, `get_phases_compound`, `backfill_mc3d_uuid` | 改 | 1 天 |
| B2 | 写 `test/test_mc3d_compound.py` | 新建 | 0.5 天 |
| B3 | 跑 backfill 把 Mg-Si-O 各相 uuid 填回 JSON | 脚本 | 0.5 天 |

### 阶段 D：2D (P,T) seg + per-phase press_grid（2-2.5 天，与 A 并行）

| # | 任务 | 文件 | 工作量 |
|---|---|---|---|
| D1 | `dpgen.py` 加 `_generate_press_grid_for_phase` | 改 | 0.5 天 |
| D2 | 升级 `_set_model_devi_jobs_from_segments` 支持 per-phase press_grid | 改 | 1 天 |
| D3 | `test/test_dpgen.py` 加兼容性用例（旧 segs 仍走全局网格） | 改 | 0.5 天 |

### 阶段 C：POTCAR 扩展（0.5 天，E 之前必做）

| # | 任务 | 文件 | 工作量 |
|---|---|---|---|
| C1 | `potcar_map.py` 补 9 个元素 | 改 | 0.5 天 |
| C2 | 用户验证本地 `~/potpaw_PBE.54/` 有这些 POTCAR | (用户) | — |

### 阶段 E：CompoundDPBuilder（3-4 天，依赖 A+B+D+C）

| # | 任务 | 文件 | 工作量 |
|---|---|---|---|
| E1 | 写 `campaign/compound_dp.py` CompoundDPBuilder 主体 | 新建 | 1.5 天 |
| E2 | override `_get_tm`, `get_phase_segments`, `generate_poscars` | | 1 天 |
| E3 | `build_deep_earth_dataset` 批量入口 | | 0.5 天 |
| E4 | 写 `test/test_compound_dp.py` dry-run MgO | 新建 | 0.5 天 |

### 阶段 F：DS v4 Flash 调研（持续，用户跑）

| # | 任务 | 执行者 | 工作量 |
|---|---|---|---|
| F1 | 用户跑 DS v4 Flash 调研 Mg-Si-O 全部相边界 | 用户 | 持续 |
| F2 | 我聚合结果到 `mineral_phases.json` | 我 | 0.5 天/批 |

### 阶段 G1：FeO 试点（5-7 天，集群，依赖 E）

- 真正的新化合物端到端验证
- 用 `CompoundDPBuilder('FeO')` 跑通完整流程
- 与已知相边界对比

### 阶段 G2：扩展到 16 元素（持续）

| 群组 | 元素 | 关键化合物 |
|---|---|---|
| Mg-Si-O | Mg, Si, O | MgO, SiO2, MgSiO3, Mg2SiO4 |
| Fe 群 | Fe, Ni, S | Fe, FeO, Fe2SiO4, FeS2, Fe-Ni |
| Ca-Al | Ca, Al, Si | CaSiO3, Al2O3, CaAl2Si2O8 |
| Na-K | Na, K, Al, Si | NaAlSi3O8, KAlSi3O8 |
| 挥发性 | H, C, N, He | (在化合物中掺杂) |
| Ti-P-S | Ti, P, S | TiO2, SiP2, FeS |

---

## 十四、关键约束（已确认）

| 约束 | 决策 |
|---|---|
| 气体元素 (H/He/N/P/S/O) | 仅在化合物中出现，不跑纯气相 |
| Fe/Ni 磁性 | 非极化 (VASP ISPIN=1，默认) |
| 数据库存储 | Repo 内 JSON，git-tracked，手动精修留痕 |
| 已有 Mg-Si-O 数据 | 用户已有势函数；通过 `extra_init_root` 增量训练 |
| 试点化合物 | FeO (真正新化合物，启动 Fe 群) |
| 类设计 | 新建 `CompoundDPBuilder`，薄继承 `DPBuilder` |
| 数据边界 | 16 元素（含高压相）+ 化合物进 JSON；`ELEMENT_PHASE_DATA` 不动 |
| 合并工作流 | 独立 DPGEN 后合并，本轮不实现 |

---

## 十五、关键风险与对策

| 风险 | 对策 |
|---|---|
| POTCAR 库缺 9 元素 | 阶段 C2 验证；缺失则补 VASP 官方 pp |
| MC3D 缺 post-perovskite 等 HP 相 uuid | `structure_source='manual'` fallback，本地 POSCAR 一等公民 |
| 已有 Mg-Si-O 数据格式不兼容 | `extra_init_root` 增量训练，不重跑 |
| DS v4 Flash 调研结果不准 | `refs` 强制；至少 2 源交叉验证 |
| 0-400 GPa 跨 6 个量级 | per-phase log-spaced press_grid，自适应 n_points |
| Fe/Ni 非极化在低压低温不准 | 接受偏差（KPI 关注深地 T>1000K）；后续可补磁性 refinement |
| 16 元素统一大模型训练成本高 | DPA-2 + `training_reuse_iter` + 分阶段合并数据 |

---

## 十六、验收标准

### 阶段 A+B+D 完成时

- [ ] `extrempy/data/mineral_phases.json` 含 Mg-Si-O 4 化合物 + 16 元素骨架
- [ ] `extrempy/lazy/minerals.py` 通过 `test_minerals.py`
- [ ] `extrempy/lazy/mc3d.py` `get_phases_compound('MgO')` 返回正确相数
- [ ] `extrempy/lazy/dpgen.py` per-phase press_grid 通过 `test_dpgen.py`
- [ ] `extrempy/data/literature_mining/prompts/` 含 2 个 prompt 模板
- [ ] `lib.py` 的 Mg2SiO4 硬编码已删除
- [ ] 现有 `test_lib.py` / `test_dpgen.py` 全部通过（不破坏）

### 阶段 E 完成时

- [ ] `extrempy/campaign/compound_dp.py` CompoundDPBuilder 可实例化
- [ ] `CompoundDPBuilder('MgO').get_phase_segments()` 返回 2D segs
- [ ] `CompoundDPBuilder('MgO').generate_poscars(segs)` 成功（MC3D fetch）
- [ ] `test_compound_dp.py` dry-run 通过

### 阶段 G1 完成时

- [ ] FeO 完整流程跑通：POSCAR → AIMD → DPGEN → frozen_model
- [ ] 训练后的 DP 能复现已知 FeO 相边界（B1/B8 转变）

---

## 十七、启动顺序

**并行启动**：阶段 A + B + D（这三阶段不需 POTCAR 验证）

**串行依赖**：
```
A (数据) ──┐
B (MC3D) ──┼──→ E (CompoundDPBuilder) ──→ G1 (FeO 试点) ──→ G2 (扩展)
D (dpgen) ─┘                                    ↑
C (POTCAR) ──────────────────────────────────────┘
F (DS 调研) ── 持续并行 ──→ 聚合到 A 的 JSON
```

**总周期估算**：
- 2-3 周完成阶段 A-E（代码 + Mg-Si-O 数据库骨架）
- 持续阶段 F（DS v4 Flash 调研并行）
- 1-2 个月完成 G1 (FeO 试点) + 启动 G2 扩展
