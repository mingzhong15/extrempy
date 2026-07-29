# Prompt: 化合物清单调研

## 任务

给定一组元素，请列出这些元素组合形成的**全部已知深地矿物/化合物**及其多形。

## 输入

- 元素集合：{{ELEMENTS}}
- 压力范围：0-400 GPa
- 温度范围：300-10000 K

## 输出要求

请按以下 JSON 格式输出，每个化合物一条：

```json
[
  {
    "formula": "<化学式，如 MgSiO3>",
    "common_name": "<通用名，如 bridgmanite>",
    "Tm_at_1bar_K": <常压熔点，整数>,
    "phases": [
      {
        "label": "<FORMULA>-<phase-name>",
        "phase_name": "<相名>",
        "spacegroup": <整数或 null>,
        "spg_intl": "<国际符号或 null>",
        "P_range_GPa": [<下界>, <上界>],
        "T_range_K": [<下界>, <上界>],
        "refs": ["<文献>"]
      }
    ]
  }
]
```

## 调研范围

- 二元化合物（如 MgO, FeO, SiO2）
- 三元化合物（如 MgSiO3, Mg2SiO4, FeSiO3）
- 纯元素（如 Fe, Ni, Si）的深地多形
- 包括：硅酸盐、氧化物、碳酸盐、硫化物、氮化物、氢化物

## 深地矿物学优先级

1. **下地幔主要矿物**：bridgmanite, ferropericlase, Ca-perovskite
2. **过渡带矿物**：wadsleyite, ringwoodite, majorite
3. **上地幔矿物**：olivine, pyroxene, garnet
4. **D''层**：post-perovskite
5. **外核/内核**：Fe-Ni 合金相（hcp, bcc, fcc）
6. **挥发分载体**：H2O, CH4, CO2, NH3 等在矿物中的存在形式

## 注意

- 只列出在 0-400 GPa 范围内有稳定性的化合物
- 每个化合物至少列出 2 个相（固相 + 液相）
- `refs` 至少 1 篇
- 不要编造；不确定的标 `null`

## 示例（Mg-Si-O 体系）

```json
[
  {
    "formula": "MgO",
    "common_name": "periclase",
    "Tm_at_1bar_K": 3098,
    "phases": [
      {"label": "MgO-B1", "phase_name": "periclase",
       "spacegroup": 225, "spg_intl": "Fm-3m",
       "P_range_GPa": [0, 400], "T_range_K": [300, 10000],
       "refs": ["Duffy 1995"]}
    ]
  }
]
```

## 当前待调研元素集合

**{{ELEMENTS}}**

请列出这些元素组合的全部深地矿物，并输出 JSON。
