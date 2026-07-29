# Prompt: 单化合物相边界调研

## 任务

调研化合物 **{{FORMULA}}** 在深地条件下的全部多形（polymorph）及其热力学稳定区间。

## 输出要求

请按以下 JSON 格式输出，每个相一条。**只输出 JSON，不要解释文字**：

```json
[
  {
    "label": "{{FORMULA}}-<phase-name>",
    "phase_name": "<通用相名，如 bridgmanite, post-perovskite, B1, B2>",
    "spacegroup": <空间群号，整数，如 225>,
    "spg_intl": "<国际符号，如 Fm-3m, Pnma>",
    "structure_source": "mc3d",
    "P_range_GPa": [<下界>, <上界>],
    "T_range_K": [<下界>, <上界>],
    "refs": ["<文献1>", "<文献2>"]
  }
]
```

## 调研范围

- 压力：0 到 400 GPa（覆盖地幔 + 外核 + 内核）
- 温度：300 到 10000 K
- 相：包括所有实验已确认和理论预测的高压/高温相
- 液相：包含一个 LIQ 相，label 为 `{{FORMULA}}-LIQ`，structure_source 为 `from_aimd_contcar`

## 文献优先级

1. 实验相图（DAC, multi-anvil）
2. 第一性原理相图计算（USPEX, phonopy, AIMD）
3. 综述文章（Stixrude & Lithgow-Bertelloni, Duffy, Oganov）

## 注意

- `P_range_GPa` 和 `T_range_K` 是该相**热力学稳定**的区间（非采样区间）
- 如果文献给出相边界曲线（如 Clapeyron slope），请给出该相稳定区的包围盒
- 如果某相在低温高压稳定但高温 destabilize，请如实反映
- `refs` 至少 1 篇；优先 2 篇交叉验证
- 不要编造数据；如不确定，写 `null` 并在 refs 中标 "TODO"

## 示例（MgO 的参考输出）

```json
[
  {
    "label": "MgO-B1",
    "phase_name": "periclase (B1, NaCl-type)",
    "spacegroup": 225,
    "spg_intl": "Fm-3m",
    "structure_source": "mc3d",
    "P_range_GPa": [0, 400],
    "T_range_K": [300, 10000],
    "refs": ["Duffy et al. 1995 JGR", "Oganov et al. 2003 JCP"]
  },
  {
    "label": "MgO-B2",
    "phase_name": "CsCl-type (B2)",
    "spacegroup": 221,
    "spg_intl": "Pm-3m",
    "structure_source": "mc3d",
    "P_range_GPa": [400, 1500],
    "T_range_K": [300, 10000],
    "refs": ["Oganov et al. 2003 JCP"]
  },
  {
    "label": "MgO-LIQ",
    "phase_name": "liquid",
    "spacegroup": null,
    "spg_intl": null,
    "structure_source": "from_aimd_contcar",
    "P_range_GPa": [0, 400],
    "T_range_K": [3098, 10000],
    "refs": []
  }
]
```

## 当前待调研化合物

**{{FORMULA}}**

请开始调研并输出 JSON。
