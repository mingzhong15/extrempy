# extrempy

post-processing code for atomistic modeling (at extreme conditions)


## Installation
---

### from pip 

```
pip install extrempy
```

### from repository

to install the `extrempy` package directly from the repository, clone it from GitHub and use `pip` to install it:

```
git clone https://github.com/mingzhong15/extrempy.git
cd extrempy
pip install .
```

## Simple Guide
---
### visuallization of dataset

```python
from extrempy.dataset import SetSys

fig, ax1 = plt.subplots(figsize=(3,2),dpi=200)

ss = SetSys( SET_DIR, is_printf=False )
ss._read_thermo()

ax.plot( ss.pres, ss.temp, 'o', ms=3, mew=0.2, color='#1f77b4',alpha=0.6, mfc='none')

```

### visuallization of data sampling

we can use `sys._plot_model_devi` to visuallize model deviation for different iterations

```python
from extrempy.dpsample import SampleSys

for case_idx in case_list:

  fig, ax = plt.subplots(figsize=(3,1),dpi=200)

  for iter_idx in [0]:
      
      print("Iter.%.3d Case.%.3d"%(iter_idx, case_idx))

      for sys_idx in [0]:
          sys._plot_model_devi(ax, iter_idx = iter_idx, sys_idx = sys_idx, case_idx = case_idx)

```

we can use `sys._plot_all_sampling` to visuallize data sampling in (p,T) space

```python

sys = SampleSys(DIR, printf=False)  

fig, ax = plt.subplots(figsize=(5,3),dpi=200)

color_list = ['coral','crimson','firebrick']

sys._plot_all_sampling(ax, color = color_list)

```


we can collect the sampled data from each iterations (containing `fparam.npy`, `aparam.npy`)

```python

sys = SampleSys(DIR)
sys._collect_data(OUT_DIR, exe_path='/personal/raw_to_set.sh')

```


## Structure & Campaign Workflow
---

### unified structure resolver

`resolve_poscar` is the single entry point for obtaining a POSCAR file. it handles "already exists / LIQ placeholder / write file"; the actual structure production is delegated to a `source` (callable / `ase.Atoms` / path str).

```python
from extrempy import resolve_poscar, ase_source, mc3d_source

# ASE standard structure (fcc/bcc/hcp/...) — auto-detects rt_structure
resolve_poscar('Al', 'Al-FCC',
               confs_dir='./confs',
               source=ase_source('Al', structure_type='fcc'))

# MC3D fetch — requires structure_uuid
resolve_poscar('Al', 'Al-225',
               confs_dir='./confs',
               source=mc3d_source('uuid-from-mc3d', target_atoms=100))

# existing POSCAR file path
resolve_poscar('Al', 'Al-FCC',
               confs_dir='./confs',
               source='./external/Al-fcc.POSCAR')
```

labels ending in `-LIQ` are treated as liquid placeholders and skipped (returns `None`); they are populated later from AIMD CONTCAR by `DPBuilder.collect_init_data()`.

### DP potential construction — `ElementDPBuilder`

```python
from extrempy import ElementDPBuilder

b = ElementDPBuilder('Al',
    work_root    = '/share/zeng/metals',
    potcar_lib   = '~/potpaw_PBE.54',
    machine_template = '~/template/dpgen-machine.json',
    mc3d_mode    = 'ambient',     # fetch structures from MC3D
)

segs = b.get_phase_segments()     # phase segments from RT to 2*Tm
b.generate_poscars(segs)          # → confs/{label}.POSCAR (via resolve_poscar)
b.generate_init_aimd(segs)        # → init_vasp/{STRUCTURE}-{T}K/
b.submit_init_aimd()              # submit AIMD jobs

# after AIMD completes:
b.collect_init_data(segs)         # → init_data/ + confs/{el}-LIQ.POSCAR
b.generate_dpgen(segs)            # → dpgen/param.json + machine.json
b.submit_dpgen()                  # launch DPGEN

# after DPGEN converges:
b.collect_dpgen()                 # → dpgen/frozen_model.pb symlink + collected/
```

### EOS / melt calculation — `ElementEOSCalculator`

```python
from extrempy import ElementEOSCalculator

calc = ElementEOSCalculator('Al',
    work_root    = '/share/zeng/metals/dpmd',
    # element-internal DPGEN dir (consistent with DPBuilder.dpgen_dir):
    dpgen_dir    = '/share/zeng/metals/Al/dpgen',
    machine_template = '~/template/dpgen-machine.json',
)

calc.generate_two_phase()         # 3 candidate T around Tm
calc.submit_two_phase()
calc.analyze_two_phase()          # → Tm interval via Q4/Q6

calc.generate_npt()               # solid + liquid NPT scan
calc.submit_npt()
calc.analyze_npt()                # → DataFrame (T, volume, density, energy)

calc.generate_nvt_traj()          # NVT trajectory at Tm
calc.submit_nvt_traj()
```

POSCAR selection is role-based: `solid_rt` looks up `{element}-{rt_structure}.POSCAR` (e.g. `Al-FCC.POSCAR`); `liquid` looks up `{element}-LIQ.POSCAR` and silently falls back to solid_rt if absent. the DP model is found via the `collect_dpgen` symlink at `dpgen_dir/frozen_model.pb`, then `iter.*/00.train/000/`.

