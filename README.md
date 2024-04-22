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




