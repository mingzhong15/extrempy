import os
import shutil
import subprocess
import glob

import numpy as np
import dpdata


def bootstrap_init_data(aimd_dirs, sample_root, raw_to_set_script,
                        drop_first=200, low_T_stride=50, high_T_stride=30,
                        high_T_threshold=1500):
    """Read VASP AIMD OUTCAR, extract frames, write deepmd/npy + fparam.raw.

    Parameters
    ----------
    aimd_dirs : list[(label, aimd_dir)]
    sample_root : str  (init_data root)
    raw_to_set_script : str  (path to raw_to_set.sh)
    drop_first : int  steps to discard at start
    low_T_stride : int  frame stride for T < high_T_threshold
    high_T_stride : int  frame stride for T >= high_T_threshold
    high_T_threshold : float  K

    Returns
    -------
    init_data_sys : list[str]  relative paths of generated sets
    """
    init_data_sys = []
    for label, aimd_dir in aimd_dirs:
        outcar = os.path.join(aimd_dir, 'OUTCAR')
        incar = os.path.join(aimd_dir, 'INCAR')
        if not os.path.exists(outcar):
            print(f"  SKIP {label}: OUTCAR not found in {aimd_dir}")
            continue
        ss = dpdata.LabeledSystem(outcar, fmt='vasp/outcar')
        nframe = ss.get_nframes()
        if nframe <= drop_first:
            print(f"  SKIP {label}: nframe={nframe} <= drop_first={drop_first}")
            continue
        fparam_K = 300.0
        if os.path.exists(incar):
            with open(incar, 'r') as f:
                for line in f:
                    if 'TEBEG' in line:
                        parts = line.split('=')
                        if len(parts) > 1:
                            try:
                                fparam_K = float(parts[1].strip().split()[0])
                            except ValueError:
                                pass
                        break
        stride = high_T_stride if fparam_K >= high_T_threshold else low_T_stride
        indices = list(range(drop_first, nframe, stride))
        sub = ss[indices]
        n_selected = len(indices)
        out_dir = os.path.join(sample_root, label)
        os.makedirs(out_dir, exist_ok=True)
        sub.to_deepmd_raw(out_dir)
        fpar = np.ones(n_selected).reshape(-1, 1) * fparam_K
        np.savetxt(os.path.join(out_dir, 'fparam.raw'), fpar)
        raw_list = glob.glob(os.path.join(out_dir, '*', '*.raw'))
        if raw_list:
            set_numb = 100000
            cmd = '%s %.d' % (raw_to_set_script, set_numb)
            ret = subprocess.run(cmd, shell=True, cwd=out_dir,
                                 capture_output=True)
            if ret.returncode != 0:
                stderr = ret.stderr.decode() if ret.stderr else ''
                print(f"  WARNING: {raw_to_set_script} returned "
                      f"{ret.returncode} for {label}")
                if stderr:
                    print(f"    stderr: {stderr[:200]}")
        init_data_sys.append(label)
        print(f"  {label}: T_ref={fparam_K:.0f}K, "
              f"{n_selected}/{nframe} frames, stride={stride}")
    return init_data_sys


def generate_liquid_poscar_from_contcar(aimd_dir_label, poscar_lib):
    """Copy CONTCAR (last frame) from AIMD dir to poscar_lib as LIQ POSCAR.

    Parameters
    ----------
    aimd_dir_label : tuple (label, aimd_dir)
        label should end with '-LIQ' for liquid phase
    poscar_lib : str  POSCAR library dir

    Returns
    -------
    liquid_poscar_path : str or None
    """
    label, aimd_dir = aimd_dir_label
    if not label.endswith('-LIQ'):
        return None
    contcar = os.path.join(aimd_dir, 'CONTCAR')
    if not os.path.exists(contcar):
        print(f"  WARNING: CONTCAR not found in {aimd_dir}, "
              f"cannot generate {label} POSCAR")
        return None
    dst = os.path.join(poscar_lib, label + '.POSCAR')
    shutil.copy(contcar, dst)
    print(f"  LIQ POSCAR generated: {dst}")
    return dst


def scale_poscar_volume(poscar_path, scale, out_path=None):
    """Scale POSCAR lattice vectors by factor."""
    ss = dpdata.System(poscar_path, fmt='vasp/poscar')
    coords = ss['coords'][0]
    orig_cell = ss['cells'][0]
    new_cell = orig_cell * scale
    atom_names = ss['atom_names']
    atom_numbs = ss['atom_numbs']
    symbols = []
    for name, count in zip(atom_names, atom_numbs):
        symbols.extend([name] * count)
    from ase import Atoms
    atoms = Atoms(symbols=symbols, positions=coords, cell=new_cell, pbc=True)
    if out_path is None:
        out_path = poscar_path
    atoms.write(out_path, format='vasp', direct=True)
    return out_path
