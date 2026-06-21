import os
import shutil

import numpy as np
import dpdata


def _raw_to_set(out_dir, nline_per_set=100000):
    raw_names = ['box.raw', 'coord.raw', 'energy.raw', 'force.raw',
                 'virial.raw', 'atom_ener.raw', 'fparam.raw', 'aparam.raw']
    for raw_name in raw_names:
        raw_path = os.path.join(out_dir, raw_name)
        if not os.path.exists(raw_path):
            continue
        data = np.loadtxt(raw_path)
        nframe = data.shape[0] if data.ndim > 0 else 1
        for i in range(0, nframe, nline_per_set):
            chunk = data[i:i + nline_per_set]
            set_dir = os.path.join(out_dir, f'set.{i // nline_per_set:03d}')
            os.makedirs(set_dir, exist_ok=True)
            base = raw_name.replace('.raw', '')
            np.save(os.path.join(set_dir, base), chunk.astype(np.float32))
        os.remove(raw_path)


def _is_valid_data_dir(path):
    if not os.path.isdir(path):
        return False
    if not os.path.isfile(os.path.join(path, 'type.raw')):
        return False
    for entry in os.listdir(path):
        if entry.startswith('set.') and os.path.isdir(os.path.join(path, entry)):
            return True
    return False


def bootstrap_init_data(aimd_dirs, sample_root,
                        drop_first=200, solid_stride=50, liq_stride=30):
    """Read VASP AIMD OUTCAR, extract frames, write deepmd/npy + fparam.raw.

    Parameters
    ----------
    aimd_dirs : list[(label, aimd_dir)] or list[(label, aimd_dir, is_liquid)]
                if 2-tuple, liquid detection by label.endswith('-LIQ')
    sample_root : str  (init_data root)
    drop_first : int  steps to discard at start
    solid_stride : int  frame stride for solid phases
    liq_stride : int  frame stride for liquid phases

    Returns
    -------
    init_data_sys : list[str]  relative paths of generated sets
    """
    init_data_sys = []
    for item in aimd_dirs:
        if len(item) == 3:
            label, aimd_dir, is_liquid = item
        else:
            label, aimd_dir = item
            is_liquid = label.upper().endswith('-LIQ')
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
        out_dir = os.path.join(sample_root, label)
        if _is_valid_data_dir(out_dir):
            print(f"  SKIP {label}: data already exists in {out_dir}")
            init_data_sys.append(label)
            continue
        stride = liq_stride if is_liquid else solid_stride
        indices = list(range(drop_first, nframe, stride))
        sub = ss[indices]
        n_selected = len(indices)
        os.makedirs(out_dir, exist_ok=True)
        sub.to_deepmd_raw(out_dir)
        fpar = np.ones(n_selected).reshape(-1, 1) * fparam_K
        np.savetxt(os.path.join(out_dir, 'fparam.raw'), fpar)
        _raw_to_set(out_dir)
        init_data_sys.append(label)
        phase = 'LIQ' if is_liquid else 'SOL'
        print(f"  {label}: T_ref={fparam_K:.0f}K, "
              f"{n_selected}/{nframe} frames, {phase} stride={stride}")
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
