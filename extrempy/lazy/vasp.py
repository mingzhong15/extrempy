import numpy as np
import os
import shutil
import json

import dpdata

from extrempy.constant import kb, J2eV, kb_eV
from .base import InputGenerator, fermi_dirac

# KNOWN ISSUE: use_ele_temp=1 expects fparam in eV, but this project stores
# fparam in K (read from INCAR TEBEG). This is deferred to a future fix.
# If DP model quality is poor, check this first.


def _incar_dict(preset_name, *, encut, nbands, ele_temp_K, latt_temp_K=None,
                md_steps=500, dt=1, is_magnetic=False):
    if latt_temp_K is None:
        latt_temp_K = ele_temp_K
    base = {
        'ISTART': 0, 'ICHARG': 2,
        'LWAVE': '.FALSE.', 'LCHARG': '.FALSE.',
    }
    base['ENCUT'] = encut
    base['NELM'] = 100
    base['ALGO'] = 'Normal'
    base['PREC'] = 'High'
    base['ISMEAR'] = -1
    base['SIGMA'] = '%.8f' % (ele_temp_K * kb_eV)
    base['EDIFF'] = '1E-6'
    base['NBANDS'] = nbands
    base['GGA'] = 'PS'
    base['LREAL'] = 'Auto'
    base['KPAR'] = 4
    base['NPAR'] = 4
    base['KSPACING'] = 0.5
    base['KGAMMA'] = '.TRUE.'
    if preset_name == 'scf':
        base['IBRION'] = -1
        base['NSW'] = 0
    elif preset_name == 'relax':
        base['IBRION'] = 1
        base['NSW'] = md_steps
        base['ISIF'] = 2
    elif preset_name == 'aimd-ttm':
        base['IBRION'] = 0
        base['ISIF'] = 2
        base['NSW'] = md_steps
        base['POTIM'] = '%.2f' % dt
        base['TEBEG'] = int(latt_temp_K)
        base['TEEND'] = int(latt_temp_K)
        base['SMASS'] = 0
        base['MDALGO'] = 2
    else:
        raise ValueError(f"Unknown preset: {preset_name}")
    if is_magnetic:
        base['ISPIN'] = 2
    return base


def _render_incar(d):
    d = dict(d)  # work on a copy
    lines = []
    lines.append('# CONTROL')
    for k in ('ISTART', 'ICHARG', 'LWAVE', 'LCHARG'):
        if k in d:
            lines.append('%s = %s' % (k, d.pop(k)))
    lines.append('')
    lines.append('# ELECTRON')
    for k in ('ENCUT', 'NELM', 'ALGO', 'PREC', 'ISMEAR', 'SIGMA',
              'EDIFF', 'NBANDS'):
        if k in d:
            lines.append('%s = %s' % (k, d.pop(k)))
    lines.append('')
    lines.append('# XC FUNCTIONAL')
    for k in ('GGA',):
        if k in d:
            lines.append('%s = %s' % (k, d.pop(k)))
    lines.append('')
    lines.append('# ION')
    for k in ('IBRION', 'NSW', 'ISIF', 'POTIM', 'TEBEG', 'TEEND',
              'SMASS', 'MDALGO'):
        if k in d:
            lines.append('%s = %s' % (k, d.pop(k)))
    lines.append('')
    lines.append('# PARALLIZATION')
    for k in ('LREAL', 'KPAR', 'NPAR'):
        if k in d:
            lines.append('%s = %s' % (k, d.pop(k)))
    lines.append('')
    lines.append('# K-POINTS')
    for k in ('KSPACING', 'KGAMMA'):
        if k in d:
            lines.append('%s = %s' % (k, d.pop(k)))
    lines.append('')
    for k, v in d.items():
        lines.append('%s = %s' % (k, v))
    lines.append('')
    return '\n'.join(lines)


class VASPGenerator(InputGenerator):

    def __init__(self, *arg, work_path=None, poscar_file=None,
                 potcar_lib_path=None, potcar_map=None, potcar_set='PBE54',
                 **kwargs):
        if work_path is not None:
            kwargs['work_path'] = work_path
        super().__init__(*arg, **kwargs)
        self.poscar_file = poscar_file
        self.potcar_lib_path = potcar_lib_path
        self.zval_list = []
        self.numb_atom = []
        self.atom_names = []
        if potcar_map is not None:
            self.potcar_map = potcar_map
        elif potcar_lib_path is not None:
            from .potcar_map import PotcarMap
            self.potcar_map = PotcarMap(potcar_set, potcar_lib_path)
        else:
            self.potcar_map = None
        if poscar_file is not None:
            self._generate_poscar()
            self._generate_potcar()

    def _generate_poscar(self):
        try:
            sys = dpdata.System(self.poscar_file, fmt='vasp/poscar')
            self.numb_atom = sys.get_atom_numbs()
            self.atom_names = sys.get_atom_names()
            shutil.copy(self.poscar_file,
                        os.path.join(self.work_path, 'POSCAR'))
        except Exception:
            print('POSCAR file is not found or unreadable')

    def _generate_potcar(self):
        if self.potcar_map is None:
            raise RuntimeError(
                "potcar_map or potcar_lib_path must be provided "
                "to generate POTCAR.")
        self.zval_list, _ = self.potcar_map.write_potcar(
            self.atom_names,
            os.path.join(self.work_path, 'POTCAR'))

    def set_params(self, encut=600, ele_temp=300, nbands=None,
                   scale=1.2, nband_min=5):
        self.encut = encut
        self.ele_temp = ele_temp
        total_electrons = 0.0
        total_atoms = 0
        for zval, nat in zip(self.zval_list, self.numb_atom):
            total_electrons += zval * nat
            total_atoms += nat
        nbands_0 = int(total_electrons / 2 * scale) + nband_min * total_atoms
        if nbands is None:
            self.nbands = nbands_0
        else:
            self.nbands = nbands

    def generate_incar(self, md_steps=500, dt=1, latt_temp=None,
                       mode='scf', preset=None, override=None,
                       is_magnetic=False):
        if latt_temp is None:
            latt_temp = self.ele_temp
        if preset is None:
            if mode in ('scf', 'relax'):
                preset = mode
            elif mode == 'md':
                preset = 'aimd-ttm'
            else:
                preset = mode
        d = _incar_dict(preset, encut=self.encut, nbands=self.nbands,
                        ele_temp_K=self.ele_temp, latt_temp_K=latt_temp,
                        md_steps=md_steps, dt=dt, is_magnetic=is_magnetic)
        if override:
            d.update(override)
        content = _render_incar(d)
        with open(os.path.join(self.work_path, 'INCAR'), 'w') as f:
            f.write(content)

    @classmethod
    def from_existing(cls, work_path, poscar_path, potcar_path):
        gen = cls(work_path=work_path, poscar_file=None)
        gen.poscar_file = poscar_path
        sys = dpdata.System(poscar_path, fmt='vasp/poscar')
        gen.numb_atom = sys.get_atom_numbs()
        gen.atom_names = sys.get_atom_names()
        gen.zval_list = []
        with open(potcar_path, 'r') as f:
            for line in f:
                if 'ZVAL' in line:
                    val_str = line.split('=')[-1].strip().rstrip(';')
                    gen.zval_list.append(float(val_str))
        if len(gen.zval_list) != len(gen.atom_names):
            gen.zval_list = [0] * len(gen.atom_names)
            print("WARNING: ZVAL count mismatch in from_existing, "
                  "using 0. set_params will produce wrong NBANDS.")
        shutil.copy(poscar_path, os.path.join(work_path, 'POSCAR'))
        shutil.copy(potcar_path, os.path.join(work_path, 'POTCAR'))
        return gen


class VASPReader:

    def __init__(self, work_path, is_prinft=False):
        self.work_path = work_path
        self.is_prinft = is_prinft

    def _read_outcar(self):
        with open(os.path.join(self.work_path, 'OUTCAR'), 'r') as f:
            count = 0
            while True:
                line = f.readline()
                if not line:
                    break
                if 'band No.  band energies     occupation' in line \
                        and count == 0:
                    band = []
                    band_energy = []
                    band_occupation = []
                    while True:
                        line = f.readline()
                        if len(line.split()) == 3:
                            tmp = line.split()
                            band.append(int(tmp[0]))
                            band_energy.append(float(tmp[1]))
                            band_occupation.append(float(tmp[2]))
                        else:
                            break
                    count += 1
                    self.band = np.vstack(
                        (band, band_energy, band_occupation)).T
                    if self.is_prinft:
                        print('%.d Band data read' % (self.band.shape[0]))
                if 'SIGMA = ' in line:
                    self.sigma = float(line.split()[-1])
                    self.ele_temp = self.sigma / kb / J2eV
                    if self.is_prinft:
                        print('Fermi-Dirac smearing: %.12f eV = %.2f K'
                              % (self.sigma, self.ele_temp))
                if 'E-fermi' in line:
                    self.efermi = float(line.split()[2])
                    if self.is_prinft:
                        print('Fermi energy: %.12f eV' % self.efermi)
                if 'free  energy' in line:
                    self.free_energy = float(line.split()[-2])
                    if self.is_prinft:
                        print('Free energy: %.12f eV' % self.free_energy)
                if 'energy  without entropy' in line:
                    self.internal_energy = float(line.split()[-4])
                    if self.is_prinft:
                        print('Internal energy: %.12f eV'
                              % self.internal_energy)
                if 'external pressure' in line:
                    self.press = float(line.split()[3])
                    if self.is_prinft:
                        print('Pressure: %.12f kBar' % self.press)
                if 'Total CPU time used' in line:
                    self.cpu_time = float(line.split()[-1])
                    if self.is_prinft:
                        print('CPU time: %.12f s' % self.cpu_time)
                    break

        if hasattr(self, 'internal_energy') and hasattr(self, 'free_energy') and hasattr(self, 'ele_temp'):
            self.entropy_product = self.internal_energy - self.free_energy
            self.ele_entropy = self.entropy_product / self.ele_temp
            if self.is_prinft:
                print('Entropy product Te*Se: %.12f eV' % self.entropy_product)
                print('Electronic entropy Se: %.12e eV/K' % self.ele_entropy)

    @staticmethod
    def parse_outcar_frames(outcar_path):
        """Parse OUTCAR, return (free_energy, internal_energy, sigma_eV, ele_temp_K) per ionic step.

        Uses 'FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)' as the ionic step
        boundary marker, which matches dpdata's frame separation.
        """
        frames = []
        sigma = None
        try:
            fh = open(outcar_path)
        except (FileNotFoundError, IOError):
            return frames
        with fh:
            for line in fh:
                if 'SIGMA = ' in line and sigma is None:
                    sigma = float(line.split()[-1])
                if 'FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)' not in line:
                    continue
                try:
                    line = next(fh)
                    line = next(fh)  # free energy TOTEN line
                    if 'free  energy' not in line:
                        continue
                    A = float(line.split()[-2])
                    line = next(fh)  # blank line
                    line = next(fh)  # energy without entropy line
                    if 'without entropy' not in line:
                        continue
                    U = float(line.split()[-4])
                except StopIteration:
                    break
                except ValueError:
                    continue
                if sigma is not None:
                    ele_temp_K = sigma / kb_eV
                else:
                    ele_temp_K = 0.0
                frames.append((A, U, sigma or 0.0, ele_temp_K))
        return frames

    def _plot_band(self, ax, cc='dimgray', ele_temp=None):
        ax.plot(self.band[:, 1], self.band[:, 2] / 2, 'o', ms=4,
                mew=0.5, color=cc, mfc='none')
        ax.axvline(self.efermi, color=cc, ls=':', lw=1.0)
        xx = np.linspace(self.band[:, 1].min(), self.band[:, 1].max(), 1000)
        if ele_temp is None:
            ele_temp = self.ele_temp
        yy = fermi_dirac(xx, self.efermi, ele_temp)
        ax.plot(xx, yy, '-', color=cc, lw=1.0)
        y0 = np.zeros_like(xx)
        ax.fill_between(xx, y0, yy, facecolor=cc, alpha=0.1)
        ax.set_xlabel('$E$ (eV)')
        ax.set_ylabel('$g(E)$')
        ax.set_xlim(xx.min(), )
        cri = self.band[:, 2] > 0.0
        self.min_nband = (self.band[:, 1][cri]).shape[0]
        print(self.min_nband,
              'bands are needed to be occupied with > 1e-5 fraction')
