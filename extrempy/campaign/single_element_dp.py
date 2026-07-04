import os
import json
import shutil
import glob

import numpy as np

from extrempy.lazy.vasp import VASPGenerator, VASPReader, _incar_dict, _render_incar
from extrempy.lazy.dpgen import (DPGENGenerator,
                                 _generate_dpgen_machine_from_file)
from extrempy.lazy.lib import (get_phase_segments, get_viable_elements,
                               ELEMENT_PHASE_DATA, SUPPORTED_STRUCTURES)
from extrempy.lazy.potcar_map import PotcarMap
from extrempy.lazy.init_data import (bootstrap_init_data,
                                     generate_liquid_poscar_from_contcar,
                                     scale_poscar_volume)


class DPBuilder:
    """Base class for automated DP potential construction across wide T-P.

    Subclasses implement get_phase_segments(), generate_poscars(), build_potcar().
    """

    def __init__(self, work_root, potcar_lib, potcar_set='PBE54',
                 machine_template=None, job_template=None,
                 platform='bh',
                 encut=600, nband_scale=1.2, nband_min=5,
                 press_grid=None,
                 nsteps_per_phase=5,
                 init_steps=None,
                 f_trust=None, model_devi_skip=0,
                  training_reuse_iter=99, numb_frame_per_iter_per_PT=5,
                  trj_freq=20, aimd_steps=500, aimd_dt=1,
                  liquid_T_factor=1.8, liquid_V_scale=1.10):
        self.work_root = work_root
        self.potcar_lib = potcar_lib
        self.potcar_set = potcar_set
        self.machine_template = machine_template
        self.job_template = job_template
        self.encut = encut
        self.nband_scale = nband_scale
        self.nband_min = nband_min
        self.press_grid = press_grid or [1, 10, 100, 1000, 10000]  # bar; adjust upward for extreme high-P
        self.nsteps_per_phase = nsteps_per_phase
        self.init_steps = init_steps or [1000, 2000, 4000, 8000, 16000]
        self.f_trust = f_trust or [0.005, 0.3]
        self.model_devi_skip = model_devi_skip
        self.training_reuse_iter = training_reuse_iter
        self.numb_frame_per_iter_per_PT = numb_frame_per_iter_per_PT
        self.trj_freq = trj_freq
        self.aimd_steps = aimd_steps
        self.aimd_dt = aimd_dt
        self.liquid_T_factor = liquid_T_factor
        self.liquid_V_scale = liquid_V_scale
        self.platform = platform
        self.potcar_map = PotcarMap(potcar_set, potcar_lib)
        self.jparam = None
        self.work_dir = work_root
        self.element = None
        self.elements = []

    @property
    def confs_dir(self):
        return os.path.join(self.work_dir, 'confs')

    @property
    def init_vasp_dir(self):
        return os.path.join(self.work_dir, 'init_vasp')

    @property
    def init_data_dir(self):
        return os.path.join(self.work_dir, 'init_data')

    @property
    def dpgen_dir(self):
        return os.path.join(self.work_dir, 'dpgen')

    def _ensure_dirs(self):
        for d in [self.work_dir, self.confs_dir, self.init_vasp_dir,
                  self.init_data_dir, self.dpgen_dir]:
            os.makedirs(d, exist_ok=True)

    # ---- abstract hooks ----
    def get_phase_segments(self):
        raise NotImplementedError

    def generate_poscars(self, segs):
        raise NotImplementedError

    def build_potcar(self, elements):
        text, zvals, _ = self.potcar_map.build(elements)
        potcar_path = os.path.join(self.dpgen_dir, 'POTCAR')
        with open(potcar_path, 'wb') as f:
            f.write(text)
        return zvals

    # ---- init AIMD ----
    def generate_init_aimd(self, segs, elements=None, aimd_temps=None):
        """
        Generate AIMD input files for each phase segment.

        Parameters
        ----------
        segs : list[dict]
            Phase segments from get_phase_segments().
        elements : list[str] or None
        aimd_temps : list[int] or None
            Explicit AIMD temperature (K) for each seg.  Length must
            equal ``len(segs)``.  When None, the temperature is derived
            from each seg: solid phases use the T_core midpoint, liquid
            phases use ``liquid_T_factor * Tm`` (overheated to ensure
            melting).  When specifying manually for a liquid seg, ensure
            the temperature is high enough to melt the structure
            (typically >= 1.5*Tm).
        """
        if elements is None:
            elements = self.elements
        if aimd_temps is not None:
            if len(aimd_temps) != len(segs):
                raise ValueError(
                    f"aimd_temps length {len(aimd_temps)} != "
                    f"len(segs) {len(segs)}")
        self.potcar_map.build(elements, quiet=True)  # validate ZVAL, no file write
        self._ensure_dirs()
        Tm = self._get_tm()
        # POTCAR summary
        path, variant_used = self.potcar_map._resolve_path(elements[0])
        zval = PotcarMap._parse_zval(path)
        potcar_label = f"{elements[0]}{variant_used}" if variant_used else elements[0]
        print("-- POTCAR --")
        preferred = self.potcar_map.map[elements[0]].get('variant', '')
        if variant_used != preferred and preferred:
            print(f"  {potcar_label} (ZVAL={zval:g})  "
                  f"[!] prefers '{preferred}' variant, using '{variant_used}' (auto fallback)")
        else:
            print(f"  {potcar_label} (ZVAL={zval:g})")
        # AIMD init header
        print(f"-- AIMD Init ({len(segs)} jobs) --")
        self._aimd_dirs = []
        solid_label = None
        for i, seg in enumerate(segs):
            label = seg['label']
            is_liquid = label.endswith('-LIQ')
            # --- structure-related: poscar_ref and solid_label ---
            if is_liquid:
                if solid_label is None:
                    raise RuntimeError(
                        f"Liquid phase '{label}' has no preceding solid phase "
                        f"to seed the melt POSCAR from.")
                poscar_ref = solid_label
            else:
                poscar_ref = label
                solid_label = label
            # --- temperature: user-specified > derived ---
            if aimd_temps is not None:
                T_ref = int(aimd_temps[i])
            elif is_liquid:
                T_ref = int(self.liquid_T_factor * Tm)
            else:
                T_ref = int((seg['T_core'][0] + seg['T_core'][1]) / 2)
            job_label = f"{seg['structure'].upper()}-{T_ref}K"
            work_dir = os.path.join(self.init_vasp_dir, job_label)
            os.makedirs(work_dir, exist_ok=True)
            gen = VASPGenerator(work_path=work_dir,
                                poscar_file=os.path.join(self.confs_dir,
                                                         poscar_ref + '.POSCAR'),
                                potcar_map=self.potcar_map)
            if is_liquid:
                shutil.copy(os.path.join(self.confs_dir, poscar_ref + '.POSCAR'),
                            os.path.join(work_dir, 'POSCAR'))
                scale_poscar_volume(os.path.join(work_dir, 'POSCAR'),
                                    self.liquid_V_scale)
            gen.set_params(encut=self.encut, ele_temp=T_ref,
                           scale=self.nband_scale, nband_min=self.nband_min)
            gen.generate_incar(md_steps=self.aimd_steps, dt=self.aimd_dt,
                               latt_temp=T_ref, mode='aimd-ttm')
            self._aimd_dirs.append((job_label, work_dir, label))
            # Print temperature source for traceability.
            if aimd_temps is not None:
                src = "user-specified"
            elif is_liquid:
                src = f"1.8Tm={1.8*Tm:.0f}"
            else:
                src = "T_core midpoint"
            print(f"  [{i}] {job_label}: T_ref={T_ref}K ({src})")

    def submit_init_aimd(self, submit=True):
        status = "submitted" if submit else "not submitted"
        print(f"-- sbatch ({status}) --")
        for i, (jl, wd, *_) in enumerate(self._aimd_dirs):
            gen = VASPGenerator(work_path=wd, poscar_file=None)
            job_name = (self.element or 'system') + '-' + jl
            if self.platform == 'slurm':
                if self.machine_template and os.path.exists(self.machine_template):
                    gen.generate_submit(self.machine_template, job_name,
                                        platform='slurm')
                    if submit:
                        gen.submit()
                        print(f"  [{i}] {jl}: sbatch submitted")
                    else:
                        print(f"  [{i}] {jl} \u2192 job.sbatch")
                else:
                    print(f"  [{i}] {jl}: no machine_template, skip")
            elif self.platform == 'bh':
                if self.job_template and os.path.exists(self.job_template):
                    gen.generate_submit(self.job_template, job_name,
                                        platform='bh')
                    if submit:
                        gen.submit()
                        print(f"  [{i}] {jl}: bh submitted")
                    else:
                        print(f"  [{i}] {jl} \u2192 job.json")
                else:
                    print(f"  [{i}] {jl}: no job_template, skip")
            else:
                print(f"  [{i}] {jl}: unknown platform '{self.platform}', skip")

    # ---- collect init data ----
    def collect_init_data(self, segs, drop_first=200,
                          solid_stride=50, liq_stride=30):
        _dirs = [(l, w, orig.endswith('-LIQ'))
                 for l, w, orig in self._aimd_dirs]
        results = bootstrap_init_data(
            _dirs, self.init_data_dir,
            drop_first=drop_first,
            solid_stride=solid_stride,
            liq_stride=liq_stride)
        for _, wd, orig in self._aimd_dirs:
            generate_liquid_poscar_from_contcar((orig, wd), self.confs_dir)
        self._init_data_sys = results
        return results

    @staticmethod
    def _is_valid_data_dir(path):
        if not os.path.isdir(path):
            return False
        if not os.path.isfile(os.path.join(path, 'type.raw')):
            return False
        for entry in os.listdir(path):
            if entry.startswith('set.') and os.path.isdir(os.path.join(path, entry)):
                return True
        return False

    # ---- DPGEN ----
    def generate_dpgen(self, segs, elements=None,
                       extra_init_sys=None, extra_init_root=None,
                       phase_ids=None, phase_labels=None,
                       extra_params=None):
        if elements is None:
            elements = self.elements

        # Compute sub_indices for exploration (sys numbering stays consistent)
        sub_indices = None
        if phase_ids is not None:
            sub_indices = set(phase_ids)
        elif phase_labels is not None:
            sub_indices = {i for i, s in enumerate(segs)
                           if s['label'] in phase_labels}
        if sub_indices is not None and not sub_indices:
            raise ValueError("No phase segments match the given filter")

        self._ensure_dirs()
        zvals = self.build_potcar(elements)
        type_map = list(elements)
        g = DPGENGenerator(work_path=self.dpgen_dir, type_map=type_map)

        # Init data: AIMD collected (if any) + extra
        for entry in sorted(os.listdir(self.init_data_dir)):
            sub = os.path.join(self.init_data_dir, entry)
            if self._is_valid_data_dir(sub):
                rel = os.path.abspath(sub)
                if rel not in g.jparam['init_data_sys']:
                    g.jparam['init_data_sys'].append(rel)
        if extra_init_sys:
            for path in extra_init_sys:
                rel = os.path.abspath(path)
                if rel not in g.jparam['init_data_sys']:
                    g.jparam['init_data_sys'].append(rel)
        if extra_init_root:
            root = os.path.abspath(extra_init_root)
            for entry in sorted(os.listdir(root)):
                sub = os.path.join(root, entry)
                if self._is_valid_data_dir(sub):
                    rel = os.path.abspath(sub)
                    if rel not in g.jparam['init_data_sys']:
                        g.jparam['init_data_sys'].append(rel)

        # Sys configs: use segs ordering (not glob dict order) for sys_idx consistency
        sys_configs = []
        missing = []
        for seg in segs:
            p = os.path.abspath(os.path.join(self.confs_dir, seg['label'] + '.POSCAR'))
            if os.path.exists(p):
                sys_configs.append([p])
            else:
                missing.append(seg['label'])
        if missing:
            raise FileNotFoundError(
                f"Missing POSCAR for phases {missing}; run generate_poscars "
                f"and (for LIQ) collect_init_data first.")
        g.jparam['sys_configs'] = sys_configs
        g.jparam['sys_configs_prefix'] = ''

        g._set_model_training_settings(stop_batch=200000, is_ele_temp=True)
        g._set_model_devi_settings(dt=0.001, f_trust=self.f_trust,
                                   is_relative=False, epsilon=1.0)
        g.jparam["model_devi_skip"] = self.model_devi_skip
        g.jparam["training_reuse_iter"] = self.training_reuse_iter
        g.jparam["fp_accurate_threshold"] = 0.99
        g._set_model_devi_jobs_from_segments(
            segs, nsteps_per_phase=self.nsteps_per_phase,
            init_steps=self.init_steps, press_grid=self.press_grid,
            trj_freq=self.trj_freq,
            numb_frame_per_iter_per_PT=self.numb_frame_per_iter_per_PT,
            ensemble='npt', sub_indices=sub_indices)
        self._write_fp_incar(g)
        potcar_path = os.path.join(self.dpgen_dir, 'POTCAR')
        if not os.path.exists(potcar_path):
            self.potcar_map.write_potcar(elements, potcar_path)
        g.jparam["fp_pp_path"] = self.dpgen_dir
        g.jparam["fp_pp_files"] = ['POTCAR']
        g.jparam["fp_incar"] = os.path.join(self.dpgen_dir, 'INCAR')
        g._optimize_prefix()
        if extra_params:
            g.jparam.update(extra_params)
        param_path = os.path.join(self.dpgen_dir, 'param.json')
        with open(param_path, 'w') as f:
            json.dump(g.jparam, f, indent=4)
        print(f"  param.json written: {param_path}")
        if self.machine_template and os.path.exists(self.machine_template):
            prefix = self.element or 'system'
            mparam = _generate_dpgen_machine_from_file(
                self.machine_template, prefix, is_pimd=False, nbeads=0)
            with open(os.path.join(self.dpgen_dir, 'machine.json'), 'w') as f:
                json.dump(mparam, f, indent=4)
            print(f"  machine.json written")
        self.jparam = g.jparam
        self._dpgen_gen = g
        return g

    def _write_fp_incar(self, g):
        d = _incar_dict('scf', encut=self.encut, nbands=2,
                        ele_temp_K=300)
        d['SIGMA'] = '0.5'
        d['NBANDS'] = '_AUTO_'
        content = _render_incar(d)
        content = content.replace('NBANDS = _AUTO_',
                                  '# NBANDS set by dpgen per task')
        with open(os.path.join(self.dpgen_dir, 'INCAR'), 'w') as f:
            f.write(content)

    def submit_dpgen(self):
        if self.platform == 'slurm':
            import subprocess, os
            log_path = os.path.join(self.dpgen_dir, 'dpgen.log')
            with open(log_path, 'w') as flog:
                proc = subprocess.Popen(
                    ['dpgen', 'run', 'param.json', 'machine.json'],
                    cwd=self.dpgen_dir,
                    stdout=flog, stderr=subprocess.STDOUT,
                )
            print(f"  DPGEN submitted: {self.dpgen_dir}/dpgen.log (PID {proc.pid})")
            return
        job_name = (self.element or 'system') + '_dpgen'
        g = self._dpgen_gen
        g.generate_submit(job_template_path=self.job_template, job_name=job_name,
                          platform='bh')
        g.submit()
        print(f"  DPGEN submitted: {job_name}")

    def collect_dpgen(self, collected_dir=None, set_numb=20000):
        """Harvest DPGEN results after completion.

        1. Symlink the final frozen model to ``{dpgen_dir}/frozen_model.pb``.
        2. Collect all FP-labeled structures from all iterations into
           *collected_dir* (deepmd/npy format, ready for ``extra_init_root``).

        Parameters
        ----------
        collected_dir : str, optional
            Output path.  Default: ``{dpgen_dir}/collected/``.
        set_numb : int
            Max frames per ``set.*`` subdirectory.

        Returns
        -------
        dict
            ``{'model_path': ..., 'collected_dir': ...,
               'n_iters': ..., 'n_systems': ..., 'n_frames_total': ...}``
        """
        import dpdata
        import numpy as np
        from extrempy.lazy.init_data import raw_to_set

        if collected_dir is None:
            collected_dir = os.path.join(self.dpgen_dir, 'collected')

        iter_dirs = sorted(glob.glob(os.path.join(self.dpgen_dir, 'iter.*')))
        if not iter_dirs:
            raise FileNotFoundError(
                f"No iter.* directories in {self.dpgen_dir}; "
                f"DPGEN may not have started yet.")

        print(f"-- Collect DPGEN ({len(iter_dirs)} iterations) --")

        # 1. Final frozen model (last iteration, model 000)
        model_path = None
        for name in ['frozen_model_compressed.pb', 'frozen_model.pb']:
            files = sorted(glob.glob(
                os.path.join(self.dpgen_dir, 'iter.*', '00.train', '000', name)))
            if files:
                model_path = files[-1]
                break

        if model_path:
            dst = os.path.join(self.dpgen_dir, 'frozen_model.pb')
            if os.path.lexists(dst):
                os.remove(dst)
            os.symlink(model_path, dst)
            print(f"  model: {dst}")

        # 2. Read sys labels from param.json
        param_path = os.path.join(self.dpgen_dir, 'param.json')
        labels = []
        if os.path.exists(param_path):
            with open(param_path) as f:
                jdata = json.load(f)
            labels = [
                os.path.basename(c[0]).rsplit('.', 1)[0]
                for c in jdata.get('sys_configs', [])
            ]

        # 3. Discover system indices from 02.fp/ directories
        sys_indices = set()
        for d in glob.glob(os.path.join(
                self.dpgen_dir, 'iter.*', '02.fp', 'data.*')):
            sys_indices.add(int(os.path.basename(d).split('.')[1]))

        if not sys_indices:
            print("  (no FP data found in any iteration)")
            return {
                'model_path': model_path,
                'collected_dir': collected_dir,
                'n_iters': len(iter_dirs),
                'n_systems': 0,
                'n_frames_total': 0,
            }

        os.makedirs(collected_dir, exist_ok=True)
        collected_paths = {}
        n_frames_total = 0

        for sys_idx in sorted(sys_indices):
            label = labels[sys_idx] if sys_idx < len(labels) else f'sys.{sys_idx}'
            data_dirs = sorted(glob.glob(
                os.path.join(self.dpgen_dir, 'iter.*', '02.fp',
                             f'data.{sys_idx:03d}')))
            if not data_dirs:
                continue

            ms = dpdata.MultiSystems()
            fparam_vals = []
            internal_vals = []
            last_sys = None

            for dd in data_dirs:
                subdirs = sorted([
                    os.path.join(dd, d) for d in os.listdir(dd)
                    if os.path.isdir(os.path.join(dd, d))
                ])
                sys_dirs = subdirs if subdirs else [dd]
                for sys_path in sys_dirs:
                    sys = dpdata.LabeledSystem(sys_path, fmt='deepmd/raw')
                    ms.append(sys)
                    last_sys = sys
                    fp = os.path.join(sys_path, 'fparam.raw')
                    if os.path.exists(fp):
                        fparam_vals.extend(np.loadtxt(fp).ravel().tolist())
                    ip = os.path.join(sys_path, 'internal.raw')
                    if os.path.exists(ip):
                        internal_vals.extend(np.loadtxt(ip).ravel().tolist())

            out_sub = os.path.join(collected_dir, label)
            ms.to_deepmd_raw(out_sub)

            # Flatten: MultiSystems writes to {out_sub}/{formula}/ - move raws up
            for item in os.listdir(out_sub):
                item_path = os.path.join(out_sub, item)
                if os.path.isdir(item_path):
                    for f in os.listdir(item_path):
                        shutil.move(os.path.join(item_path, f), out_sub)
                    os.rmdir(item_path)

            nf = sum(len(sub) for sub in ms.systems.values())
            print(f"  {label}: {nf} frames", end='')

            # 从 task OUTCAR 提取内能 U 和电子熵 Se
            task_outcars = sorted(glob.glob(
                os.path.join(self.dpgen_dir, 'iter.*', '02.fp',
                             f'task.{sys_idx:03d}.*', 'OUTCAR')))
            if task_outcars:
                U_list = []
                Te_list = []
                for oc in task_outcars:
                    frames = VASPReader.parse_outcar_frames(oc)
                    if frames:
                        _, U, _, Te = frames[-1]
                        U_list.append(U)
                        Te_list.append(Te)
                if len(U_list) == nf:
                    U_arr = np.array(U_list)
                    Te_arr = np.array(Te_list)
                    A_arr = np.loadtxt(os.path.join(out_sub, 'energy.raw'))
                    Se_arr = (U_arr - A_arr) / Te_arr.clip(min=1.0)
                    np.savetxt(os.path.join(out_sub, 'internal_energy.raw'), U_arr.reshape(-1, 1))
                    np.savetxt(os.path.join(out_sub, 'ele_entropy.raw'), Se_arr.reshape(-1, 1))
                    print(f' (U+Se)', end='')
                else:
                    print(f' (WARN: OUTCAR frames {len(U_list)} != {nf})', end='')

            if fparam_vals and last_sys is not None:
                natom = last_sys.get_natoms()
                fparr = np.array(fparam_vals)
                aparam = fparr.reshape(-1, 1).repeat(natom, axis=1)
                np.savetxt(os.path.join(out_sub, 'fparam.raw'), fparr)
                np.savetxt(os.path.join(out_sub, 'aparam.raw'), aparam)
                print(' (fparam)', end='')

            if internal_vals:
                np.savetxt(os.path.join(out_sub, 'internal.raw'),
                           np.array(internal_vals))
                print(' (internal)', end='')

            raw_to_set(out_sub, set_numb)
            collected_paths[label] = out_sub
            n_frames_total += nf
            print()

        print(f"  total: {n_frames_total} frames, "
              f"{len(collected_paths)} systems")
        return {
            'model_path': model_path,
            'collected_dir': collected_dir,
            'n_iters': len(iter_dirs),
            'n_systems': len(collected_paths),
            'n_frames_total': n_frames_total,
        }

    def inspect(self):
        from extrempy.dpsample import SampleSys
        return SampleSys(self.dpgen_dir, printf=True)

    def _get_tm(self):
        return 1000.0


class ElementDPBuilder(DPBuilder):
    """v1: single-element DP potential across phases from RT to 2*Tm."""

    def __init__(self, element, *,
                 mc3d_mode=None,
                 mc3d_method="pbesol-v2",
                 target_atoms=100,
                 supercell=None,
                 **kwargs):
        super().__init__(**kwargs)
        self.element = element
        self.elements = [element]
        self.work_dir = os.path.join(self.work_root, element)
        self.mc3d_mode = mc3d_mode
        self.mc3d_method = mc3d_method
        self.target_atoms = target_atoms
        self.supercell = supercell

    def _get_tm(self):
        data = ELEMENT_PHASE_DATA.get(self.element, {})
        return data.get('Tm', 1000)

    def list_mc3d_phases(self, mode=None, top_n=None):
        """
        Print MC3D phase table for this element.

        Parameters
        ----------
        mode : str or None
            None → self.mc3d_mode → 'ambient'.
        top_n : int or None
            Show only first N entries.

        Returns
        -------
        list[dict]
        """
        from extrempy.lazy.mc3d import list_phases as _list
        mode = mode or self.mc3d_mode or "ambient"
        return _list(self.element, method=self.mc3d_method, mode=mode)

    def get_phase_segments(self):
        Tm = self._get_tm()
        T_min = 200 if Tm <= 500 else 300

        if self.mc3d_mode is not None:
            from extrempy.lazy.mc3d import make_phase_segments as _mc3d
            segs = _mc3d(self.element, method=self.mc3d_method,
                         mode=self.mc3d_mode, Tm=Tm, T_min=T_min)
        else:
            segs = get_phase_segments(self.element, T_range=(T_min, None),
                                      skip_unsupported=True)

        print(f"-- Phase Segments -- {self.element} "
              f"(Tm={Tm}K, {len(segs)} phases)")
        for i, s in enumerate(segs):
            label = s['label']
            if label.endswith('-LIQ'):
                print(f"  [{i}] {label}")
                continue
            if "sg" in s:
                intl = s.get('spg_intl', '')
                sg_str = f"SG#{s['sg']}" + (f" ({intl})" if intl else "")
            else:
                sg_str = s.get('structure', '').upper()
            n_cell = s.get('n_atoms_cell')
            n_str = f"{n_cell} atoms/cell" if n_cell is not None else ""
            e_pa = s.get('energy_per_atom')
            e_str = f"{e_pa:.4f} eV/atom" if e_pa is not None else ""
            parts = [f"[{i}]", label, sg_str, n_str, e_str]
            print("  " + "  ".join(p for p in parts if p))
        return segs

    def _resolve_mc3d_poscar(self, seg, label):
        """Resolve an MC3D seg to (final_label, atoms_or_None).

        Three cases (checked in order):
          1. re-run: ``{label}-N.POSCAR`` already exists → extract N,
             return ``(label-N, None)`` (no download, no write).
          2. legacy: ``{label}.POSCAR`` exists (pre-natoms naming) →
             read N, rename to ``{label}-N.POSCAR``, return
             ``(label-N, None)``.  One-time migration; safe to remove
             after all projects migrated to the new naming.
          3. first run: download via mc3d_source, return
             ``(label-N, atoms)``; caller writes the file.

        ``atoms`` is None in cases 1-2 (file already on disk); caller
        skips ``resolve_poscar`` and just prints.
        """
        from ase.io import read

        # 1. re-run: {label}-N.POSCAR already exists
        existing = sorted(glob.glob(
            os.path.join(self.confs_dir, f'{label}-*.POSCAR')))
        if existing:
            natoms = os.path.basename(existing[0])[:-len('.POSCAR')].rsplit('-', 1)[-1]
            return f'{label}-{natoms}', None

        # 2. legacy: {label}.POSCAR exists (one-time migration)
        legacy = os.path.join(self.confs_dir, f'{label}.POSCAR')
        if os.path.exists(legacy):
            n = len(read(legacy, format='vasp'))
            os.rename(legacy, os.path.join(self.confs_dir, f'{label}-{n}.POSCAR'))
            return f'{label}-{n}', None

        # 3. first run: download + supercell
        from extrempy.structure import mc3d_source
        uid = seg.get('structure_uuid')
        if not uid:
            raise ValueError(
                f"No structure_uuid in seg for {label} "
                f"(make_phase_segments should populate it)")
        atoms = mc3d_source(uid,
                            target_atoms=self.target_atoms,
                            method=self.mc3d_method)()
        return f'{label}-{len(atoms)}', atoms

    def generate_poscars(self, segs):
        from extrempy.structure import resolve_poscar, ase_source

        self._ensure_dirs()
        print("-- POSCAR --")

        for seg in segs:
            label = seg['label']
            st = seg['structure']

            # LIQ segs: resolve_poscar handles the LIQ placeholder
            # (returns None, prints "LIQ (placeholder)").
            if label.endswith('-LIQ'):
                resolve_poscar(self.element, label,
                               confs_dir=self.confs_dir,
                               source=lambda: None)
                continue

            if st == 'mc3d':
                # Resolve label (with natoms) and get atoms if first run.
                final_label, atoms = self._resolve_mc3d_poscar(seg, label)
                seg['label'] = final_label
                natoms = final_label.rsplit('-', 1)[-1]
                if atoms is not None:
                    # First run: write the POSCAR with the final label.
                    resolve_poscar(self.element, final_label,
                                   confs_dir=self.confs_dir, source=atoms)
                print(f"  \u2713 {final_label}  ({natoms} atoms)")
                continue

            # ASE standard structures.
            if st in SUPPORTED_STRUCTURES:
                print(f'  [{label}] {st.upper()} (ASE)')
                src = ase_source(self.element, structure_type=st,
                                 supercell=self.supercell)
                resolve_poscar(self.element, label,
                               confs_dir=self.confs_dir, source=src)
            else:
                raise FileNotFoundError(
                    f"Structure type '{st}' not supported for {label}.\n"
                    f"  Use mc3d_mode='ambient' to fetch from MC3D.")


def build_all_elements(work_root, elements=None, **kwargs):
    """Batch-submit AIMD for all viable elements.

    Covers: generate_poscars + generate_init_aimd + submit_init_aimd.
    After AIMD completes, manually run collect_init_data then generate_dpgen/submit_dpgen.
    """
    if elements is None:
        elements = get_viable_elements()
    results = {}
    for el in elements:
        try:
            b = ElementDPBuilder(el, work_root=work_root, **kwargs)
            segs = b.get_phase_segments()
            b.generate_poscars(segs)
            b.generate_init_aimd(segs)
            b.submit_init_aimd()
            results[el] = {'status': 'aimd_submitted', 'builder': b}
            print(f"[{el}] AIMD submitted\n")
        except Exception as e:
            results[el] = {'status': 'error', 'msg': str(e)}
            print(f"[{el}] ERROR: {e}\n")
    return results
