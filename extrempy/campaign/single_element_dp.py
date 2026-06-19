import os
import json
import shutil
import glob

from extrempy.lazy.vasp import VASPGenerator, _incar_dict, _render_incar
from extrempy.lazy.dpgen import (DPGENGenerator,
                                 _generate_dpgen_machine_from_file,
                                 _generate_temp_list)
from extrempy.lazy.lib import (get_phase_segments, get_viable_elements,
                               ELEMENT_PHASE_DATA)
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
                 liquid_T_factor=1.8, liquid_V_scale=1.10,
                 drop_first_aimd=200, low_T_stride=50, high_T_stride=30,
                 high_T_threshold=1500):
        self.work_root = work_root
        self.potcar_lib = potcar_lib
        self.potcar_set = potcar_set
        self.machine_template = machine_template
        self.job_template = job_template
        self.encut = encut
        self.nband_scale = nband_scale
        self.nband_min = nband_min
        self.press_grid = press_grid or [1, 10, 100, 1000, 10000]
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
        self.drop_first_aimd = drop_first_aimd
        self.low_T_stride = low_T_stride
        self.high_T_stride = high_T_stride
        self.high_T_threshold = high_T_threshold
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
    def generate_init_aimd(self, segs, elements=None):
        if elements is None:
            elements = self.elements
        self.potcar_map.build(elements)  # validate ZVAL, no file write
        self._ensure_dirs()
        Tm = self._get_tm()
        self._aimd_dirs = []
        solid_label = None
        for seg in segs:
            label = seg['label']
            is_liquid = label.endswith('-LIQ')
            if is_liquid:
                poscar_ref = solid_label
                T_ref = int(self.liquid_T_factor * Tm)
            else:
                poscar_ref = label
                solid_label = label
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
            print(f"  {job_label}: T_ref={T_ref}K (1.8Tm={1.8*Tm:.0f})"
                  if is_liquid else
                  f"  {job_label}: T_ref={T_ref}K (T_core midpoint)")

    def submit_init_aimd(self, submit=True):
        for jl, wd, *_ in self._aimd_dirs:
            gen = VASPGenerator(work_path=wd, poscar_file=None)
            job_name = (self.element or 'system') + '-' + jl
            if self.platform == 'slurm':
                if self.machine_template and os.path.exists(self.machine_template):
                    gen.generate_submit(self.machine_template, job_name,
                                        platform='slurm')
                    if submit:
                        gen.submit()
                        print(f"  Submitted: {jl} ({wd})")
                    else:
                        print(f"  Script generated (not submitted): {jl} ({wd})")
                else:
                    print(f"  SKIP submit {jl}: no machine_template")
            elif self.platform == 'bh':
                if self.job_template and os.path.exists(self.job_template):
                    gen.generate_submit(self.job_template, job_name,
                                        platform='bh')
                    if submit:
                        gen.submit()
                        print(f"  Submitted: {jl} ({wd})")
                    else:
                        print(f"  Script generated (not submitted): {jl} ({wd})")
                else:
                    print(f"  SKIP submit {jl}: no job_template")
            else:
                print(f"  SKIP submit {jl}: unknown platform '{self.platform}'")

    # ---- collect init data ----
    def collect_init_data(self, segs):
        _dirs = [(l, w) for l, w, _ in self._aimd_dirs]
        results = bootstrap_init_data(
            _dirs, self.init_data_dir,
            drop_first=self.drop_first_aimd,
            low_T_stride=self.low_T_stride,
            high_T_stride=self.high_T_stride,
            high_T_threshold=self.high_T_threshold)
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
                       phase_ids=None, phase_labels=None):
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

        # Sys configs: always all phases (consistent sys numbering)
        g._set_sys_configs(set_dir=self.confs_dir,
                           prefix=self.element + '-*.POSCAR')

        g._set_model_traninig_settings(stop_batch=200000, is_ele_temp=True)
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
            print(f"  DPGEN configurations ready: {self.dpgen_dir}")
            print(f"  Run manually: cd {self.dpgen_dir} && dpgen run param.json")
            return
        job_name = (self.element or 'system') + '_dpgen'
        g = self._dpgen_gen
        g.generate_submit(job_template_path=self.job_template, job_name=job_name,
                          platform='bh')
        g.submit()
        print(f"  DPGEN submitted: {job_name}")

    def inspect(self):
        from extrempy.dpsample import SampleSys
        return SampleSys(self.dpgen_dir, printf=True)

    def _get_tm(self):
        return 1000.0


class ElementDPBuilder(DPBuilder):
    """v1: single-element DP potential across phases from RT to 2*Tm."""

    def __init__(self, element, **kwargs):
        super().__init__(**kwargs)
        self.element = element
        self.elements = [element]
        self.work_dir = os.path.join(self.work_root, element)

    def _get_tm(self):
        data = ELEMENT_PHASE_DATA.get(self.element, {})
        return data.get('Tm', 1000)

    def get_phase_segments(self):
        segs = get_phase_segments(self.element)
        Tm = self._get_tm()
        print(f"Element: {self.element}, Tm={Tm}K, {len(segs)} phase segment(s)")
        for i, s in enumerate(segs):
            tc = s['T_core']
            te = s['T_explore']
            nT = len(_generate_temp_list(te[0], te[1]))
            print(f"  [{i}] {s['label']} ({s['structure']})  "
                  f"T_core=[{tc[0]:.0f},{tc[1]:.0f}]K  "
                  f"T_explore=[{te[0]:.0f},{te[1]:.0f}]K  "
                  f"{nT} T-points")
        return segs

    def generate_poscars(self, segs):
        from extrempy.structure import generate_element_structure

        DEFAULT_SUPERCELL = {
            'fcc': (2, 2, 2),
            'bcc': (3, 3, 3),
            'hcp': (3, 3, 4),
            'diamond': (2, 2, 2),
            'dhcp': (3, 3, 2),
            'sc': (3, 3, 3),
            'bct': (3, 3, 3),
        }
        self._ensure_dirs()
        self._segs = segs
        solid_poscar_path = None
        for seg in segs:
            label = seg['label']
            st = seg['structure']
            out_path = os.path.join(self.confs_dir, label + '.POSCAR')
            if seg['label'].endswith('-LIQ'):
                print(f"  SKIP {label}: generated from AIMD CONTCAR")
                continue
            sc = DEFAULT_SUPERCELL.get(st, (3, 3, 3))
            if os.path.exists(out_path):
                print(f"  SKIP {label}: POSCAR already exists")
                solid_poscar_path = out_path
                continue
            generate_element_structure(
                element=self.element,
                output_dir=self.confs_dir,
                supercell=sc,
                structure_type=st,
                verbose=True)
            generated = glob.glob(os.path.join(
                self.confs_dir, self.element + '-' + st.upper() + '*.POSCAR'))
            if generated and os.path.basename(generated[0]) != label + '.POSCAR':
                os.rename(generated[0], out_path)
            solid_poscar_path = out_path
            print(f"  POSCAR: {out_path}")


def build_all_elements(work_root, elements=None, **kwargs):
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
            print(f"[{el}] AIMD submitted")
        except Exception as e:
            results[el] = {'status': 'error', 'msg': str(e)}
            print(f"[{el}] ERROR: {e}")
    return results
