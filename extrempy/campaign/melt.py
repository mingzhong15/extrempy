import os
import glob
import subprocess
import numpy as np
import ase.io

from extrempy.lazy.lammps import LAMMPSGenerator
from extrempy.lazy.lib import _get_mass_map, ELEMENT_PHASE_DATA
from extrempy.lazy.base import parse_machine_json

_TEMPLATE_DIR = os.path.join(os.path.dirname(__file__), 'templates')


def recommend_mpi_layout(natoms, nodes, ntasks_per_node):
    """Check atoms-per-task ratio and suggest MPI layout if too low.

    Returns ``(ok, atoms_per_task, suggestion)``.
    """
    ntasks = nodes * ntasks_per_node
    apt = natoms / max(ntasks, 1)
    if apt >= 500:
        return True, apt, ''
    ideal_tasks = max(1, natoms // 800)
    ideal_nodes = max(1, (ideal_tasks + ntasks_per_node - 1) // ntasks_per_node)
    return False, apt, (
        f'atoms/task={apt:.0f} < 500, '
        f'suggest ~{ideal_tasks} MPI tasks ({ideal_nodes} node(s))')


class EOSCalculator:
    """LAMMPS-based Equation-of-State and melt determination pipeline.

    The constructor captures only infrastructure configuration (paths, Slurm
    resources, template directory).  All simulation parameters are passed
    directly to the ``generate_*`` step methods.

    Directory layout per element::

        {work_root}/{element}/
            melt/{T}k/        # two-phase result (one dir per candidate T)
            npt/{T}k_{phase}/ # NPT property scan
            traj/{T}k_{phase}/# NVT trajectory dump

    Parameters
    ----------
    work_root : str
        Root directory for all element data.
    dpgen_dir : str or None
        DPGEN project root.  Auto-discovers DP model from
        ``{dpgen_dir}/{element}/dpgen/iter.*/00.train/000/`` and POSCAR
        from ``{dpgen_dir}/{element}/confs/*.POSCAR``.
    dp_model_path : str or None
        Explicit path to ``frozen_model.pb``.  Overrides all other search.
    poscar_path : str or None
        Explicit path to a POSCAR file.  Overrides all other search.
    poscar_dir : str or None
        Directory to glob ``{element}-*POSCAR``.  Falls back to
        ``dpgen_dir/{element}/confs/`` if not set.
    machine_template : str or None
        dpgen-style ``machine.json`` — primary source for slurm resources
        (partition, nodes, ``source_list`` for env setup, etc.).
        Constructor parameters with a non-``None`` value override these.
    partition : str or None
        Slurm partition (overrides ``machine_template``).
    nodes : int or None
        Number of nodes (overrides ``machine_template``).
    ntasks_per_node : int or None
        Tasks (MPI ranks) per node (overrides ``machine_template``).
    wall_time : str or None
        Wall time in ``'HH:MM:SS'`` format (overrides ``machine_template``).
    gres : str or None
        Generic resource scheduling, e.g. ``'gpu:1'`` (overrides ``machine_template``).
    lmp_command : str or None
        LAMMPS command, e.g. ``'lmp -in run.in > log.run'``.
        Defaults to ``'lmp -in run.in > log.run'`` if neither
        constructor nor ``machine_template`` specify it.
    template_dir : str or None
        Jinja2 template directory.  Defaults to the built-in templates
        shipped with the package (``extrempy/campaign/templates/``).
    """

    def __init__(self, work_root,
                 # model & POSCAR (dpgen_dir auto-links, explicit paths override)
                 dpgen_dir=None,
                 dp_model_path=None,
                 poscar_path=None,
                 poscar_dir=None,

                 # slurm resources (machine_template provides defaults)
                 machine_template=None,
                 partition=None,
                 nodes=None,
                 ntasks_per_node=None,
                 wall_time=None,
                 gres=None,
                 lmp_command=None,

                 # templates (built-in by default)
                 template_dir=None):

        self.work_root = work_root

        # model & POSCAR sources
        self.dpgen_dir = dpgen_dir
        self.dp_model_path = dp_model_path
        self.poscar_path = poscar_path
        self.poscar_dir = poscar_dir

        # slurm
        self.machine_template = machine_template
        self.partition = partition
        self.nodes = nodes
        self.ntasks_per_node = ntasks_per_node
        self.wall_time = wall_time
        self.gres = gres
        self.lmp_command = lmp_command

        # templates
        self.template_dir = template_dir if template_dir else _TEMPLATE_DIR

        # state (set by subclass or batch runner)
        self.element = None

    # ---- hooks (overridable) ------------------------------------------------

    def _get_tm(self):
        raise NotImplementedError

    def _find_pot(self):
        """Return path to the DP frozen model.

        Priority:
          1. ``self.dp_model_path`` (explicit)
          2. ``self.dpgen_dir`` → auto-detect DPGEN output
        """
        if self.dp_model_path:
            if not os.path.exists(self.dp_model_path):
                raise FileNotFoundError(
                    f'dp_model_path not found: {self.dp_model_path}')
            return self.dp_model_path

        if self.dpgen_dir:
            for name in ['frozen_model_compressed.pb', 'frozen_model.pb']:
                pat = os.path.join(
                    self.dpgen_dir,
                    f'{self.element}/dpgen/iter.*/00.train/000/',
                    name)
                files = sorted(glob.glob(pat))
                if files:
                    return files[-1]

        raise FileNotFoundError(
            f'No DP model found for {self.element}. '
            'Set dp_model_path or dpgen_dir.')

    def _find_poscar(self, idx=0):
        """Return path to a POSCAR file.

        Priority:
          1. ``self.poscar_path`` (explicit file)
          2. ``self.poscar_dir`` → glob ``{element}-*POSCAR``
          3. ``self.dpgen_dir`` → glob ``{dpgen_dir}/{element}/confs/*.POSCAR``
          4. auto-generate from ``ELEMENT_PHASE_DATA`` (primitive cell)
        """
        if self.poscar_path:
            return self.poscar_path

        if self.poscar_dir:
            pat = os.path.join(self.poscar_dir, f'{self.element}-*POSCAR')
            files = sorted(glob.glob(pat))
            if files:
                return files[idx]

        if self.dpgen_dir:
            pat = os.path.join(self.dpgen_dir, self.element, 'confs', '*.POSCAR')
            files = sorted(glob.glob(pat))
            if files:
                return files[idx]

        # 4th priority: auto-generate from local ELEMENT_PHASE_DATA
        from extrempy.structure import _get_local_candidates, save_structures
        cands = _get_local_candidates(self.element)
        if cands:
            out_dir = os.path.join(self._element_dir, 'confs')
            saved = save_structures(
                cands[:1], out_dir, supercell=(1, 1, 1))
            return next(iter(saved.values()))

        raise FileNotFoundError(
            f'No POSCAR found for {self.element}. '
            'Set poscar_path, poscar_dir, or dpgen_dir.')

    def _get_natoms(self, poscar_path, nx, ny, nz):
        atoms = ase.io.read(poscar_path, format='vasp')
        natoms_uc = len(atoms)
        natoms_total = natoms_uc * nx * ny * nz
        print(f'[{self.element}] Unit cell atoms : {natoms_uc}')
        print(f'[{self.element}] Supercell       : {nx} x {ny} x {nz}')
        print(f'[{self.element}] Total atoms     : {natoms_total}')
        return natoms_total

    def _get_npt_temps(self, npt_n=5, npt_dT=100, npt_shift=-600):
        """Return NPT temperature series."""
        Tm = self._get_tm()
        half = npt_n // 2
        offsets = np.arange(npt_n) - half
        temps = Tm + offsets * npt_dT + npt_shift
        return [int(t) for t in temps if t >= 250]

    # ---- paths --------------------------------------------------------------

    @property
    def _element_dir(self):
        return os.path.join(self.work_root, self.element)

    @property
    def melt_dir(self):
        return os.path.join(self._element_dir, 'melt')

    @property
    def npt_dir(self):
        return os.path.join(self._element_dir, 'npt')

    @property
    def traj_dir(self):
        return os.path.join(self._element_dir, 'traj')

    # ---- helpers ------------------------------------------------------------

    def _build_gen(self, work_dir, template_file):
        return LAMMPSGenerator(
            work_path=work_dir,
            template_path=self.template_dir,
            template_file=template_file)

    def _base_params(self, dt=0.001, pressure=0.0001):
        return dict(
            pressure=pressure,
            dt=dt,
            elements=[dict(
                id=1, name=self.element,
                mass=_get_mass_map([self.element])[0])])

    def _resolve_slurm_config(self):
        """Build Slurm config with priority: machine_template -> constructor -> fallback.

        Returns dict with keys: partition, nodes, ntasks_per_node,
        cores_per_node, wall_time, gres, command, source_list,
        custom_flags, envs.
        """
        # ---- Phase 1: hardcoded fallback defaults ----
        cfg = dict(
            partition=None,
            nodes=1,
            ntasks_per_node=32,
            wall_time='24:00:00',
            gres=None,
            command='lmp -in run.in > log.run',
            source_list=[],
            custom_flags=[],
            envs={},
        )

        # ---- Phase 2: machine_template fills in higher-priority defaults ----
        if self.machine_template:
            tmpl = parse_machine_json(
                self.machine_template, section='model_devi')
            for k, v in tmpl.items():
                if k == 'command':          # env activation goes via source_list
                    continue                # lmp_command comes from fallback or user
                if v:                       # non-empty values override phase-1
                    cfg[k] = v

        # ---- Phase 3: constructor explicit values (non-None) override ----
        overrides = [
            ('partition', 'partition'),
            ('nodes', 'nodes'),
            ('ntasks_per_node', 'ntasks_per_node'),
            ('wall_time', 'wall_time'),
            ('gres', 'gres'),
            ('lmp_command', 'command'),
        ]
        for attr, key in overrides:
            val = getattr(self, attr)
            if val is not None:
                cfg[key] = val

        return cfg

    def _write_sbatch(self, cfg, job_name, work_dir):
        """Write job.sbatch in *work_dir*."""
        lines = ['#!/bin/bash']
        lines.append(f'#SBATCH -J {job_name}')
        if cfg.get('partition'):
            lines.append(f'#SBATCH -p {cfg["partition"]}')
        lines.append(f'#SBATCH -N {cfg["nodes"]}')
        lines.append(f'#SBATCH --ntasks-per-node={cfg["ntasks_per_node"]}')
        if cfg.get('wall_time'):
            lines.append(f'#SBATCH -t {cfg["wall_time"]}')
        if cfg.get('gres'):
            lines.append(f'#SBATCH --gres={cfg["gres"]}')
        for flag in cfg.get('custom_flags', []):
            flag = flag.strip()
            if flag.startswith('#SBATCH') and '--job-name' not in flag:
                lines.append(flag)
        lines.append('')
        lines.append(f'cd {work_dir}')

        # auto-inject --cpus-per-task + OMP_NUM_THREADS when user reduces ntasks
        cores = cfg.get('cores_per_node')
        ntasks = cfg.get('ntasks_per_node')
        if cores and ntasks and cores > ntasks and cores % ntasks == 0:
            cpt = cores // ntasks
            if cpt > 1:
                lines.insert(4, f'#SBATCH --cpus-per-task={cpt}')
                lines.append(f'export OMP_NUM_THREADS={cpt}')
                lines.append(f'export TF_INTRA_OP_PARALLELISM_THREADS={cpt}')
                lines.append(f'export TF_INTER_OP_PARALLELISM_THREADS=2')

        # environment setup from machine_template
        for src in cfg.get('source_list', []):
            lines.append(src)
        for k, v in cfg.get('envs', {}).items():
            lines.append(f'export {k}={v}')

        lines.append('')
        lines.append('mkdir -p traj')
        cmd = cfg['command']
        total_ranks = cfg['nodes'] * cfg['ntasks_per_node']
        if total_ranks > 1 and not any(x in cmd for x in ['mpirun', 'srun', 'mpiexec']):
            cmd = f'mpirun -np $SLURM_NTASKS {cmd}'
        lines.append(cmd)

        path = os.path.join(work_dir, 'job.sbatch')
        with open(path, 'w') as f:
            f.write('\n'.join(lines) + '\n')
        return path

    def _submit_slurm_job(self, gen, job_name, submit=True):
        """Write sbatch -> optionally submit."""
        cfg = self._resolve_slurm_config()
        self._write_sbatch(cfg, job_name, gen.work_path)
        if submit:
            subprocess.run(
                ['sbatch', 'job.sbatch'],
                cwd=gen.work_path,
                check=True)

    # ---- Phase 1: two-phase melt determination ------------------------------

    def generate_two_phase(self,
                           supercell=(5, 5, 10),
                           two_phase_temps=None,
                           two_phase_delta=100, two_phase_count=3,
                           equil_steps=100000, heat_steps=10000,
                           dt=0.001, pressure=0.0001, Q_cutoff=3.0,
                           liquid_superheat=1.9):
        """Generate two-phase LAMMPS inputs at candidate temperatures."""
        Tm = self._get_tm()

        if two_phase_temps is not None:
            temps = two_phase_temps
        else:
            half = two_phase_count // 2
            offsets = [i * two_phase_delta for i in range(-half, half + 1)]
            temps = [int(Tm + off) for off in offsets]

        print(f'[{self.element}] Melting point (Tm) : {Tm} K')
        print(f'[{self.element}] Two-phase temps   : {temps}')

        poscar = self._find_poscar()
        pot = self._find_pot()

        _nx, _ny, _nz = supercell
        self._get_natoms(poscar, _nx, _ny, _nz)

        base = self._base_params(dt=dt, pressure=pressure)
        base.update(
            Q_cutoff=Q_cutoff,
            heating_step=heat_steps,
            equilibrate_step=equil_steps,
            nx=_nx, ny=_ny,
            nz=_nz)

        self._two_phase_gens = []
        for T_est in temps:
            work_dir = os.path.join(self.melt_dir, f'{T_est}k')
            os.makedirs(work_dir, exist_ok=True)
            gen = self._build_gen(work_dir, 'two-phase.j2')
            gen.get_files(poscar_path=poscar, pot_path=pot)
            gen.update(dict(
                base,
                Tm_estimate=T_est,
                T_superheat=int(liquid_superheat * Tm)))
            gen.render()
            self._two_phase_gens.append(gen)

    def submit_two_phase(self, submit=True):
        for gen in self._two_phase_gens:
            T_est = gen.params.get('Tm_estimate')
            job_name = f'{self.element}_{T_est}k'
            self._submit_slurm_job(gen, job_name, submit=submit)

    def analyze_two_phase(self):
        """Analyze two-phase dump files and return detected Tm interval."""
        from extrempy.md.traj import read_dump_file, calculate_q4_q6

        results = {}
        for gen in self._two_phase_gens:
            T_est = gen.params['Tm_estimate']
            dump_dir = os.path.join(gen.work_path, 'traj')
            dump_files = sorted(glob.glob(os.path.join(dump_dir, 'dump.*')))
            if not dump_files:
                dump_files = sorted(glob.glob(os.path.join(
                    gen.work_path, 'dump.*')))

            verdict = 'unknown'
            for df in dump_files[-3:]:
                atoms = read_dump_file(df)
                if atoms is None:
                    continue
                q4, q6, *_ = calculate_q4_q6(atoms)
                if q4 > 0.1 and q6 > 0.3:
                    verdict = 'solid'
                elif q4 > 0.05:
                    verdict = 'partial'
                else:
                    verdict = 'liquid'

            results[T_est] = verdict

        solid_temps = [t for t, v in results.items() if v == 'solid']
        liquid_temps = [t for t, v in results.items() if v == 'liquid']
        Tm_interval = (max(solid_temps), min(liquid_temps)) \
            if solid_temps and liquid_temps else None
        return dict(results=results, Tm_interval=Tm_interval)

    # ---- Phase 2: NPT property scan ----------------------------------------

    def generate_npt(self, phases=('solid', 'liquid'),
                     supercell=(5, 5, 5),
                     npt_n=5, npt_dT=100, npt_shift=-600,
                     liquid_superheat=1.9,
                     equil_steps=100000, dt=0.001, pressure=0.0001,
                     output_internal=False,
                     latt_temp_list=None):
        """Generate NPT LAMMPS inputs at multiple temperatures.

        Parameters
        ----------
        output_internal : bool
            When True, use entropy-enabled templates that output internal
            energy (U = F + T*S) via ``out_internal_energy`` pair_style
            keyword, plus ``ele_entropy`` and ``free_energy`` computes.
        latt_temp_list : list of int, optional
            Explicit temperature series.  When given, ``npt_n`` / ``npt_dT`` /
            ``npt_shift`` are ignored.
        """
        temps = latt_temp_list if latt_temp_list is not None \
            else self._get_npt_temps(npt_n, npt_dT, npt_shift)
        print(f'[{self.element}] NPT temperature series : {temps}')
        poscar = self._find_poscar(
            -1 if 'liquid' in phases and len(phases) > 1 else 0)
        pot = self._find_pot()

        _nx, _ny, _nz = supercell
        self._get_natoms(poscar, _nx, _ny, _nz)

        base = self._base_params(dt=dt, pressure=pressure)
        base.update(
            equilibrate_step=equil_steps,
            nx=_nx, ny=_ny,
            nz=_nz)

        self._npt_gens = []
        for phase in phases:
            if output_internal:
                template = ('npt-entropy-liquid.j2' if phase == 'liquid'
                            else 'npt-entropy-solid.j2')
            else:
                template = ('npt-liquid.j2' if phase == 'liquid'
                            else 'npt-solid.j2')
            for T in temps:
                work_dir = os.path.join(
                    self.npt_dir, f'{T}k_{phase}')
                os.makedirs(work_dir, exist_ok=True)
                gen = self._build_gen(work_dir, template)
                gen.get_files(poscar_path=poscar, pot_path=pot)
                params = dict(base, temperature=T)
                if output_internal:
                    model_ext = os.path.splitext(pot)[1] or '.pb'
                    params['model_name'] = f'cp{model_ext}'
                if phase == 'liquid':
                    Tm = self._get_tm()
                    params['high_temperature'] = int(
                        liquid_superheat * Tm)
                    params['is_dump'] = False
                gen.update(params)
                gen.render()
                self._npt_gens.append(gen)

    def submit_npt(self, submit=True):
        for gen in self._npt_gens:
            T = gen.params.get('temperature')
            phase = 'liquid' if 'high_temperature' in gen.params else 'solid'
            job_name = f'{self.element}_{T}k_npt_{phase}'
            self._submit_slurm_job(gen, job_name, submit=submit)

    def analyze_npt(self):
        """Read each ``thermo.dat`` and return a summary DataFrame.

        Columns include temperature, phase, and per-column averages from
        the LAMMPS thermo output (energy, free\_energy, ele\_entropy,
        volume, density, …).
        """
        import pandas as pd
        from extrempy.md.thermo import read_thermo_dat

        rows = []
        for gen in getattr(self, '_npt_gens', []):
            T = gen.params.get('temperature')
            phase = 'liquid' if 'high_temperature' in gen.params else 'solid'
            thermo_file = os.path.join(gen.work_path, 'thermo.dat')
            if not os.path.exists(thermo_file):
                continue
            _, averages, _ = read_thermo_dat(thermo_file)
            if averages is None:
                continue
            row = {'temperature': T, 'phase': phase}
            for col, stats in averages.items():
                row[f'{col}_mean'] = stats['mean']
            rows.append(row)
        df = pd.DataFrame(rows)
        if not df.empty:
            df = df.sort_values('temperature').reset_index(drop=True)
        return df

    # ---- Phase 3: NVT trajectory ------------------------------------------

    def generate_nvt_traj(self, phases=('solid', 'liquid'),
                          supercell=(5, 5, 5),
                          liquid_superheat=1.9,
                          equil_steps=100000, dump_freq=10,
                          dt=0.001, pressure=0.0001):
        """Generate NVT trajectory LAMMPS inputs."""
        Tm = int(self._get_tm())
        print(f'[{self.element}] NVT temperature         : {Tm} K')
        poscar = self._find_poscar(
            -1 if 'liquid' in phases and len(phases) > 1 else 0)
        pot = self._find_pot()

        _nx, _ny, _nz = supercell
        self._get_natoms(poscar, _nx, _ny, _nz)

        base = self._base_params(dt=dt, pressure=pressure)
        base.update(
            equilibrate_step=equil_steps,
            dump_freq=dump_freq,
            nx=_nx, ny=_ny,
            nz=_nz)

        self._traj_gens = []
        for phase in phases:
            template = ('nvt-liquid-traj.j2' if phase == 'liquid'
                        else 'nvt-solid-traj.j2')
            work_dir = os.path.join(self.traj_dir, f'{Tm}k_{phase}')
            os.makedirs(work_dir, exist_ok=True)
            gen = self._build_gen(work_dir, template)
            gen.get_files(poscar_path=poscar, pot_path=pot)
            params = dict(base, temperature=Tm)
            if phase == 'liquid':
                params['high_temperature'] = int(
                    liquid_superheat * Tm)
            gen.update(params)
            gen.render()
            self._traj_gens.append(gen)

    def submit_nvt_traj(self, submit=True):
        for gen in self._traj_gens:
            T = gen.params['temperature']
            phase = 'liquid' if 'high_temperature' in gen.params else 'solid'
            job_name = f'{self.element}_{T}k_nvt_{phase}'
            self._submit_slurm_job(gen, job_name, submit=submit)

    # ---- full pipeline ----------------------------------------------------

    def run_all(self, submit=True):
        """Convenience: run the full generate -> submit pipeline."""
        self.generate_two_phase()
        self.submit_two_phase(submit=submit)
        self.generate_npt()
        self.submit_npt(submit=submit)
        self.generate_nvt_traj()
        self.submit_nvt_traj(submit=submit)


class ElementEOSCalculator(EOSCalculator):
    """EOSCalculator bound to a single element.

    Examples
    --------
    >>> calc = ElementEOSCalculator('Al',
    ...     work_root='/share/zeng/metals/dpmd',
    ...     dpgen_dir='/share/zeng/metals/sample',
    ...     machine_template='~/template/dpgen-machine.json')
    >>> calc.run_all(submit=False)
    """

    def __init__(self, element, **kwargs):
        super().__init__(**kwargs)
        self.element = element

    def _get_tm(self):
        data = ELEMENT_PHASE_DATA.get(self.element, {})
        return data.get('Tm', 1000)


def run_eos_all(elements, work_root, **kwargs):
    """Batch EOSCalculator for a list of elements.

    Parameters
    ----------
    elements : list of str
    work_root : str
    **kwargs
        Forwarded to each ``ElementEOSCalculator``.

    Returns
    -------
    dict
        ``{element: 'generated' | error_message}``
    """
    results = {}
    for el in elements:
        try:
            calc = ElementEOSCalculator(
                el, work_root=work_root, **kwargs)
            calc.run_all(submit=False)
            results[el] = 'generated'
        except Exception as e:
            results[el] = str(e)
    return results
