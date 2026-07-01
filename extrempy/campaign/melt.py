import os
import glob
import subprocess
from abc import ABC, abstractmethod

import numpy as np
import ase.io

from extrempy.lazy.lammps import LAMMPSGenerator
from extrempy.lazy.lib import _get_mass_map, ELEMENT_PHASE_DATA
from extrempy.lazy.base import parse_machine_json

_TEMPLATE_DIR = os.path.join(os.path.dirname(__file__), 'templates')


class EOSCalculator(ABC):
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
        Element-internal DPGEN directory, i.e.
        ``{work_root}/{element}/dpgen`` (consistent with
        :class:`DPBuilder.dpgen_dir`).  Auto-discovers DP model by
        looking for ``frozen_model*.pb`` at the directory root, then
        under ``iter.*/00.train/000/``.  POSCAR is looked up in the
        sibling ``confs/`` directory (``{work_root}/{element}/confs/``).
    dp_model_path : str or None
        Explicit path to ``frozen_model.pb``.  Overrides all other search.
    poscar_path : str or None
        Explicit path to a POSCAR file.  Overrides all other search.
    poscar_dir : str or None
        Directory to glob ``{element}-*POSCAR``.  Falls back to
        ``{dpgen_dir}/../confs/`` (i.e. ``{work_root}/{element}/confs/``)
        if not set.
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
        self.Tm_refined = None      # set by analyze_two_phase; used by NPT/NVT
        self._two_phase_gens = []
        self._npt_gens = []
        self._traj_gens = []

    # ---- hooks (overridable) ------------------------------------------------

    @abstractmethod
    def _get_tm(self):
        """Subclass must implement: return estimated melting point (K)."""

    def _get_tm_for_run(self):
        """Return Tm for NPT/NVT: prefer refined Tm, else estimated.

        ``self.Tm_refined`` is set by :meth:`analyze_two_phase` once
        two-phase results are available; until then the estimated
        :meth:`_get_tm` is used.
        """
        if self.Tm_refined is not None:
            print(f'[{self.element}] Using refined Tm = {self.Tm_refined} K')
            return self.Tm_refined
        return self._get_tm()

    def _find_pot(self):
        """Return path to the DP frozen model.

        Priority:
          1. ``self.dp_model_path`` (explicit)
          2. ``self.dpgen_dir`` (element-internal DPGEN dir, i.e.
             ``{work_root}/{element}/dpgen`` — consistent with
             :class:`DPBuilder.dpgen_dir`) →
             a. ``{dpgen_dir}/frozen_model*.pb`` at root
             b. ``{dpgen_dir}/iter.*/00.train/000/frozen_model*.pb``
        """
        if self.dp_model_path:
            if not os.path.exists(self.dp_model_path):
                raise FileNotFoundError(
                    f'dp_model_path not found: {self.dp_model_path}')
            return self.dp_model_path

        if self.dpgen_dir:
            # 1. frozen_model*.pb at dpgen_dir root
            for name in ['frozen_model_compressed.pb', 'frozen_model.pb']:
                p = os.path.join(self.dpgen_dir, name)
                if os.path.exists(p):
                    return p
            # 2. glob iter.*/00.train/000/
            for name in ['frozen_model_compressed.pb', 'frozen_model.pb']:
                files = sorted(glob.glob(
                    os.path.join(self.dpgen_dir,
                                 'iter.*/00.train/000', name)))
                if files:
                    return files[-1]

        raise FileNotFoundError(
            f'No DP model found for {self.element}. '
            'Set dp_model_path or dpgen_dir.')

    def _find_poscar(self, role='solid_rt'):
        """Return path to a POSCAR file by role.

        Parameters
        ----------
        role : {'solid_rt', 'liquid'}
            ``solid_rt`` — room-temperature solid phase, looked up by
            ``ELEMENT_PHASE_DATA[element]['rt_structure']`` (e.g.
            ``Al-FCC.POSCAR``).
            ``liquid`` — looked up by ``{element}-LIQ.POSCAR``; if not
            found, silently falls back to the solid_rt POSCAR.

        Search order:
          1. ``self.poscar_path`` (only for ``solid_rt``)
          2. ``self.poscar_dir`` / sibling ``confs/`` of ``dpgen_dir``
             → match by explicit label
          3. solid_rt only: ASE auto-generation via :func:`resolve_poscar`
        """
        if role == 'solid_rt' and self.poscar_path:
            return self.poscar_path

        label = self._label_for_role(role)
        p = self._lookup_by_label(label)
        if p:
            return p

        if role == 'solid_rt':
            return self._autogen_solid(label)

        # role == 'liquid': silent fallback to solid_rt (no recursion —
        # _autogen_solid never calls back into 'liquid').
        return self._find_poscar('solid_rt')

    def _label_for_role(self, role):
        """Return the expected POSCAR label for a role."""
        if role == 'liquid':
            return f'{self.element}-LIQ'
        rt = ELEMENT_PHASE_DATA.get(self.element, {}).get('rt_structure')
        return f'{self.element}-{rt.upper()}' if rt else None

    def _search_dirs(self):
        """Directories to look for ``{label}.POSCAR``."""
        dirs = []
        if self.poscar_dir:
            dirs.append(self.poscar_dir)
        if self.dpgen_dir:
            # dpgen_dir is the element-internal dir; confs/ is its sibling.
            dirs.append(os.path.join(os.path.dirname(self.dpgen_dir), 'confs'))
        return dirs

    def _lookup_by_label(self, label):
        """Find ``{label}.POSCAR`` in the search dirs; return path or None."""
        if not label:
            return None
        for d in self._search_dirs():
            p = os.path.join(d, f'{label}.POSCAR')
            if os.path.exists(p):
                return p
        return None

    def _autogen_solid(self, label):
        """solid_rt fallback: ASE-generate the RT-stable phase POSCAR."""
        from extrempy.structure import resolve_poscar, ase_source
        rt = ELEMENT_PHASE_DATA.get(self.element, {}).get('rt_structure')
        if not rt:
            raise FileNotFoundError(
                f'No POSCAR found for {self.element} (role=solid_rt). '
                'Set poscar_path, poscar_dir, or dpgen_dir.')
        out_dir = os.path.join(self._element_dir, 'confs')
        return resolve_poscar(
            self.element, label,
            confs_dir=out_dir,
            source=ase_source(self.element, structure_type=rt),
            verbose=False)

    def _get_natoms(self, poscar_path, nx, ny, nz):
        atoms = ase.io.read(poscar_path, format='vasp')
        natoms_uc = len(atoms)
        natoms_total = natoms_uc * nx * ny * nz
        print(f'[{self.element}] Unit cell atoms : {natoms_uc}')
        print(f'[{self.element}] Supercell       : {nx} x {ny} x {nz}')
        print(f'[{self.element}] Total atoms     : {natoms_total}')
        return natoms_total

    def _get_npt_temps(self, npt_n=5, npt_dT=100, npt_shift=0):
        """Return NPT temperature series centered on Tm."""
        Tm = self._get_tm_for_run()
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

    def _write_sbatch(self, cfg, job_name, work_dir, needs_traj_dir=False):
        """Write job.sbatch in *work_dir*.

        Parameters
        ----------
        needs_traj_dir : bool
            If True, include ``mkdir -p traj`` in the body (for NVT
            dump jobs).  NPT/two-phase jobs do not need it.
        """
        # ---- SBATCH directive block ----
        sbatch_lines = ['#!/bin/bash', f'#SBATCH -J {job_name}']
        if cfg.get('partition'):
            sbatch_lines.append(f'#SBATCH -p {cfg["partition"]}')
        sbatch_lines.append(f'#SBATCH -N {cfg["nodes"]}')
        sbatch_lines.append(f'#SBATCH --ntasks-per-node={cfg["ntasks_per_node"]}')
        if cfg.get('wall_time'):
            sbatch_lines.append(f'#SBATCH -t {cfg["wall_time"]}')
        if cfg.get('gres'):
            sbatch_lines.append(f'#SBATCH --gres={cfg["gres"]}')

        # auto-inject --cpus-per-task when cores > ntasks
        omp_env = []
        cores = cfg.get('cores_per_node')
        ntasks = cfg.get('ntasks_per_node')
        if cores and ntasks and cores > ntasks and cores % ntasks == 0:
            cpt = cores // ntasks
            if cpt > 1:
                sbatch_lines.append(f'#SBATCH --cpus-per-task={cpt}')
                omp_env = [
                    f'export OMP_NUM_THREADS={cpt}',
                    f'export TF_INTRA_OP_PARALLELISM_THREADS={cpt}',
                    f'export TF_INTER_OP_PARALLELISM_THREADS=2',
                ]

        for flag in cfg.get('custom_flags', []):
            flag = flag.strip()
            if flag.startswith('#SBATCH') and '--job-name' not in flag:
                sbatch_lines.append(flag)

        # ---- body ----
        body = ['', f'cd {work_dir}', '']
        body += omp_env
        for src in cfg.get('source_list', []):
            body.append(src)
        for k, v in cfg.get('envs', {}).items():
            body.append(f'export {k}={v}')
        if needs_traj_dir:
            body.extend(['', 'mkdir -p traj'])

        body.append('')
        cmd = cfg['command']
        total_ranks = cfg['nodes'] * cfg['ntasks_per_node']
        if total_ranks > 1 and not any(x in cmd for x in ['mpirun', 'srun', 'mpiexec']):
            cmd = f'mpirun -np $SLURM_NTASKS {cmd}'
        body.append(cmd)

        path = os.path.join(work_dir, 'job.sbatch')
        with open(path, 'w') as f:
            f.write('\n'.join(sbatch_lines + body) + '\n')
        return path

    def _submit_slurm_job(self, gen, job_name, submit=True,
                          needs_traj_dir=False):
        """Write sbatch -> optionally submit."""
        cfg = self._resolve_slurm_config()
        self._write_sbatch(cfg, job_name, gen.work_path,
                           needs_traj_dir=needs_traj_dir)
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

        poscar = self._find_poscar('solid_rt')
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
        for gen in getattr(self, '_two_phase_gens', []):
            T_est = gen.params.get('Tm_estimate')
            job_name = f'{self.element}_{T_est}k'
            self._submit_slurm_job(gen, job_name, submit=submit,
                                   needs_traj_dir=False)

    def analyze_two_phase(self, vote_window=5):
        """Analyze ``chunk.profile`` + RDF and return coexist-based Tm.

        Each candidate-temperature directory is diagnosed via
        :func:`extrempy.campaign.chunk.diagnose_case`, which reads the
        layered Q4/Q6/density from ``chunk.profile`` (split into
        upper/lower halves) and optionally cross-checks with
        ``rdf_top.txt``.

        Parameters
        ----------
        vote_window : int
            Number of trailing blocks to vote over (adaptive: capped at
            ``n_blocks // 4`` for short trajectories).

        Returns
        -------
        dict
            Keys:
            - ``results``: ``{T_est: verdict}`` where verdict ∈
              ``solid`` / ``liquid`` / ``coexist`` / ``unknown``.
            - ``coexist_temps``: sorted list of coexist temperatures.
            - ``Tm_interval``: ``(min_coexist, max_coexist)`` or None.
            - ``Tm_refined``: median of coexist temps (int) or None.
              Also written to ``self.Tm_refined`` for downstream use.
            - ``details``: ``{T_est: dict}`` with diagnostic fields.

        Notes
        -----
        Breaking change: the old ``'partial'`` verdict is gone (it was
        never produced in practice because dump files were disabled).
        """
        from extrempy.campaign.chunk import diagnose_case

        results, details = {}, {}
        for gen in getattr(self, '_two_phase_gens', []):
            T_est = gen.params['Tm_estimate']
            diag = diagnose_case(
                chunk_path=os.path.join(gen.work_path, 'chunk.profile'),
                rdf_top_path=os.path.join(gen.work_path, 'rdf_top.txt'),
                vote_window=vote_window)
            results[T_est] = diag['verdict']
            details[T_est] = diag

        coexist_temps = sorted(t for t, v in results.items()
                               if v == 'coexist')
        Tm_interval = (min(coexist_temps), max(coexist_temps)) \
            if coexist_temps else None
        Tm_refined = (int(np.median(coexist_temps))
                      if coexist_temps else None)
        self.Tm_refined = Tm_refined

        return dict(results=results, coexist_temps=coexist_temps,
                    Tm_interval=Tm_interval, Tm_refined=Tm_refined,
                    details=details)

    # ---- Phase 2: NPT property scan ----------------------------------------

    def generate_npt(self, phases=('solid', 'liquid'),
                     supercell=(5, 5, 5),
                     npt_n=5, npt_dT=100, npt_shift=0,
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
            ``npt_shift`` are ignored.  Default series is centered on Tm
            (``npt_shift=0``), spanning roughly ``Tm ± npt_dT*npt_n/2``.
        """
        temps = latt_temp_list if latt_temp_list is not None \
            else self._get_npt_temps(npt_n, npt_dT, npt_shift)
        print(f'[{self.element}] NPT temperature series : {temps}')
        pot = self._find_pot()

        _nx, _ny, _nz = supercell

        base = self._base_params(dt=dt, pressure=pressure)
        base.update(
            equilibrate_step=equil_steps,
            nx=_nx, ny=_ny,
            nz=_nz)

        self._npt_gens = []
        for phase in phases:
            poscar = (self._find_poscar('liquid') if phase == 'liquid'
                      else self._find_poscar('solid_rt'))
            # Print atom counts once per phase (per-phase POSCAR may differ).
            self._get_natoms(poscar, _nx, _ny, _nz)
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
                params = dict(base, temperature=T, phase=phase)
                params['is_dump'] = False   # explicit & symmetric
                if output_internal:
                    model_ext = os.path.splitext(pot)[1] or '.pb'
                    params['model_name'] = f'cp{model_ext}'
                if phase == 'liquid':
                    Tm = self._get_tm_for_run()
                    params['high_temperature'] = int(
                        liquid_superheat * Tm)
                gen.update(params)
                gen.render()
                self._npt_gens.append(gen)

    def submit_npt(self, submit=True):
        for gen in getattr(self, '_npt_gens', []):
            T = gen.params.get('temperature')
            phase = gen.params.get('phase', 'solid')
            job_name = f'{self.element}_{T}k_npt_{phase}'
            self._submit_slurm_job(gen, job_name, submit=submit,
                                   needs_traj_dir=False)

    def analyze_npt(self):
        """Read each ``thermo.dat`` and return a summary DataFrame.

        Columns include temperature, phase, and per-column averages from
        the LAMMPS thermo output (energy, free_energy, ele_entropy,
        volume, density, ...).
        """
        import pandas as pd
        from extrempy.md.thermo import read_thermo_dat

        rows = []
        for gen in getattr(self, '_npt_gens', []):
            T = gen.params.get('temperature')
            phase = gen.params.get('phase', 'solid')
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
        Tm = int(self._get_tm_for_run())
        print(f'[{self.element}] NVT temperature         : {Tm} K')
        pot = self._find_pot()

        _nx, _ny, _nz = supercell

        base = self._base_params(dt=dt, pressure=pressure)
        base.update(
            equilibrate_step=equil_steps,
            dump_freq=dump_freq,
            nx=_nx, ny=_ny,
            nz=_nz)

        self._traj_gens = []
        for phase in phases:
            poscar = (self._find_poscar('liquid') if phase == 'liquid'
                      else self._find_poscar('solid_rt'))
            # Print atom counts once per phase (per-phase POSCAR may differ).
            self._get_natoms(poscar, _nx, _ny, _nz)
            template = ('nvt-liquid-traj.j2' if phase == 'liquid'
                        else 'nvt-solid-traj.j2')
            work_dir = os.path.join(self.traj_dir, f'{Tm}k_{phase}')
            os.makedirs(work_dir, exist_ok=True)
            gen = self._build_gen(work_dir, template)
            gen.get_files(poscar_path=poscar, pot_path=pot)
            params = dict(base, temperature=Tm, phase=phase)
            if phase == 'liquid':
                params['high_temperature'] = int(
                    liquid_superheat * Tm)
            gen.update(params)
            gen.render()
            self._traj_gens.append(gen)

    def submit_nvt_traj(self, submit=True):
        for gen in getattr(self, '_traj_gens', []):
            T = gen.params['temperature']
            phase = gen.params.get('phase', 'solid')
            job_name = f'{self.element}_{T}k_nvt_{phase}'
            self._submit_slurm_job(gen, job_name, submit=submit,
                                   needs_traj_dir=True)

    # ---- full pipeline ----------------------------------------------------

    def run_two_phase(self, submit=True):
        """Phase 1 only: generate + submit two-phase jobs."""
        self.generate_two_phase()
        self.submit_two_phase(submit=submit)

    def run_property_scans(self, submit=True):
        """Phase 2+3: generate + submit NPT + NVT.

        Uses ``self.Tm_refined`` if available (set by
        :meth:`analyze_two_phase`), otherwise falls back to the
        estimated :meth:`_get_tm`.
        """
        self.generate_npt()
        self.submit_npt(submit=submit)
        self.generate_nvt_traj()
        self.submit_nvt_traj(submit=submit)

    def run_all(self, submit=True):
        """Convenience: two-phase + property scans using ESTIMATED Tm.

        NOTE: This does **not** call :meth:`analyze_two_phase`, so
        ``Tm_refined`` is not set — NPT/NVT use the estimated
        :meth:`_get_tm`.  For the refined-Tm workflow::

            calc.run_two_phase(submit=True)
            # ... wait for jobs, then:
            calc.analyze_two_phase()            # sets Tm_refined
            calc.run_property_scans(submit=True)
        """
        self.run_two_phase(submit=submit)
        self.run_property_scans(submit=submit)


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

    Notes
    -----
    ``element`` is positional-only in :class:`ElementEOSCalculator`;
    do not pass ``element=`` via ``**kwargs``.
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
