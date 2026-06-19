import os
import glob
import numpy as np

from extrempy.lazy.lammps import LAMMPSGenerator
from extrempy.lazy.lib import _get_mass_map, ELEMENT_PHASE_DATA
from extrempy.constant import kb, NA, m2A, s2ps

_TEMPLATE_DIR = os.path.join(os.path.dirname(__file__), 'templates')


class EOSCalculator:
    """LAMMPS-based Equation-of-State and melt determination pipeline.

    Mirrors DPBuilder design pattern:
      * Constructor captures ALL configuration.
      * Generate / submit / analyze are separate step methods.
      * Hook methods (_get_tm, _find_pot, …) let subclasses customize.
      * Element-level subclass (ElementEOSCalculator) for single-element jobs.

    Directory layout per element::

        {work_root}/{element}/
            melt/{T}k/        # two-phase result (one dir per candidate T)
            npt/{T}k_{phase}/ # NPT property scan
            traj/{T}k_{phase}/# NVT trajectory dump

    Parameters
    ----------
    work_root : str
        Root directory for all element data.
    poscar_dir : str
        Directory containing POSCAR files (``{element}-*POSCAR``).
    pot_root_dir : str
        Root directory containing ``{element}_sample/iter.*/00.train/000/``
        with ``frozen_model.pb``.
    template_dir : str or None
        Jinja2 template directory.  Defaults to the built-in templates
        shipped with the package (``extrempy/campaign/templates/``).
    job_template : str or None
        Path to a JSON job template for ``InputGenerator.generate_submit``.
    platform : str
        Submission platform (``'bh'`` or ``'slurm'``).
    """

    def __init__(self, work_root, poscar_dir, pot_root_dir,
                 template_dir=None, job_template=None, platform='bh',
                 two_phase_nz=10, npt_n=5, npt_dT=100, npt_shift=-600,
                 two_phase_frac=0.10, two_phase_shift=300,
                 supercell=(5, 5, 5),
                 liquid_superheat=1.9,
                 equil_steps=100000, heat_steps=10000,
                 dump_freq=10, dt=0.001, Q_cutoff=3.0, pressure=0.0001):

        self.work_root = work_root
        self.poscar_dir = poscar_dir
        self.pot_root_dir = pot_root_dir
        self.template_dir = template_dir if template_dir else _TEMPLATE_DIR
        self.job_template = job_template
        self.platform = platform

        # temperature settings
        self.two_phase_frac = two_phase_frac
        self.two_phase_shift = two_phase_shift
        self.two_phase_nz = two_phase_nz
        self.npt_n = npt_n
        self.npt_dT = npt_dT
        self.npt_shift = npt_shift
        self.supercell = supercell
        self.liquid_superheat = liquid_superheat

        # simulation parameters
        self.equil_steps = equil_steps
        self.heat_steps = heat_steps
        self.dump_freq = dump_freq
        self.dt = dt
        self.Q_cutoff = Q_cutoff
        self.pressure = pressure

        # state
        self.element = None

    # ---- hooks (overridable) ------------------------------------------------

    def _get_tm(self):
        raise NotImplementedError

    def _find_pot(self):
        """Return path to the latest compressed DP frozen model."""
        for compressed in [True, False]:
            pat = os.path.join(
                self.pot_root_dir,
                f'{self.element}_sample/iter.00*/00.train/000/'
                f'{"frozen_model_compressed.pb" if compressed else "frozen_model.pb"}')
            files = sorted(glob.glob(pat))
            if files:
                return files[-1]
        raise FileNotFoundError(f'No DP potential found for {self.element}')

    def _find_poscar(self, idx=0):
        """Return path to a POSCAR file for *element*."""
        pat = os.path.join(self.poscar_dir, f'{self.element}-*POSCAR')
        files = sorted(glob.glob(pat))
        if not files:
            raise FileNotFoundError(f'No POSCAR found for {self.element}')
        return files[idx]

    def _get_two_phase_temps(self):
        """Return candidate temperatures for two-phase runs."""
        Tm = self._get_tm()
        dT = (int(Tm * self.two_phase_frac / self.npt_dT) + 1) * self.npt_dT
        return [int(Tm + off + self.two_phase_shift)
                for off in [-dT, 0, dT]]

    def _get_npt_temps(self):
        """Return NPT temperature series."""
        Tm = self._get_tm()
        half = self.npt_n // 2
        offsets = np.arange(self.npt_n) - half
        temps = Tm + offsets * self.npt_dT + self.npt_shift
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

    def _base_params(self):
        return dict(
            pressure=self.pressure,
            dt=self.dt,
            elements=[dict(
                id=1, name=self.element,
                mass=_get_mass_map([self.element])[0])])

    # ---- Phase 1: two-phase melt determination ------------------------------

    def generate_two_phase(self):
        """Generate two-phase LAMMPS inputs at candidate temperatures."""
        Tm = self._get_tm()
        temps = self._get_two_phase_temps()
        poscar = self._find_poscar()
        pot = self._find_pot()

        base = self._base_params()
        base.update(
            Q_cutoff=self.Q_cutoff,
            heating_step=self.heat_steps,
            equilibrate_step=self.equil_steps,
            nx=self.supercell[0], ny=self.supercell[1],
            nz=self.two_phase_nz)

        self._two_phase_gens = []
        for T_est in temps:
            work_dir = os.path.join(self.melt_dir, f'{T_est}k')
            os.makedirs(work_dir, exist_ok=True)
            gen = self._build_gen(work_dir, 'two-phase.j2')
            gen.get_files(poscar_path=poscar, pot_path=pot)
            gen.update(dict(
                base,
                Tm_estimate=T_est,
                T_superheat=int(self.liquid_superheat * Tm)))
            gen.render()
            self._two_phase_gens.append(gen)

    def submit_two_phase(self, submit=True):
        for gen in self._two_phase_gens:
            T_est = gen.params.get('Tm_estimate')
            gen.generate_submit(
                self.job_template,
                job_name=f'{self.element}_{T_est}k',
                platform=self.platform)
            if submit:
                gen.submit()

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
            for df in dump_files[-3:]:   # last few dumps
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

    def generate_npt(self, phases=('solid', 'liquid')):
        """Generate NPT LAMMPS inputs at multiple temperatures."""
        temps = self._get_npt_temps()
        poscar = self._find_poscar(
            -1 if 'liquid' in phases and len(phases) > 1 else 0)
        pot = self._find_pot()

        base = self._base_params()
        base.update(
            equilibrate_step=self.equil_steps,
            nx=self.supercell[0], ny=self.supercell[1],
            nz=self.supercell[2])

        self._npt_gens = []
        for phase in phases:
            template = 'npt-liquid.j2' if phase == 'liquid' else 'npt-solid.j2'
            for T in temps:
                work_dir = os.path.join(
                    self.npt_dir, f'{T}k_{phase}')
                os.makedirs(work_dir, exist_ok=True)
                gen = self._build_gen(work_dir, template)
                gen.get_files(poscar_path=poscar, pot_path=pot)
                params = dict(base, temperature=T)
                if phase == 'liquid':
                    Tm = self._get_tm()
                    params['high_temperature'] = int(
                        self.liquid_superheat * Tm)
                    params['is_dump'] = False
                gen.update(params)
                gen.render()
                self._npt_gens.append(gen)

    def submit_npt(self, submit=True):
        for gen in self._npt_gens:
            T = gen.params.get('temperature')
            phase = 'liquid' if 'high_temperature' in gen.params else 'solid'
            gen.generate_submit(
                self.job_template,
                job_name=f'{self.element}_{T}k_npt_{phase}',
                platform=self.platform)
            if submit:
                gen.submit()

    def analyze_npt(self):
        """Process NPT directories and return summary DataFrame."""
        from extrempy.md.thermo import process_npt_directories
        summary, _ = process_npt_directories(
            base_dir=self.npt_dir,
            element=self.element,
            phases=['solid', 'liquid'])
        return summary

    # ---- Phase 3: NVT trajectory ------------------------------------------

    def generate_nvt_traj(self, phases=('solid', 'liquid')):
        """Generate NVT trajectory LAMMPS inputs."""
        Tm = int(self._get_tm())
        poscar = self._find_poscar(
            -1 if 'liquid' in phases and len(phases) > 1 else 0)
        pot = self._find_pot()

        base = self._base_params()
        base.update(
            equilibrate_step=self.equil_steps,
            dump_freq=self.dump_freq,
            nx=self.supercell[0], ny=self.supercell[1],
            nz=self.supercell[2])

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
                    self.liquid_superheat * Tm)
            gen.update(params)
            gen.render()
            self._traj_gens.append(gen)

    def submit_nvt_traj(self, submit=True):
        for gen in self._traj_gens:
            T = gen.params['temperature']
            phase = 'liquid' if 'high_temperature' in gen.params else 'solid'
            gen.generate_submit(
                self.job_template,
                job_name=f'{self.element}_{T}k_nvt_{phase}',
                platform=self.platform)
            if submit:
                gen.submit()

    # ---- full pipeline ----------------------------------------------------

    def run_all(self, submit=True):
        """Convenience: run the full generate → submit pipeline."""
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
    >>> calc = ElementEOSCalculator(
    ...     'Al', work_root='/tmp/eos', poscar_dir='/poscars',
    ...     pot_root_dir='/pot')
    >>> calc.run_all(submit=False)
    """

    def __init__(self, element, **kwargs):
        super().__init__(**kwargs)
        self.element = element

    def _get_tm(self):
        data = ELEMENT_PHASE_DATA.get(self.element, {})
        return data.get('Tm', 1000)


def run_eos_all(elements, work_root, poscar_dir, pot_root_dir, **kwargs):
    """Batch EOSCalculator for a list of elements.

    Parameters
    ----------
    elements : list of str
    work_root : str
    poscar_dir : str
    pot_root_dir : str
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
                el, work_root=work_root, poscar_dir=poscar_dir,
                pot_root_dir=pot_root_dir, **kwargs)
            calc.run_all(submit=False)
            results[el] = 'generated'
        except Exception as e:
            results[el] = str(e)
    return results
