import unittest
import sys, os, tempfile, shutil, json
from unittest.mock import MagicMock, patch

# Setup shared mocks for heavy deps
from test_helpers import setup_mocks
setup_mocks()
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from extrempy.lazy.lib import ELEMENT_PHASE_DATA, get_phase_segments
from extrempy.lazy.potcar_map import PotcarMap, POTCAR_MAP


class TestSingleElementDpBuilderInit(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        # Create a minimal POTCAR dir so PotcarMap doesn't fail on build
        # The builder calls potcar_map.build during init? No, only on build_potcar.
        # But PotcarMap.__init__ doesn't call build, so we don't need fake POTCARs.

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    @patch('extrempy.campaign.single_element_dp.get_viable_elements')
    @patch('extrempy.campaign.single_element_dp.ElementDPBuilder')
    def test_build_all_elements_keyword(self, mock_builder_cls, mock_viable):
        mock_viable.return_value = ['Ti', 'Al', 'W']
        mock_builder = MagicMock()
        mock_builder_cls.return_value = mock_builder
        mock_builder.get_phase_segments.return_value = []
        mock_builder.generate_poscars.return_value = None
        mock_builder.generate_init_aimd.return_value = None
        mock_builder.submit_init_aimd.return_value = None

        from extrempy.campaign.single_element_dp import build_all_elements
        results = build_all_elements(work_root=self.tmpdir,
                                     potcar_lib=self.tmpdir,
                                     elements=['Ti'])
        self.assertIn('Ti', results)

    def test_init_basic(self):
        from extrempy.campaign.single_element_dp import ElementDPBuilder
        b = ElementDPBuilder('Ti', work_root=self.tmpdir,
                                   potcar_lib=self.tmpdir)
        self.assertEqual(b.element, 'Ti')
        self.assertEqual(b.elements, ['Ti'])
        self.assertEqual(b.work_dir, os.path.join(self.tmpdir, 'Ti'))

    def test_get_tm(self):
        from extrempy.campaign.single_element_dp import ElementDPBuilder
        b = ElementDPBuilder('Ti', work_root=self.tmpdir,
                                   potcar_lib=self.tmpdir)
        self.assertEqual(b._get_tm(), ELEMENT_PHASE_DATA['Ti']['Tm'])

    def test_get_tm_al(self):
        from extrempy.campaign.single_element_dp import ElementDPBuilder
        b = ElementDPBuilder('Al', work_root=self.tmpdir,
                                   potcar_lib=self.tmpdir)
        self.assertEqual(b._get_tm(), ELEMENT_PHASE_DATA['Al']['Tm'])

    def test_phase_segments_delegates(self):
        from extrempy.campaign.single_element_dp import ElementDPBuilder
        b = ElementDPBuilder('Ti', work_root=self.tmpdir,
                                   potcar_lib=self.tmpdir)
        segs = b.get_phase_segments()
        self.assertEqual(len(segs), 3)
        self.assertEqual(segs[0]['label'], 'Ti-HCP')

    def test_confs_dir_property(self):
        from extrempy.campaign.single_element_dp import ElementDPBuilder
        b = ElementDPBuilder('Al', work_root=self.tmpdir,
                                   potcar_lib=self.tmpdir)
        self.assertEqual(b.confs_dir,
                         os.path.join(b.work_dir, 'confs'))

    def test_init_vasp_dir_property(self):
        from extrempy.campaign.single_element_dp import ElementDPBuilder
        b = ElementDPBuilder('Al', work_root=self.tmpdir,
                                   potcar_lib=self.tmpdir)
        self.assertEqual(b.init_vasp_dir,
                         os.path.join(b.work_dir, 'init_vasp'))

    def test_init_data_dir_property(self):
        from extrempy.campaign.single_element_dp import ElementDPBuilder
        b = ElementDPBuilder('Al', work_root=self.tmpdir,
                                   potcar_lib=self.tmpdir)
        self.assertEqual(b.init_data_dir,
                         os.path.join(b.work_dir, 'init_data'))

    def test_dpgen_dir_property(self):
        from extrempy.campaign.single_element_dp import ElementDPBuilder
        b = ElementDPBuilder('Al', work_root=self.tmpdir,
                                   potcar_lib=self.tmpdir)
        self.assertEqual(b.dpgen_dir,
                         os.path.join(b.work_dir, 'dpgen'))


class TestExtremeDPBaseMethods(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        # Build needs fake POTCARs for ZVAL parsing
        from extrempy.lazy.potcar_map import PotcarMap
        self.create_fake_potcar('Al', '', ZVAL=3.0)
        self.create_fake_potcar('Ti', '_pv', ZVAL=10.0)

    def create_fake_potcar(self, element, variant, ZVAL=3.0):
        name = element + variant
        os.makedirs(os.path.join(self.tmpdir, name), exist_ok=True)
        path = os.path.join(self.tmpdir, name, 'POTCAR')
        with open(path, 'w') as f:
            f.write(f""" PAW_PBE {name} 06Sep2000
 TITEL = PAW_PBE {name} 06Sep2000
 POMASS = 1.00; ZVAL = {ZVAL}
""")

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_build_potcar(self):
        from extrempy.campaign.single_element_dp import DPBuilder
        b = DPBuilder(work_root=self.tmpdir,
                              potcar_lib=self.tmpdir)
        b.work_dir = self.tmpdir
        os.makedirs(os.path.join(self.tmpdir, 'dpgen'), exist_ok=True)
        zvals = b.build_potcar(['Al'])
        self.assertEqual(len(zvals), 1)
        self.assertAlmostEqual(zvals[0], 3.0)

    def test_build_potcar_writes_file(self):
        from extrempy.campaign.single_element_dp import DPBuilder
        b = DPBuilder(work_root=self.tmpdir,
                              potcar_lib=self.tmpdir)
        b.work_dir = self.tmpdir
        os.makedirs(os.path.join(self.tmpdir, 'dpgen'), exist_ok=True)
        b.build_potcar(['Al'])
        potcar_path = os.path.join(self.tmpdir, 'dpgen', 'POTCAR')
        self.assertTrue(os.path.exists(potcar_path))

    def test_ensure_dirs_creates_all(self):
        from extrempy.campaign.single_element_dp import DPBuilder
        b = DPBuilder(work_root=self.tmpdir,
                              potcar_lib=self.tmpdir)
        b.work_dir = self.tmpdir
        b._ensure_dirs()
        for d in [b.confs_dir, b.init_vasp_dir, b.init_data_dir,
                  b.dpgen_dir]:
            self.assertTrue(os.path.isdir(d),
                            f'{d} should exist')

    def test_abstract_methods(self):
        from extrempy.campaign.single_element_dp import DPBuilder
        b = DPBuilder(work_root=self.tmpdir,
                              potcar_lib=self.tmpdir)
        with self.assertRaises(NotImplementedError):
            b.get_phase_segments()
        with self.assertRaises(NotImplementedError):
            b.generate_poscars([])


class TestBuildAllElements(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        # Create fake POTCAR for Al
        os.makedirs(os.path.join(self.tmpdir, 'Al'), exist_ok=True)
        with open(os.path.join(self.tmpdir, 'Al', 'POTCAR'), 'w') as f:
            f.write(f""" PAW_PBE Al 06Sep2000
 TITEL = PAW_PBE Al 06Sep2000
 POMASS = 26.98; ZVAL = 3.0
""")

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    @patch('extrempy.campaign.single_element_dp.get_viable_elements')
    def test_build_all_smoke(self, mock_viable):
        mock_viable.return_value = ['Al']
        from extrempy.campaign.single_element_dp import build_all_elements
        results = build_all_elements(work_root=self.tmpdir,
                                     potcar_lib=self.tmpdir)
        # Al should be in results (status may vary)
        self.assertIn('Al', results)



class TestSysConfigsOrdering(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        os.makedirs(os.path.join(self.tmpdir, 'confs'), exist_ok=True)
        os.makedirs(os.path.join(self.tmpdir, 'dpgen'), exist_ok=True)
        # Create POSCAR files in non-alphabetical order to test ordering
        # Li has phases: α-Li(bcc), β-Li(bcc) → gets dedup labels Li-BCC, Li-BCC-2
        for fname in ['Li-BCC.POSCAR', 'Li-BCC-2.POSCAR', 'Li-LIQ.POSCAR']:
            path = os.path.join(self.tmpdir, 'confs', fname)
            with open(path, 'w') as f:
                f.write("test POSCAR\n")

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_sys_configs_order_matches_segs(self):
        from extrempy.campaign.single_element_dp import ElementDPBuilder
        from extrempy.lazy.dpgen import DPGENGenerator
        # Verify that sys_configs ordering logic matches segs order
        segs = [
            {'label': 'Li-BCC', 'structure': 'bcc', 'T_core': (0, 300), 'T_explore': (100, 600)},
            {'label': 'Li-BCC-2', 'structure': 'bcc', 'T_core': (300, 454), 'T_explore': (200, 600)},
            {'label': 'Li-LIQ', 'structure': 'bcc', 'T_core': (454, 600), 'T_explore': (350, 600)},
        ]
        b = ElementDPBuilder('Li', work_root=self.tmpdir,
                             potcar_lib=self.tmpdir)
        b.work_dir = self.tmpdir
        # Build sys_configs the same way as generate_dpgen does
        sys_configs = []
        missing = []
        for seg in segs:
            p = os.path.abspath(os.path.join(b.confs_dir, seg['label'] + '.POSCAR'))
            if os.path.exists(p):
                sys_configs.append([p])
            else:
                missing.append(seg['label'])
        self.assertEqual(len(missing), 0,
                         f"Missing POSCARs: {missing}")
        # Verify ordering: each entry corresponds to the seg at the same index
        for i, seg in enumerate(segs):
            expected_suffix = seg['label'] + '.POSCAR'
            self.assertTrue(sys_configs[i][0].endswith(expected_suffix),
                            f"sys_configs[{i}]={sys_configs[i][0]} does not match seg['{expected_suffix}']")
        # Also verify with Fe (BCC, FCC, BCC-2 order)
        segs_fe = [
            {'label': 'Fe-BCC', 'structure': 'bcc', 'T_core': (0, 1185), 'T_explore': (100, 1400)},
            {'label': 'Fe-FCC', 'structure': 'fcc', 'T_core': (1185, 1667), 'T_explore': (1000, 1900)},
            {'label': 'Fe-BCC-2', 'structure': 'bcc', 'T_core': (1667, 1811), 'T_explore': (1500, 2000)},
        ]
        os.makedirs(b.confs_dir, exist_ok=True)
        for fname in ['Fe-BCC.POSCAR', 'Fe-FCC.POSCAR', 'Fe-BCC-2.POSCAR']:
            path = os.path.join(b.confs_dir, fname)
            if not os.path.exists(path):
                with open(path, 'w') as f:
                    f.write("test POSCAR\n")
        sys_configs_fe = []
        missing_fe = []
        for seg in segs_fe:
            p = os.path.abspath(os.path.join(b.confs_dir, seg['label'] + '.POSCAR'))
            if os.path.exists(p):
                sys_configs_fe.append([p])
            else:
                missing_fe.append(seg['label'])
        self.assertEqual(len(missing_fe), 0)
        self.assertEqual(len(sys_configs_fe), 3)
        self.assertTrue(sys_configs_fe[0][0].endswith('Fe-BCC.POSCAR'))
        self.assertTrue(sys_configs_fe[1][0].endswith('Fe-FCC.POSCAR'))
        self.assertTrue(sys_configs_fe[2][0].endswith('Fe-BCC-2.POSCAR'))


# ---- Helper to build synthetic chunk.profile for EOS tests -------------

def _write_chunk_profile(path, verdict_type, n_blocks=20):
    """Write a synthetic chunk.profile producing the given verdict.

    verdict_type: 'solid' | 'liquid' | 'coexist'
    """
    if verdict_type == 'solid':
        q4_lo, q4_hi = 0.18, 0.19
        q6_lo, q6_hi = 0.50, 0.50
        rho_lo, rho_hi = 2.70, 2.70
    elif verdict_type == 'liquid':
        q4_lo, q4_hi = 0.02, 0.01
        q6_lo, q6_hi = 0.01, 0.01
        rho_lo, rho_hi = 2.50, 2.50
    else:  # coexist
        q4_lo, q4_hi = 0.18, 0.02
        q6_lo, q6_hi = 0.50, 0.01
        rho_lo, rho_hi = 2.70, 2.50

    header = '# Chunk Coord1 Ncount density/mass temp v_virial_atom c_Q[1] c_Q[2]'
    with open(path, 'w') as f:
        for blk in range(n_blocks):
            ts = 1000 * (blk + 1)
            f.write(f'{ts}\n')
            f.write(header + '\n')
            for i in range(10):
                half = i // 5
                q4 = q4_lo if half == 0 else q4_hi
                q6 = q6_lo if half == 0 else q6_hi
                rho = rho_lo if half == 0 else rho_hi
                coord = 0.05 + 0.1 * i
                f.write(f'{i+1} {coord:.4f} 100 {rho:.4f} 300.0 0.0 {q4:.4f} {q6:.4f}\n')


_REAL_NUMPY = None  # cached real numpy module (loaded once; C ext can't reload)


def _restore_real_numpy():
    """Replace the mocked numpy with the real one and reload modules
    that captured the mock at import time (melt.py, chunk.py).

    Returns a dict to pass to ``_restore_mocked_numpy`` in tearDown.
    The real numpy is cached at module level because numpy's C extension
    cannot be reloaded more than once per process.

    NB: the submodule-injection loop below is conservative — it only
    registers submodules that already exist as attributes on the cached
    numpy module object.  This covers melt.py's current usage
    (``np.arange / np.median / np.mean``).  If a future change needs
    ``np.linalg.*`` or other lazily-imported submodules, add an explicit
    ``import numpy.linalg`` here to ensure it is registered.
    """
    global _REAL_NUMPY
    import importlib
    saved = {}
    for mod in list(sys.modules):
        if mod == 'numpy' or mod.startswith('numpy.'):
            saved[mod] = sys.modules.pop(mod)
    if _REAL_NUMPY is None:
        _REAL_NUMPY = importlib.import_module('numpy')
    # Inject the cached real numpy (and its submodules) into sys.modules.
    sys.modules['numpy'] = _REAL_NUMPY
    for name in dir(_REAL_NUMPY):
        sub = getattr(_REAL_NUMPY, name, None)
        if isinstance(sub, type(_REAL_NUMPY)):
            sys.modules[f'numpy.{name}'] = sub
    # Force reload of melt (captured mock np) and chunk (same).
    for mod_name in ('extrempy.campaign.chunk', 'extrempy.campaign.melt'):
        if mod_name in sys.modules:
            importlib.reload(sys.modules[mod_name])
    return saved


def _restore_mocked_numpy(saved):
    sys.modules.update(saved)


class TestEOSAnalyzer(unittest.TestCase):
    """Test analyze_two_phase with synthetic chunk.profile + Tm_refined.

    Note: analyze_two_phase / _get_npt_temps use numpy for real, so we
    restore the real numpy module (test_helpers.setup_mocks mocks it).
    """

    def setUp(self):
        self._np_saved = _restore_real_numpy()
        self.tmpdir = tempfile.mkdtemp()
        from extrempy.campaign.melt import ElementEOSCalculator
        # Al Tm=933; create 3 candidate-T dirs: 1133k(solid),1233k(coexist),1333k(liquid)
        self.temps = [1133, 1233, 1333]
        self.verdicts = ['solid', 'coexist', 'liquid']
        for T, v in zip(self.temps, self.verdicts):
            d = os.path.join(self.tmpdir, 'Al', 'melt', f'{T}k')
            os.makedirs(d)
            _write_chunk_profile(os.path.join(d, 'chunk.profile'), v)
        # Build calc with fake gens
        self.calc = ElementEOSCalculator('Al', work_root=self.tmpdir, dpgen_dir=None)
        self.calc._two_phase_gens = []
        for T in self.temps:
            gen = MagicMock()
            gen.params = {'Tm_estimate': T}
            gen.work_path = os.path.join(self.tmpdir, 'Al', 'melt', f'{T}k')
            self.calc._two_phase_gens.append(gen)

    def tearDown(self):
        shutil.rmtree(self.tmpdir)
        _restore_mocked_numpy(self._np_saved)

    def test_coexist_detection(self):
        result = self.calc.analyze_two_phase()
        self.assertEqual(result['coexist_temps'], [1233])
        self.assertEqual(result['Tm_interval'], (1233, 1233))
        self.assertEqual(result['Tm_refined'], 1233)
        self.assertEqual(self.calc.Tm_refined, 1233)

    def test_no_coexist(self):
        # Make all temps solid → no coexist
        for T in self.temps:
            d = os.path.join(self.tmpdir, 'Al', 'melt', f'{T}k')
            _write_chunk_profile(os.path.join(d, 'chunk.profile'), 'solid')
        result = self.calc.analyze_two_phase()
        self.assertEqual(result['coexist_temps'], [])
        self.assertIsNone(result['Tm_interval'])
        self.assertIsNone(result['Tm_refined'])

    def test_submit_without_generate(self):
        # BUG-5: submit_two_phase without generate should not raise AttributeError
        from extrempy.campaign.melt import ElementEOSCalculator
        calc = ElementEOSCalculator('Al', work_root=self.tmpdir, dpgen_dir=None)
        calc.submit_two_phase(submit=False)  # no-op, no crash
        calc.submit_npt(submit=False)
        calc.submit_nvt_traj(submit=False)

    def test_tm_refined_used_by_npt(self):
        # Set Tm_refined; _get_npt_temps should center on it
        from extrempy.campaign.melt import ElementEOSCalculator
        calc = ElementEOSCalculator('Al', work_root=self.tmpdir, dpgen_dir=None)
        calc.Tm_refined = 1100
        temps = calc._get_npt_temps(npt_n=5, npt_dT=100, npt_shift=0)
        # 1100 + [-200,-100,0,100,200] = [900,1000,1100,1200,1300]
        self.assertEqual(temps, [900, 1000, 1100, 1200, 1300])


class TestSbatchGeneration(unittest.TestCase):
    """Test _write_sbatch: cpus-per-task, needs_traj_dir."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def _make_calc(self):
        from extrempy.campaign.melt import ElementEOSCalculator
        return ElementEOSCalculator('Al', work_root=self.tmpdir, dpgen_dir=None)

    def test_cpus_per_task_appended(self):
        calc = self._make_calc()
        cfg = calc._resolve_slurm_config()
        cfg['cores_per_node'] = 40
        cfg['ntasks_per_node'] = 20
        path = calc._write_sbatch(cfg, 'testjob', self.tmpdir, needs_traj_dir=False)
        with open(path) as f:
            content = f.read()
        self.assertIn('#SBATCH --cpus-per-task=2', content)
        self.assertIn('export OMP_NUM_THREADS=2', content)

    def test_no_mkdir_traj_for_npt(self):
        calc = self._make_calc()
        cfg = calc._resolve_slurm_config()
        path = calc._write_sbatch(cfg, 'Al_933k_npt_solid', self.tmpdir,
                                  needs_traj_dir=False)
        with open(path) as f:
            content = f.read()
        self.assertNotIn('mkdir -p traj', content)

    def test_mkdir_traj_for_nvt(self):
        calc = self._make_calc()
        cfg = calc._resolve_slurm_config()
        path = calc._write_sbatch(cfg, 'Al_933k_nvt_solid', self.tmpdir,
                                  needs_traj_dir=True)
        with open(path) as f:
            content = f.read()
        self.assertIn('mkdir -p traj', content)

    def test_machine_template_command_used(self):
        # D-6: machine_template command should now be respected
        from extrempy.campaign.melt import ElementEOSCalculator
        # Create a minimal machine.json
        mj_path = os.path.join(self.tmpdir, 'machine.json')
        with open(mj_path, 'w') as f:
            json.dump({
                'model_devi': {
                    'command': '/usr/local/bin/lmp_custom -in run.in',
                    'resources': {'number_node': 2, 'cpu_per_node': 32}
                }
            }, f)
        calc = ElementEOSCalculator('Al', work_root=self.tmpdir, dpgen_dir=None,
                                    machine_template=mj_path)
        cfg = calc._resolve_slurm_config()
        self.assertEqual(cfg['command'], '/usr/local/bin/lmp_custom -in run.in')
        self.assertEqual(cfg['nodes'], 2)


class TestNptShiftDefault(unittest.TestCase):
    def test_npt_shift_zero(self):
        saved = _restore_real_numpy()
        try:
            from extrempy.campaign.melt import ElementEOSCalculator
            tmpdir = tempfile.mkdtemp()
            try:
                calc = ElementEOSCalculator('Al', work_root=tmpdir, dpgen_dir=None)
                # Al Tm=933: [733,833,933,1033,1133]
                temps = calc._get_npt_temps()
                self.assertEqual(temps, [733, 833, 933, 1033, 1133])
            finally:
                shutil.rmtree(tmpdir)
        finally:
            _restore_mocked_numpy(saved)


class TestRunAllSplit(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def _make_calc_with_mocks(self):
        from extrempy.campaign.melt import ElementEOSCalculator
        calc = ElementEOSCalculator('Al', work_root=self.tmpdir, dpgen_dir=None)
        calc.generate_two_phase = MagicMock()
        calc.submit_two_phase = MagicMock()
        calc.generate_npt = MagicMock()
        calc.submit_npt = MagicMock()
        calc.generate_nvt_traj = MagicMock()
        calc.submit_nvt_traj = MagicMock()
        return calc

    def test_run_two_phase_only(self):
        calc = self._make_calc_with_mocks()
        calc.run_two_phase(submit=False)
        calc.generate_two_phase.assert_called_once()
        calc.submit_two_phase.assert_called_once_with(submit=False)
        calc.generate_npt.assert_not_called()

    def test_run_property_scans_only(self):
        calc = self._make_calc_with_mocks()
        calc.run_property_scans(submit=False)
        calc.generate_npt.assert_called_once()
        calc.submit_npt.assert_called_once_with(submit=False)
        calc.generate_nvt_traj.assert_called_once()
        calc.generate_two_phase.assert_not_called()

    def test_run_all_equals_both(self):
        calc = self._make_calc_with_mocks()
        calc.run_all(submit=False)
        calc.generate_two_phase.assert_called_once()
        calc.generate_npt.assert_called_once()
        calc.generate_nvt_traj.assert_called_once()


class TestEOSAbstract(unittest.TestCase):
    def test_base_class_not_instantiable(self):
        from extrempy.campaign.melt import EOSCalculator
        with self.assertRaises(TypeError):
            EOSCalculator(work_root='/tmp')

    def test_subclass_works(self):
        from extrempy.campaign.melt import ElementEOSCalculator
        tmpdir = tempfile.mkdtemp()
        try:
            calc = ElementEOSCalculator('Al', work_root=tmpdir, dpgen_dir=None)
            self.assertEqual(calc.element, 'Al')
        finally:
            shutil.rmtree(tmpdir)


if __name__ == '__main__':
    unittest.main()
