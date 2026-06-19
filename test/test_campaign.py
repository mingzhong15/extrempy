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
    @patch('extrempy.campaign.single_element_dp.SingleElementDPBuilder')
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
        from extrempy.campaign.single_element_dp import SingleElementDPBuilder
        b = SingleElementDPBuilder('Ti', work_root=self.tmpdir,
                                   potcar_lib=self.tmpdir)
        self.assertEqual(b.element, 'Ti')
        self.assertEqual(b.elements, ['Ti'])
        self.assertEqual(b.work_dir, os.path.join(self.tmpdir, 'Ti'))

    def test_get_tm(self):
        from extrempy.campaign.single_element_dp import SingleElementDPBuilder
        b = SingleElementDPBuilder('Ti', work_root=self.tmpdir,
                                   potcar_lib=self.tmpdir)
        self.assertEqual(b._get_tm(), ELEMENT_PHASE_DATA['Ti']['Tm'])

    def test_get_tm_al(self):
        from extrempy.campaign.single_element_dp import SingleElementDPBuilder
        b = SingleElementDPBuilder('Al', work_root=self.tmpdir,
                                   potcar_lib=self.tmpdir)
        self.assertEqual(b._get_tm(), ELEMENT_PHASE_DATA['Al']['Tm'])

    def test_phase_segments_delegates(self):
        from extrempy.campaign.single_element_dp import SingleElementDPBuilder
        b = SingleElementDPBuilder('Ti', work_root=self.tmpdir,
                                   potcar_lib=self.tmpdir)
        segs = b.get_phase_segments()
        self.assertEqual(len(segs), 3)
        self.assertEqual(segs[0]['label'], 'Ti-HCP')

    def test_confs_dir_property(self):
        from extrempy.campaign.single_element_dp import SingleElementDPBuilder
        b = SingleElementDPBuilder('Al', work_root=self.tmpdir,
                                   potcar_lib=self.tmpdir)
        self.assertEqual(b.confs_dir,
                         os.path.join(b.work_dir, 'confs'))

    def test_init_vasp_dir_property(self):
        from extrempy.campaign.single_element_dp import SingleElementDPBuilder
        b = SingleElementDPBuilder('Al', work_root=self.tmpdir,
                                   potcar_lib=self.tmpdir)
        self.assertEqual(b.init_vasp_dir,
                         os.path.join(b.work_dir, 'init_vasp'))

    def test_init_data_dir_property(self):
        from extrempy.campaign.single_element_dp import SingleElementDPBuilder
        b = SingleElementDPBuilder('Al', work_root=self.tmpdir,
                                   potcar_lib=self.tmpdir)
        self.assertEqual(b.init_data_dir,
                         os.path.join(b.work_dir, 'init_data'))

    def test_dpgen_dir_property(self):
        from extrempy.campaign.single_element_dp import SingleElementDPBuilder
        b = SingleElementDPBuilder('Al', work_root=self.tmpdir,
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
        from extrempy.campaign.single_element_dp import ExtremeDPBuilder
        b = ExtremeDPBuilder(work_root=self.tmpdir,
                              potcar_lib=self.tmpdir)
        b.work_dir = self.tmpdir
        os.makedirs(os.path.join(self.tmpdir, 'dpgen'), exist_ok=True)
        zvals = b.build_potcar(['Al'])
        self.assertEqual(len(zvals), 1)
        self.assertAlmostEqual(zvals[0], 3.0)

    def test_build_potcar_writes_file(self):
        from extrempy.campaign.single_element_dp import ExtremeDPBuilder
        b = ExtremeDPBuilder(work_root=self.tmpdir,
                              potcar_lib=self.tmpdir)
        b.work_dir = self.tmpdir
        os.makedirs(os.path.join(self.tmpdir, 'dpgen'), exist_ok=True)
        b.build_potcar(['Al'])
        potcar_path = os.path.join(self.tmpdir, 'dpgen', 'POTCAR')
        self.assertTrue(os.path.exists(potcar_path))

    def test_ensure_dirs_creates_all(self):
        from extrempy.campaign.single_element_dp import ExtremeDPBuilder
        b = ExtremeDPBuilder(work_root=self.tmpdir,
                              potcar_lib=self.tmpdir)
        b.work_dir = self.tmpdir
        b._ensure_dirs()
        for d in [b.confs_dir, b.init_vasp_dir, b.init_data_dir,
                  b.dpgen_dir]:
            self.assertTrue(os.path.isdir(d),
                            f'{d} should exist')

    def test_abstract_methods(self):
        from extrempy.campaign.single_element_dp import ExtremeDPBuilder
        b = ExtremeDPBuilder(work_root=self.tmpdir,
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


if __name__ == '__main__':
    unittest.main()
