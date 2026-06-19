import unittest
import sys, os, tempfile, shutil
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from extrempy.lazy.potcar_map import PotcarMap, POTCAR_MAP, FALLBACK_VARIANTS


def _make_fake_potcar(dirpath, element, variant='', ZVAL=10.0, mass=47.867):
    """Create a minimal VASP-format POTCAR file for testing."""
    name = element + variant
    os.makedirs(os.path.join(dirpath, name), exist_ok=True)
    path = os.path.join(dirpath, name, 'POTCAR')
    content = f""" PAW_PBE {element}{variant} 06Sep2000
 TITEL = PAW_PBE {element}{variant} 06Sep2000
 the quick brown fox
 POMASS = {mass}; ZVAL = {ZVAL}
 some other data
 """
    with open(path, 'w') as f:
        f.write(content)
    return path


class TestPotcarMapInit(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_init_valid_set(self):
        pm = PotcarMap('PBE54', self.tmpdir)
        self.assertEqual(pm.set_name, 'PBE54')
        self.assertEqual(pm.potcar_dir, self.tmpdir)

    def test_init_invalid_set(self):
        with self.assertRaises(KeyError):
            PotcarMap('NONEXIST', self.tmpdir)

    def test_map_contains_key_elements(self):
        pm = PotcarMap('PBE54', self.tmpdir)
        for el in ['Ti', 'Al', 'Cu', 'W', 'Li']:
            self.assertIn(el, pm.map, f'{el} should be in PBE54 map')

    def test_pbe52_is_empty(self):
        self.assertEqual(POTCAR_MAP['PBE52'], {})
        with self.assertRaises(KeyError):
            PotcarMap('PBE52', self.tmpdir).build(['Ti'])


class TestPotcarMapResolve(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        _make_fake_potcar(self.tmpdir, 'Al', '', ZVAL=3.0, mass=26.98)
        _make_fake_potcar(self.tmpdir, 'Ti', '_pv', ZVAL=10.0, mass=47.87)
        _make_fake_potcar(self.tmpdir, 'Ti', '', ZVAL=4.0, mass=47.87)
        self.pm = PotcarMap('PBE54', self.tmpdir)

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_resolve_exact_match(self):
        path, variant = self.pm._resolve_path('Al')
        self.assertTrue(os.path.exists(path))
        self.assertEqual(variant, '')
        self.assertIn('Al', path)

    def test_resolve_fallback(self):
        # Ti prefers '_pv', should find it even before trying ''
        path, variant = self.pm._resolve_path('Ti')
        self.assertEqual(variant, '_pv')

    def test_resolve_missing_element(self):
        with self.assertRaises(KeyError):
            self.pm._resolve_path('NoSuchElement')

    def test_resolve_no_potcar_file(self):
        pm2 = PotcarMap('PBE54', tempfile.mkdtemp())
        with self.assertRaises(FileNotFoundError):
            pm2._resolve_path('Al')

    def test_resolve_custom_variant(self):
        path, variant = self.pm._resolve_path('Ti')
        self.assertIn('Ti_pv', path)

    def test_resolve_fallback_order(self):
        # Only Ti_pv exists (no Ti), should find Ti_pv
        path, variant = self.pm._resolve_path('Ti')
        self.assertEqual(variant, '_pv')


class TestPotcarMapParsing(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.path = _make_fake_potcar(self.tmpdir, 'Ti', '_pv',
                                       ZVAL=10.0, mass=47.87)
        self.pm = PotcarMap('PBE54', self.tmpdir)

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_parse_zval(self):
        zval = self.pm._parse_zval(self.path)
        self.assertAlmostEqual(zval, 10.0)

    def test_parse_zval_integer_returned_as_float(self):
        zval = self.pm._parse_zval(self.path)
        self.assertIsInstance(zval, float)

    def test_parse_titel(self):
        titel = self.pm._parse_titel(self.path)
        self.assertIn('TITEL', titel)
        self.assertIn('Ti_pv', titel)


class TestPotcarMapBuild(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        _make_fake_potcar(self.tmpdir, 'Al', '', ZVAL=3.0, mass=26.98)
        _make_fake_potcar(self.tmpdir, 'Ti', '_pv', ZVAL=10.0, mass=47.87)
        _make_fake_potcar(self.tmpdir, 'Cu', '_pv', ZVAL=11.0, mass=63.55)
        self.pm = PotcarMap('PBE54', self.tmpdir)

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_build_single_element(self):
        text, zvals, paths = self.pm.build(['Ti'])
        self.assertEqual(len(zvals), 1)
        self.assertAlmostEqual(zvals[0], 10.0)

    def test_build_multi_element(self):
        text, zvals, paths = self.pm.build(['Al', 'Ti', 'Cu'])
        self.assertEqual(len(zvals), 3)
        self.assertAlmostEqual(zvals[0], 3.0)
        self.assertAlmostEqual(zvals[1], 10.0)
        self.assertAlmostEqual(zvals[2], 11.0)

    def test_build_zval_fill(self):
        # ZVAL=0 in dict should be auto-filled
        self.assertEqual(self.pm.map['Al']['ZVAL'], 0)
        self.pm.build(['Al'])
        self.assertAlmostEqual(self.pm.map['Al']['ZVAL'], 3.0)

    def test_build_zval_mismatch(self):
        # Force a wrong ZVAL in map
        self.pm.map['Al']['ZVAL'] = 99.0
        with self.assertRaises(ValueError):
            self.pm.build(['Al'])

    def test_build_missing_element_in_map(self):
        with self.assertRaises(KeyError):
            self.pm.build(['NoSuchElement'])

    def test_build_no_file_found(self):
        # Element in map but no POTCAR file
        with self.assertRaises(FileNotFoundError):
            self.pm.build(['Ag'])  # no Ag POTCAR in tmpdir

    def test_build_concatenated_output(self):
        text, zvals, paths = self.pm.build(['Al', 'Ti'])
        content = text.decode()
        self.assertIn('PAW_PBE Al', content)
        self.assertIn('PAW_PBE Ti_pv', content)
        # Al's content should come before Ti's
        self.assertLess(content.index('PAW_PBE Al'),
                        content.index('PAW_PBE Ti_pv'))

    def test_write_potcar(self):
        outpath = os.path.join(self.tmpdir, 'combined_POTCAR')
        zvals, paths = self.pm.write_potcar(['Al', 'Ti'], outpath)
        self.assertTrue(os.path.exists(outpath))
        with open(outpath, 'rb') as f:
            content = f.read()
        self.assertIn(b'PAW_PBE Al', content)
        self.assertIn(b'PAW_PBE Ti_pv', content)


if __name__ == '__main__':
    unittest.main()
