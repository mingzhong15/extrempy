import unittest
import sys, os
from unittest.mock import MagicMock

# Setup shared mocks for heavy deps
from test_helpers import setup_mocks
setup_mocks()
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from extrempy.lazy.vasp import _incar_dict, _render_incar
from extrempy.constant import kb_eV


class TestIncarDict(unittest.TestCase):
    def setUp(self):
        self.common = dict(encut=600, nbands=200, ele_temp_K=300,
                           latt_temp_K=500, md_steps=500, dt=1)

    def test_scf_preset(self):
        d = _incar_dict('scf', **self.common)
        self.assertEqual(d['IBRION'], -1)
        self.assertEqual(d['NSW'], 0)
        self.assertNotIn('ISIF', d)
        self.assertNotIn('SMASS', d)
        self.assertNotIn('MDALGO', d)

    def test_relax_preset(self):
        d = _incar_dict('relax', **self.common)
        self.assertEqual(d['IBRION'], 1)
        self.assertEqual(d['NSW'], self.common['md_steps'])
        self.assertEqual(d['ISIF'], 2)

    def test_aimd_ttm_preset(self):
        d = _incar_dict('aimd-ttm', **self.common)
        self.assertEqual(d['IBRION'], 0)
        self.assertEqual(d['ISIF'], 2)
        self.assertEqual(d['NSW'], self.common['md_steps'])
        self.assertEqual(d['SMASS'], 0)
        self.assertEqual(d['MDALGO'], 2)
        self.assertEqual(d['TEBEG'], 500)
        self.assertEqual(d['TEEND'], 500)

    def test_aimd_ttm_latt_temp_none(self):
        # latt_temp is passed, not None, so always explicit
        d = _incar_dict('aimd-ttm', **self.common)
        self.assertEqual(d['TEBEG'], 500)

    def test_sigma_uses_kb_eV(self):
        T_high = 11600
        d = _incar_dict('scf', encut=600, nbands=100, ele_temp_K=T_high,
                        latt_temp_K=T_high, md_steps=0, dt=1)
        expected_sigma = T_high * kb_eV
        parsed = float(d['SIGMA'])
        self.assertAlmostEqual(parsed, expected_sigma, delta=1e-6)

    def test_unknown_preset(self):
        with self.assertRaises(ValueError):
            _incar_dict('bad_preset', **self.common)

    def test_magnetic_adds_ispin(self):
        d = _incar_dict('scf', encut=600, nbands=100, ele_temp_K=300,
                        md_steps=0, dt=1, is_magnetic=True)
        self.assertEqual(d['ISPIN'], 2)

    def test_magnetic_false(self):
        d = _incar_dict('scf', encut=600, nbands=100, ele_temp_K=300,
                        md_steps=0, dt=1, is_magnetic=False)
        self.assertNotIn('ISPIN', d)

    def test_encut_passed(self):
        d = _incar_dict('scf', **self.common)
        self.assertEqual(d['ENCUT'], 600)

    def test_nbands_passed(self):
        d = _incar_dict('scf', **self.common)
        self.assertEqual(d['NBANDS'], 200)

    def test_prec_default(self):
        d = _incar_dict('scf', **self.common)
        self.assertEqual(d['PREC'], 'High')

    def test_kspacing_default(self):
        d = _incar_dict('scf', **self.common)
        self.assertEqual(d['KSPACING'], 0.5)

    def test_lwave_false(self):
        d = _incar_dict('scf', **self.common)
        self.assertEqual(d['LWAVE'], '.FALSE.')

    def test_lcharg_false(self):
        d = _incar_dict('scf', **self.common)
        self.assertEqual(d['LCHARG'], '.FALSE.')


class TestRenderIncar(unittest.TestCase):
    def setUp(self):
        self.d = _incar_dict('scf', encut=600, nbands=200,
                             ele_temp_K=300, md_steps=0, dt=1)

    def test_contains_sections(self):
        text = _render_incar(self.d)
        for section in ['CONTROL', 'ELECTRON', 'XC FUNCTIONAL',
                        'ION', 'PARALLIZATION', 'K-POINTS']:
            self.assertIn(f'# {section}', text)

    def test_key_value_format(self):
        text = _render_incar(self.d)
        lines = [l for l in text.split('\n') if '=' in l]
        for line in lines:
            self.assertIn(' = ', line,
                          f'Line missing " = " format: {line}')
            parts = line.split(' = ')
            self.assertEqual(len(parts), 2,
                             f'Line has unexpected format: {line}')

    def test_does_not_mutate_input(self):
        original = dict(self.d)
        _render_incar(self.d)
        self.assertEqual(self.d, original)

    def test_trailing_newline(self):
        text = _render_incar(self.d)
        self.assertTrue(text.endswith('\n'))

    def test_iblank_line_separator(self):
        text = _render_incar(self.d)
        # Sections separated by blank lines
        self.assertIn('\n\n', text)

    def test_override_encut(self):
        d2 = _incar_dict('scf', encut=600, nbands=100,
                         ele_temp_K=300, md_steps=0, dt=1)
        text = _render_incar(d2)
        self.assertIn('ENCUT = 600', text)

    def test_aimd_ttm_has_tebeg(self):
        d2 = _incar_dict('aimd-ttm', encut=600, nbands=100,
                         ele_temp_K=300, latt_temp_K=800,
                         md_steps=500, dt=1)
        text = _render_incar(d2)
        self.assertIn('TEBEG = 800', text)
        self.assertIn('TEEND = 800', text)

    def test_extra_keys_go_at_bottom(self):
        d2 = dict(self.d, EXTRA='VALUE')
        text = _render_incar(d2)
        lines = text.strip().split('\n')
        self.assertEqual(lines[-1], 'EXTRA = VALUE')


if __name__ == '__main__':
    unittest.main()
