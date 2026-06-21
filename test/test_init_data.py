import unittest
import sys, os, tempfile, shutil
from unittest.mock import MagicMock, patch, call

# Setup shared mocks for heavy deps
from test_helpers import setup_mocks
setup_mocks()
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))


class TestGenerateLiquidPoscar(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def _make_contcar(self, label, content='POSCAR content\n'):
        subdir = os.path.join(self.tmpdir, label)
        os.makedirs(subdir, exist_ok=True)
        path = os.path.join(subdir, 'CONTCAR')
        with open(path, 'w') as f:
            f.write(content)
        return subdir

    def test_liquid_copy_contcar(self):
        from extrempy.lazy.init_data import generate_liquid_poscar_from_contcar
        aimd_dir = self._make_contcar('Ti-LIQ', 'test data')
        result = generate_liquid_poscar_from_contcar(
            ('Ti-LIQ', aimd_dir), self.tmpdir)
        expected = os.path.join(self.tmpdir, 'Ti-LIQ.POSCAR')
        self.assertEqual(result, expected)
        self.assertTrue(os.path.exists(expected))
        with open(expected) as f:
            self.assertEqual(f.read(), 'test data')

    def test_non_liquid_returns_none(self):
        from extrempy.lazy.init_data import generate_liquid_poscar_from_contcar
        aimd_dir = self._make_contcar('Ti-HCP', 'data')
        result = generate_liquid_poscar_from_contcar(
            ('Ti-HCP', aimd_dir), self.tmpdir)
        self.assertIsNone(result)

    def test_missing_contcar(self):
        from extrempy.lazy.init_data import generate_liquid_poscar_from_contcar
        result = generate_liquid_poscar_from_contcar(
            ('Ti-LIQ', '/nonexistent'), self.tmpdir)
        self.assertIsNone(result)

    def test_overwrite_existing(self):
        from extrempy.lazy.init_data import generate_liquid_poscar_from_contcar
        # Create an existing LIQ POSCAR
        old_path = os.path.join(self.tmpdir, 'Ti-LIQ.POSCAR')
        with open(old_path, 'w') as f:
            f.write('old')
        # New CONTCAR should overwrite
        aimd_dir = self._make_contcar('Ti-LIQ', 'new data')
        generate_liquid_poscar_from_contcar(('Ti-LIQ', aimd_dir), self.tmpdir)
        with open(old_path) as f:
            self.assertEqual(f.read(), 'new data')


class TestScalePoscarVolume(unittest.TestCase):
    def _make_system_mock(self, atom_names, atom_numbs, n_atoms):
        from unittest.mock import MagicMock
        class Cell3x3(list):
            def __mul__(self, other):
                return Cell3x3([row[:] for row in self])
        sys_mock = MagicMock()
        cell = Cell3x3([[4.0, 0.0, 0.0],
                        [0.0, 4.0, 0.0],
                        [0.0, 0.0, 4.0]])
        coords = Cell3x3([[1.0] * 3 for _ in range(n_atoms)])
        sys_mock.__getitem__.side_effect = lambda k: {
            'coords': [[list(c) for c in coords]],
            'cells': [cell],
            'atom_names': atom_names,
            'atom_numbs': atom_numbs,
        }[k]
        return sys_mock

    def test_scales_lattice(self):
        from unittest.mock import patch, MagicMock
        with patch('extrempy.lazy.init_data.dpdata.System') as mock_dpsys:
            mock_dpsys.return_value = self._make_system_mock(
                ['Ti'], [54], 54)
            from extrempy.lazy.init_data import scale_poscar_volume
            scale_poscar_volume('/fake/POSCAR', 1.10, '/fake/out')
            from ase import Atoms
            Atoms.assert_called_once()
            call_kw = Atoms.call_args[1]
            self.assertEqual(call_kw['symbols'], ['Ti'] * 54)
            self.assertTrue(call_kw['pbc'])

    def test_writes_to_out_path(self):
        from unittest.mock import patch, MagicMock
        from ase import Atoms
        Atoms.reset_mock()
        with patch('extrempy.lazy.init_data.dpdata.System') as mock_dpsys:
            mock_dpsys.return_value = self._make_system_mock(
                ['Al'], [4], 4)
            from extrempy.lazy.init_data import scale_poscar_volume
            out = scale_poscar_volume('/fake/POSCAR', 1.0,
                                      '/tmp/fake_out_path')
            Atoms.return_value.write.assert_any_call(
                '/tmp/fake_out_path', format='vasp', direct=True)


class TestBootstrapInitDataSmoke(unittest.TestCase):
    def test_skip_missing_outcar(self):
        from unittest.mock import patch
        with patch('extrempy.lazy.init_data.dpdata.LabeledSystem') as mock_ls:
            from extrempy.lazy.init_data import bootstrap_init_data
            result = bootstrap_init_data(
                [('Ti-HCP', '/nonexistent')],
                '/tmp/fake_root', 'raw_to_set.sh')
            self.assertEqual(result, [])
            mock_ls.assert_not_called()

    def test_skip_too_few_frames(self):
        from unittest.mock import patch, MagicMock
        with patch('extrempy.lazy.init_data.dpdata.LabeledSystem') as mock_ls:
            mock_sys = MagicMock()
            mock_sys.get_nframes.return_value = 100
            mock_ls.return_value = mock_sys
            import tempfile
            with tempfile.TemporaryDirectory() as tmp:
                from extrempy.lazy.init_data import bootstrap_init_data
                incar_dir = os.path.join(tmp, 'Ti-HCP')
                os.makedirs(incar_dir, exist_ok=True)
                with open(os.path.join(incar_dir, 'INCAR'), 'w') as f:
                    f.write('TEBEG = 800\n')
                result = bootstrap_init_data(
                    [('Ti-HCP', incar_dir)],
                    os.path.join(tmp, 'init_data'),
                    drop_first=200)

    def test_liquid_uses_liq_stride(self):
        from unittest.mock import patch, MagicMock
        with patch('extrempy.lazy.init_data.dpdata.LabeledSystem') as mock_ls, \
             patch('extrempy.lazy.init_data.np.savetxt'):
            mock_sys = MagicMock()
            mock_sys.get_nframes.return_value = 500
            mock_sys.__getitem__.return_value = mock_sys
            mock_ls.return_value = mock_sys
            import tempfile
            with tempfile.TemporaryDirectory() as tmp:
                aimd_dir = os.path.join(tmp, 'Ti-LIQ')
                os.makedirs(aimd_dir, exist_ok=True)
                with open(os.path.join(aimd_dir, 'OUTCAR'), 'w') as f:
                    f.write("dummy outcar content")
                with open(os.path.join(aimd_dir, 'INCAR'), 'w') as f:
                    f.write('TEBEG = 3494\n')
                from extrempy.lazy.init_data import bootstrap_init_data
                init_sys = os.path.join(tmp, 'init_data')
                result = bootstrap_init_data(
                    [('Ti-LIQ', aimd_dir)], init_sys,
                    drop_first=200, solid_stride=50, liq_stride=30)
                self.assertIn('Ti-LIQ', result)


if __name__ == '__main__':
    unittest.main()
