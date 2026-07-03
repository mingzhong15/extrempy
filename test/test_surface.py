"""Tests for extrempy.lazy.surface (make_slab / make_slabs).

Mock-based tests stub heavy deps (numpy, ase, scipy) via
test_helpers.setup_mocks() so they run without those packages
installed.  Integration tests at the bottom swap in the real ASE
modules and are skipped when ASE is unavailable.
"""
import sys, os, tempfile, shutil, importlib, importlib.util, unittest
from unittest.mock import MagicMock, patch

# Check for real ASE/numpy/scipy availability.  We must probe BEFORE
# test_helpers.setup_mocks() replaces sys.modules with mock objects,
# but we also need to be robust to other test files (e.g. test_campaign)
# having run before us and left sys.modules['ase'].__spec__ = None.
def _spec_available(name):
    """Return True iff *name* is importable on this system.

    Bypasses sys.modules entirely so mock-replaced packages don't
    confuse the probe.
    """
    try:
        # Use PathFinder to look up the spec directly without consulting
        # sys.modules (which may have mock objects with __spec__ = None).
        from importlib.machinery import PathFinder
        spec = PathFinder.find_spec(name)
        return spec is not None
    except (ImportError, ValueError, AttributeError):
        return False

HAS_ASE = (_spec_available('ase')
           and _spec_available('numpy')
           and _spec_available('scipy'))

# Stub heavy deps before importing extrempy.
from test_helpers import setup_mocks
setup_mocks()

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import extrempy.lazy.surface as surf_mod
from extrempy.lazy.surface import make_slab, make_slabs


class _FakeAtoms:
    """Lightweight stand-in for ase.Atoms (has get_positions, __len__)."""
    def __init__(self, n=4):
        self._n = n
    def __len__(self):
        return self._n
    def get_positions(self):
        return [(0, 0, 0)] * self._n


class TestMakeSlabLogic(unittest.TestCase):
    """Mock-based tests for make_slab argument validation and the
    _to_surface_supercell_matrix helper (no real ASE)."""

    # ── argument validation ──
    def test_layers_and_z_min_both_none_raises(self):
        """Omitting both layers and z_min must raise ValueError."""
        with self.assertRaises(ValueError):
            make_slab(_FakeAtoms(), [1, 0, 0])

    def test_layers_and_z_min_both_set_raises(self):
        """Specifying both layers and z_min must raise ValueError."""
        with self.assertRaises(ValueError):
            make_slab(_FakeAtoms(), [1, 0, 0], layers=2, z_min=8.0)

    def test_invalid_supercell_raises(self):
        """Non int/tuple/list and wrong-length tuples must raise ValueError."""
        with self.assertRaises(ValueError):
            surf_mod._to_surface_supercell_matrix("abc")
        with self.assertRaises(ValueError):
            surf_mod._to_surface_supercell_matrix((1, 2, 3, 4))

    # ── _to_surface_supercell_matrix ──
    def test_to_surface_supercell_matrix_int(self):
        """int n expands to diag([n, n, 1]) (in-plane n×n, z unchanged)."""
        fake_np = MagicMock()
        with patch.object(surf_mod, 'np', fake_np):
            surf_mod._to_surface_supercell_matrix(2)
        fake_np.diag.assert_called_once_with([2, 2, 1])

    def test_to_surface_supercell_matrix_2tuple(self):
        """(nx, ny) expands to diag([nx, ny, 1]) (in-plane only)."""
        fake_np = MagicMock()
        with patch.object(surf_mod, 'np', fake_np):
            surf_mod._to_surface_supercell_matrix((2, 3))
        fake_np.diag.assert_called_once_with([2, 3, 1])

    def test_to_surface_supercell_matrix_3tuple(self):
        """(nx, ny, nz) expands to diag([nx, ny, nz]) (full control)."""
        fake_np = MagicMock()
        with patch.object(surf_mod, 'np', fake_np):
            surf_mod._to_surface_supercell_matrix((2, 3, 4))
        # Source passes the tuple straight through: np.diag(supercell).
        fake_np.diag.assert_called_once_with((2, 3, 4))

    def test_to_surface_supercell_matrix_identity(self):
        """(1,1,1) yields identity; make_slab must NOT call make_supercell."""
        fake_np = MagicMock()
        fake_np.array_equal.return_value = True  # smat == eye(3)
        with patch.object(surf_mod, 'np', fake_np), \
             patch.object(surf_mod, 'general_surface') as gs, \
             patch.object(surf_mod, 'make_supercell') as ms:
            gs.surface.return_value = _FakeAtoms()
            make_slab(_FakeAtoms(), [1, 0, 0],
                      layers=2, vacuum_min=10, supercell=(1, 1, 1))
        # (1,1,1) is a 3-tuple → np.diag((1, 1, 1)) (passed through).
        fake_np.diag.assert_called_with((1, 1, 1))
        ms.assert_not_called()


class TestResolveSource(unittest.TestCase):
    """Mock-based tests for the _resolve_source helper."""

    def test_callable_source(self):
        """A callable source is invoked and its return value used."""
        atoms = _FakeAtoms(n=5)
        result = surf_mod._resolve_source(lambda: atoms)
        self.assertIs(result, atoms)

    def test_atoms_source(self):
        """An object with get_positions is returned unchanged."""
        atoms = _FakeAtoms(n=6)
        result = surf_mod._resolve_source(atoms)
        self.assertIs(result, atoms)

    def test_path_source(self):
        """A string path is dispatched to ase.io.read (mocked here)."""
        atoms = _FakeAtoms(n=8)
        with patch.object(surf_mod, 'read', return_value=atoms) as rd:
            result = surf_mod._resolve_source('/fake/path.POSCAR')
        self.assertIs(result, atoms)
        rd.assert_called_once_with('/fake/path.POSCAR', format='vasp')

    def test_path_source_cif(self):
        """A .cif path is read with format='cif'."""
        atoms = _FakeAtoms(n=8)
        with patch.object(surf_mod, 'read', return_value=atoms) as rd:
            result = surf_mod._resolve_source('/fake/struct.cif')
        self.assertIs(result, atoms)
        rd.assert_called_once_with('/fake/struct.cif', format='cif')

    def test_invalid_source_raises(self):
        """Unsupported source types must raise TypeError."""
        with self.assertRaises(TypeError):
            surf_mod._resolve_source(12345)


class TestMakeSlabsLogic(unittest.TestCase):
    """Mock-based tests for the make_slabs batch wrapper."""

    def test_returns_atoms_when_no_out_dir(self):
        """Without out_dir, returns {label: Atoms} for each miller."""
        fake_slab = _FakeAtoms(n=10)
        bulk = _FakeAtoms()
        with patch.object(surf_mod, 'make_slab', return_value=fake_slab) as ms:
            result = make_slabs(bulk, [[1, 0, 0], [1, 1, 1]], layers=2)
        self.assertEqual(len(result), 2)
        self.assertIn('surf-100', result)
        self.assertIn('surf-111', result)
        self.assertIs(result['surf-100'], fake_slab)
        self.assertIs(result['surf-111'], fake_slab)
        self.assertEqual(ms.call_count, 2)

    def test_writes_files_when_out_dir_given(self):
        """With out_dir, writes POSCAR files and returns {label: path}."""
        fake_slab = _FakeAtoms(n=10)
        bulk = _FakeAtoms()
        tmpdir = tempfile.mkdtemp()
        try:
            with patch.object(surf_mod, 'make_slab', return_value=fake_slab), \
                 patch.object(surf_mod, 'write') as w:
                result = make_slabs(bulk, [[1, 0, 0]],
                                    out_dir=tmpdir, layers=2)
            expected_path = os.path.join(tmpdir, 'surf-100.POSCAR')
            self.assertEqual(result, {'surf-100': expected_path})
            w.assert_called_once()
            # First positional arg is the path; format=vasp, direct=True.
            call_args, call_kwargs = w.call_args
            self.assertEqual(call_args[0], expected_path)
            self.assertEqual(call_kwargs.get('format'), 'vasp')
            self.assertTrue(call_kwargs.get('direct'))
        finally:
            shutil.rmtree(tmpdir)

    def test_label_format(self):
        """Miller [1,1,1] produces label 'surf-111' (no separators)."""
        fake_slab = _FakeAtoms(n=10)
        with patch.object(surf_mod, 'make_slab', return_value=fake_slab):
            result = make_slabs(_FakeAtoms(), [[1, 1, 1]], layers=2)
        self.assertEqual(list(result.keys()), ['surf-111'])

    def test_supercell_forwarded_to_make_slab(self):
        """The supercell kwarg is forwarded to each make_slab call."""
        fake_slab = _FakeAtoms(n=10)
        with patch.object(surf_mod, 'make_slab', return_value=fake_slab) as ms:
            make_slabs(_FakeAtoms(), [[1, 0, 0]],
                       layers=2, supercell=(2, 2, 1))
        _, kwargs = ms.call_args
        self.assertEqual(kwargs.get('supercell'), (2, 2, 1))


@unittest.skipUnless(HAS_ASE, "requires real ASE")
class TestMakeSlabIntegration(unittest.TestCase):
    """End-to-end tests using real ASE (skipped if ASE not installed).

    Each test runs in a fresh subprocess via ``python3 -c "..."`` so
    that real numpy/ase/scipy load cleanly without interacting with
    the mock objects installed by :func:`test_helpers.setup_mocks`.
    This sidesteps the well-known numpy-C-extension-cannot-reload
    problem that arises when test files swap mocks in and out of
    ``sys.modules`` within a single process (see test_campaign.py for
    the in-process approach and its caveats).
    """

    @staticmethod
    def _run(script):
        """Run *script* (a Python source string) in a subprocess.

        Returns stdout (stripped).  Raises AssertionError on non-zero
        exit so the calling test fails informatively.
        """
        import subprocess
        proc = subprocess.run(
            [sys.executable, '-c', script],
            capture_output=True, text=True, timeout=60)
        if proc.returncode != 0:
            raise AssertionError(
                f"subprocess failed (rc={proc.returncode}):\n"
                f"--- stdout ---\n{proc.stdout}\n"
                f"--- stderr ---\n{proc.stderr}")
        return proc.stdout.strip()

    # Common preamble for subprocess scripts: import surface.py directly
    # (bypassing extrempy/__init__.py which pulls in dpdata etc.).
    _preamble = (
        "import sys, importlib.util\n"
        f"spec = importlib.util.spec_from_file_location("
        f"'surface', {os.path.join(os.path.dirname(__file__), '..', 'extrempy', 'lazy', 'surface.py')!r})\n"
        "surf = importlib.util.module_from_spec(spec)\n"
        "spec.loader.exec_module(surf)\n"
    )

    # ── layers mode ──
    def test_fcc_al_100_layers(self):
        """FCC Al (100) with layers=4 yields a non-empty slab with vacuum."""
        out = self._run(
            self._preamble
            + "from ase.build import bulk\n"
            "bulk_al = bulk('Al', 'fcc', a=4.05, cubic=True)\n"
            "slab = surf.make_slab(bulk_al, [1, 0, 0], layers=4, vacuum_min=10)\n"
            "c = slab.cell.lengths()[2]\n"
            "assert len(slab) > 0, 'slab is empty'\n"
            "assert c > 4.05 * 4, f'c={c} too small'\n"
            "print(f'{len(slab)} {c:.2f}')\n"
        )
        n, c = out.split()
        self.assertGreater(int(n), 0)
        self.assertGreater(float(c), 4.05 * 4)

    def test_fcc_al_111_layers(self):
        """FCC Al (111) with layers=4 yields a non-empty slab."""
        out = self._run(
            self._preamble
            + "from ase.build import bulk\n"
            "bulk_al = bulk('Al', 'fcc', a=4.05, cubic=True)\n"
            "slab = surf.make_slab(bulk_al, [1, 1, 1], layers=4, vacuum_min=10)\n"
            "assert len(slab) > 0, 'slab is empty'\n"
            "print(len(slab))\n"
        )
        self.assertGreater(int(out), 0)

    # ── z_min mode ──
    def test_z_min_mode(self):
        """z_min mode produces a slab whose c-length reaches z_min."""
        out = self._run(
            self._preamble
            + "from ase.build import bulk\n"
            "bulk_al = bulk('Al', 'fcc', a=4.05, cubic=True)\n"
            "slab = surf.make_slab(bulk_al, [1, 0, 0], z_min=8.0, vacuum_min=10)\n"
            "c = slab.cell.lengths()[-1]\n"
            "assert c >= 8.0, f'c={c} < z_min=8.0'\n"
            "print(f'{c:.2f}')\n"
        )
        self.assertGreaterEqual(float(out), 8.0)

    # ── supercell ──
    def test_supercell_int(self):
        """supercell=2 quadruples the in-plane atom count (2×2×1)."""
        out = self._run(
            self._preamble
            + "from ase.build import bulk\n"
            "bulk_al = bulk('Al', 'fcc', a=4.05, cubic=True)\n"
            "base = surf.make_slab(bulk_al, [1, 0, 0], layers=2, vacuum_min=10)\n"
            "expanded = surf.make_slab(bulk_al, [1, 0, 0], layers=2, "
            "vacuum_min=10, supercell=2)\n"
            "assert len(expanded) == len(base) * 4, "
            "f'{len(expanded)} != {len(base)} * 4'\n"
            "print(f'{len(base)} {len(expanded)}')\n"
        )
        base_n, exp_n = out.split()
        self.assertEqual(int(exp_n), int(base_n) * 4)

    def test_supercell_3tuple(self):
        """supercell=(2,2,1) quadruples the in-plane atom count."""
        out = self._run(
            self._preamble
            + "from ase.build import bulk\n"
            "bulk_al = bulk('Al', 'fcc', a=4.05, cubic=True)\n"
            "base = surf.make_slab(bulk_al, [1, 0, 0], layers=2, vacuum_min=10)\n"
            "expanded = surf.make_slab(bulk_al, [1, 0, 0], layers=2, "
            "vacuum_min=10, supercell=(2, 2, 1))\n"
            "assert len(expanded) == len(base) * 4, "
            "f'{len(expanded)} != {len(base)} * 4'\n"
            "print(f'{len(base)} {len(expanded)}')\n"
        )
        base_n, exp_n = out.split()
        self.assertEqual(int(exp_n), int(base_n) * 4)

    # ── file writing ──
    def test_writes_poscar_file(self):
        """make_slabs with out_dir writes a POSCAR readable by ase.io.read."""
        tmpdir = tempfile.mkdtemp()
        try:
            out = self._run(
                self._preamble
                + "from ase.build import bulk\n"
                "from ase.io import read\n"
                "import os\n"
                "bulk_al = bulk('Al', 'fcc', a=4.05, cubic=True)\n"
                f"result = surf.make_slabs(bulk_al, [[1, 0, 0]], "
                f"out_dir={tmpdir!r}, layers=2, vacuum_min=10)\n"
                "path = result['surf-100']\n"
                "assert os.path.isfile(path), f'missing: {path}'\n"
                "reloaded = read(path)\n"
                "assert len(reloaded) > 0, 'reloaded is empty'\n"
                "print(len(reloaded))\n"
            )
            self.assertGreater(int(out), 0)
            self.assertTrue(os.path.isfile(
                os.path.join(tmpdir, 'surf-100.POSCAR')))
        finally:
            shutil.rmtree(tmpdir)


if __name__ == '__main__':
    unittest.main()
