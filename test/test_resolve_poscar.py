"""Tests for extrempy.structure.resolve_poscar and the source factories.

Uses test_helpers.setup_mocks() to stub heavy deps (numpy, dpdata,
ase, scipy) so the tests run without those packages installed.  The
source callables are exercised with light mock Atoms-like objects.
"""
import sys, os, tempfile, shutil, unittest
from unittest.mock import MagicMock

# Stub heavy deps before importing extrempy.
from test_helpers import setup_mocks
setup_mocks()
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from extrempy.structure import (
    resolve_poscar, ase_source, mc3d_source, DEFAULT_SUPERCELL,
)


class _FakeAtoms:
    """Lightweight stand-in for ase.Atoms (has get_positions, __len__)."""
    def __init__(self, n=4):
        self._n = n
    def __len__(self):
        return self._n
    def get_positions(self):
        return [(0, 0, 0)] * self._n


class TestResolvePoscar(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.confs = os.path.join(self.tmpdir, 'confs')

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    # ── LIQ label ──
    def test_liq_returns_none_and_writes_nothing(self):
        src = MagicMock(return_value=_FakeAtoms())
        p = resolve_poscar('Al', 'Al-LIQ', confs_dir=self.confs,
                           source=src, verbose=False)
        self.assertIsNone(p)
        src.assert_not_called()  # must short-circuit before calling source
        self.assertFalse(os.path.exists(self.confs))

    # ── existing file → skip ──
    def test_existing_skipped_without_force(self):
        target = os.path.join(self.confs, 'Al-FCC.POSCAR')
        os.makedirs(self.confs)
        open(target, 'w').close()  # empty file
        src = MagicMock(return_value=_FakeAtoms())
        p = resolve_poscar('Al', 'Al-FCC', confs_dir=self.confs,
                           source=src, force=False, verbose=False)
        self.assertEqual(p, target)
        src.assert_not_called()  # must not rebuild when file exists

    def test_force_overwrites_existing(self):
        target = os.path.join(self.confs, 'Al-FCC.POSCAR')
        os.makedirs(self.confs)
        open(target, 'w').close()
        src = MagicMock(return_value=_FakeAtoms(n=8))
        p = resolve_poscar('Al', 'Al-FCC', confs_dir=self.confs,
                           source=src, force=True, verbose=False)
        self.assertEqual(p, target)
        src.assert_called_once()

    # ── source: callable ──
    def test_callable_source_writes_file(self):
        atoms = _FakeAtoms(n=4)
        # ase.io.write is mocked via sys.modules — just ensure no exception
        # and the source is invoked exactly once.
        p = resolve_poscar('Al', 'Al-FCC', confs_dir=self.confs,
                           source=lambda: atoms, verbose=False)
        self.assertEqual(p, os.path.join(self.confs, 'Al-FCC.POSCAR'))

    # ── source: Atoms object ──
    def test_atoms_source(self):
        atoms = _FakeAtoms(n=6)
        p = resolve_poscar('Al', 'Al-BCC', confs_dir=self.confs,
                           source=atoms, verbose=False)
        self.assertEqual(p, os.path.join(self.confs, 'Al-BCC.POSCAR'))

    # ── source: path str ──
    def test_path_source_missing_raises(self):
        with self.assertRaises(FileNotFoundError):
            resolve_poscar('Al', 'Al-FCC', confs_dir=self.confs,
                           source=os.path.join(self.tmpdir, 'nope.POSCAR'),
                           verbose=False)

    # ── source: invalid type ──
    def test_invalid_source_type_raises(self):
        with self.assertRaises(TypeError):
            resolve_poscar('Al', 'Al-FCC', confs_dir=self.confs,
                           source=12345, verbose=False)


class TestSourceFactories(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_ase_source_factory_signature(self):
        s = ase_source('Al', structure_type='fcc', supercell=(2, 2, 2))
        self.assertTrue(callable(s))
        # _generate_atoms is real (only numpy is mocked) — but bulk() is
        # mocked via sys.modules['ase'].build.bulk, so calling s() will
        # engage the mock chain. Just verify it doesn't crash on shutdown.

    def test_mc3d_source_factory_signature(self):
        s = mc3d_source('uuid-abc', target_atoms=100)
        self.assertTrue(callable(s))

    def test_default_supercell_table_completeness(self):
        # Every structure type DPBuilder/EOS relies on must have a default.
        for st in ('fcc', 'bcc', 'hcp', 'diamond', 'dhcp', 'sc', 'bct'):
            self.assertIn(st, DEFAULT_SUPERCELL, f'{st} missing from table')
            sc = DEFAULT_SUPERCELL[st]
            self.assertEqual(len(sc), 3)
            self.assertTrue(all(isinstance(x, int) and x >= 1 for x in sc))


class TestDpBuilderUsesResolvePoscar(unittest.TestCase):
    """Smoke test: ElementDPBuilder.generate_poscars dispatches to
    resolve_poscar with the right source type per seg['structure']."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_ase_seg_uses_ase_source(self):
        from extrempy.campaign.single_element_dp import ElementDPBuilder
        import extrempy.campaign.single_element_dp as sed

        b = ElementDPBuilder(
            'Al', work_root=self.tmpdir, potcar_lib=self.tmpdir,
            potcar_set='PBE54')
        segs = [{'label': 'Al-FCC', 'structure': 'fcc',
                 'T_core': (300, 600), 'T_explore': (200, 1200)}]

        calls = []
        orig = sed.resolve_poscar if hasattr(sed, 'resolve_poscar') else None

        # Patch resolve_poscar where generate_poscars imports it from
        # (extrempy.structure.resolve_poscar).
        import extrempy.structure as struct_mod
        saved = struct_mod.resolve_poscar

        def fake_resolve(element, label, *, confs_dir, source, **kw):
            calls.append((element, label, type(source).__name__
                          if not callable(source) else 'callable'))
            return os.path.join(confs_dir, f'{label}.POSCAR')

        struct_mod.resolve_poscar = fake_resolve
        try:
            b.generate_poscars(segs)
        finally:
            struct_mod.resolve_poscar = saved

        self.assertEqual(len(calls), 1)
        self.assertEqual(calls[0][0], 'Al')
        self.assertEqual(calls[0][1], 'Al-FCC')
        # ase_source returns a closure (callable), not an Atoms object.
        self.assertTrue(callable(calls[0][2]) or calls[0][2] == 'callable')

    def test_liq_seg_is_skipped_by_resolve_poscar(self):
        from extrempy.campaign.single_element_dp import ElementDPBuilder

        b = ElementDPBuilder(
            'Al', work_root=self.tmpdir, potcar_lib=self.tmpdir,
            potcar_set='PBE54')
        segs = [
            {'label': 'Al-FCC', 'structure': 'fcc',
             'T_core': (300, 600), 'T_explore': (200, 1200)},
            {'label': 'Al-LIQ', 'structure': 'fcc',
             'T_core': (600, 1200), 'T_explore': (200, 1200)},
        ]
        import extrempy.structure as struct_mod
        saved = struct_mod.resolve_poscar
        calls = []

        def fake_resolve(element, label, *, confs_dir, source, **kw):
            calls.append(label)
            if label.endswith('-LIQ'):
                return None  # resolve_poscar returns None for LIQ
            return os.path.join(confs_dir, f'{label}.POSCAR')

        struct_mod.resolve_poscar = fake_resolve
        try:
            b.generate_poscars(segs)
        finally:
            struct_mod.resolve_poscar = saved
        # LIQ seg still goes through resolve_poscar (it returns None there),
        # which matches the unified contract.
        self.assertIn('Al-LIQ', calls)


class TestEosFindPoscarRoleBased(unittest.TestCase):
    """Verify EOSCalculator._find_poscar dispatches by role and falls
    back liquid → solid_rt silently."""

    def _make_calc(self, element='Al'):
        from extrempy.campaign.melt import ElementEOSCalculator
        return ElementEOSCalculator(
            element, work_root=self.tmpdir, dpgen_dir=None)

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        # Create an Al-fcc POSCAR in poscar_dir.
        self.poscar_dir = os.path.join(self.tmpdir, 'poscars')
        os.makedirs(self.poscar_dir)
        self.solid_path = os.path.join(self.poscar_dir, 'Al-FCC.POSCAR')
        open(self.solid_path, 'w').close()
        self.liquid_path = os.path.join(self.poscar_dir, 'Al-LIQ.POSCAR')
        open(self.liquid_path, 'w').close()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_solid_rt_resolves_by_label(self):
        calc = self._make_calc()
        calc.poscar_dir = self.poscar_dir
        # Al's rt_structure is 'fcc' per ELEMENT_PHASE_DATA.
        p = calc._find_poscar('solid_rt')
        self.assertEqual(p, self.solid_path)

    def test_liquid_resolves_when_liq_exists(self):
        calc = self._make_calc()
        calc.poscar_dir = self.poscar_dir
        p = calc._find_poscar('liquid')
        self.assertEqual(p, self.liquid_path)

    def test_liquid_falls_back_to_solid_silently(self):
        calc = self._make_calc()
        calc.poscar_dir = self.poscar_dir
        os.remove(self.liquid_path)  # LIQ missing
        p = calc._find_poscar('liquid')
        self.assertEqual(p, self.solid_path)  # silent fallback


if __name__ == '__main__':
    unittest.main()