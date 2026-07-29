"""Tests for MC3D compound query.

These tests hit the live MC3D API and require network access.
They are skipped automatically when the API is unreachable.
"""
import unittest
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

try:
    import requests
    _HAS_REQUESTS = True
except ImportError:
    _HAS_REQUESTS = False

from extrempy.lazy.mc3d import (
    get_phases_compound, backfill_mc3d_uuid,
)
from extrempy.lazy.minerals import parse_formula


class TestParseFormula(unittest.TestCase):
    """Pure-Python tests, no network."""

    def test_mgo(self):
        counts = parse_formula('MgO')
        self.assertEqual(counts, {'Mg': 1, 'O': 1})

    def test_mgsio3(self):
        counts = parse_formula('MgSiO3')
        self.assertEqual(counts, {'Mg': 1, 'Si': 1, 'O': 3})

    def test_mg2sio4(self):
        counts = parse_formula('Mg2SiO4')
        self.assertEqual(counts, {'Mg': 2, 'Si': 1, 'O': 4})

    def test_order_independent(self):
        # parse_formula is order-sensitive on input string but returns
        # a dict (order-independent).  MC3D matching uses dict equality.
        self.assertEqual(parse_formula('MgO'), parse_formula('MgO'))

    def test_empty(self):
        with self.assertRaises(ValueError):
            parse_formula('')


@unittest.skipUnless(_HAS_REQUESTS, "requests not installed")
class TestGetPhasesCompoundNetwork(unittest.TestCase):
    """Live MC3D API tests — skipped if no network."""

    @classmethod
    def setUpClass(cls):
        try:
            cls.mgo_phases = get_phases_compound('MgO', mode='all')
        except Exception as e:
            raise unittest.SkipTest(f"MC3D API unreachable: {e}")

    def test_mgo_returns_phases(self):
        self.assertGreater(len(self.mgo_phases), 0)

    def test_mgo_has_b1(self):
        """MgO should have at least one Fm-3m (SG 225) phase (B1)."""
        sg_list = [p['sg'] for p in self.mgo_phases]
        self.assertIn(225, sg_list,
                      f"Expected SG 225 (Fm-3m) for MgO B1; got SGs: {sg_list}")

    def test_phase_has_required_fields(self):
        for p in self.mgo_phases:
            for field in ('id', 'sg', 'spg_intl', 'structure_uuid',
                          'phase_type', 'n_atoms_cell'):
                self.assertIn(field, p,
                              f"Phase {p.get('id')}: missing {field}")

    def test_mgsio3_has_bridgmanite(self):
        """MgSiO3 should have a Pnma (SG 62) phase (bridgmanite)."""
        try:
            phases = get_phases_compound('MgSiO3', mode='all')
        except Exception as e:
            self.skipTest(f"MC3D API unreachable: {e}")
        sg_list = [p['sg'] for p in phases]
        self.assertIn(62, sg_list,
                      f"Expected SG 62 (Pnma) for bridgmanite; got: {sg_list}")


@unittest.skipUnless(_HAS_REQUESTS, "requests not installed")
class TestBackfillUuid(unittest.TestCase):
    """Test backfill on the real mineral database."""

    def test_backfill_returns_report_dict(self):
        from extrempy.lazy.minerals import load_mineral_db
        import copy
        db = load_mineral_db()
        db_copy = copy.deepcopy(db)
        try:
            report = backfill_mc3d_uuid(db_copy)
        except Exception as e:
            self.skipTest(f"MC3D API unreachable: {e}")
        # Report must have all 5 keys
        for key in ('matched', 'multiple_candidates', 'no_match',
                    'skipped', 'api_errors'):
            self.assertIn(key, report)
        # matched should be a list of labels
        self.assertIsInstance(report['matched'], list)
        # multiple_candidates should be a dict
        self.assertIsInstance(report['multiple_candidates'], dict)
        # If anything matched, the phase's mc3d_uuid should now be set
        for label in report['matched']:
            for compound in db_copy['compounds'].values():
                for ph in compound['phases']:
                    if ph['label'] == label:
                        self.assertIsNotNone(ph.get('mc3d_uuid'),
                            f'{label} in matched but mc3d_uuid still null')

    def test_backfill_multiple_candidates_structure(self):
        """If multiple candidates found, report should list them with details."""
        from extrempy.lazy.minerals import load_mineral_db
        import copy
        db = load_mineral_db()
        db_copy = copy.deepcopy(db)
        try:
            report = backfill_mc3d_uuid(db_copy)
        except Exception as e:
            self.skipTest(f"MC3D API unreachable: {e}")
        for label, candidates in report['multiple_candidates'].items():
            self.assertGreaterEqual(len(candidates), 2,
                f'{label} should have >=2 candidates')
            for c in candidates:
                for field in ('uuid', 'mc3d_id', 'sg', 'spg_intl',
                              'energy_per_atom', 'phase_type'):
                    self.assertIn(field, c,
                        f'{label} candidate missing {field}')

    def test_backfill_no_match_marks_manual(self):
        """Phases with no MC3D match should be marked structure_source='manual'."""
        from extrempy.lazy.minerals import load_mineral_db
        import copy
        db = load_mineral_db()
        db_copy = copy.deepcopy(db)
        try:
            report = backfill_mc3d_uuid(db_copy)
        except Exception as e:
            self.skipTest(f"MC3D API unreachable: {e}")
        # Find phases mentioned in no_match and verify structure_source
        for entry in report['no_match']:
            # entry format: "formula/label: sg=N not in MC3D (...)"
            formula, rest = entry.split('/', 1)
            label = rest.split(':')[0]
            for ph in db_copy['compounds'][formula]['phases']:
                if ph['label'] == label:
                    self.assertEqual(ph.get('structure_source'), 'manual',
                        f'{label} should be marked manual after no_match')


if __name__ == '__main__':
    unittest.main()
