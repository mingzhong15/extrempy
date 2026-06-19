import unittest
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

# Pure Python imports — no numpy/dpdata needed for lib.py
from extrempy.lazy.lib import (ELEMENT_PHASE_DATA, ELEMENTS_BY_STRUCTURE,
    get_phase_segments, get_viable_elements, _compute_overlap,
    _deduplicate_labels, _get_mass_map, SUPPORTED_STRUCTURES,
    SKIP_MAGNETIC, SKIP_GAS, SKIP_LANTHANIDE, SKIP_RADIOACTIVE,
    SKIP_ACTINIDE)


class TestViableElements(unittest.TestCase):
    def test_count(self):
        viable = get_viable_elements()
        self.assertGreaterEqual(len(viable), 35)
        self.assertLessEqual(len(viable), 50)

    def test_contains_key_elements(self):
        viable = set(get_viable_elements())
        for el in ['Ti', 'Al', 'Cu', 'W', 'Li', 'Be', 'Si', 'Ge']:
            self.assertIn(el, viable, f'{el} should be viable')

    def test_excludes_magnetic(self):
        viable = set(get_viable_elements())
        for el in SKIP_MAGNETIC:
            self.assertNotIn(el, viable, f'{el} (magnetic) should be excluded')

    def test_excludes_gas(self):
        viable = set(get_viable_elements())
        for el in SKIP_GAS:
            self.assertNotIn(el, viable, f'{el} (gas) should be excluded')

    def test_excludes_radioactive(self):
        viable = set(get_viable_elements())
        for el in SKIP_RADIOACTIVE:
            self.assertNotIn(el, viable, f'{el} (radioactive) should be excluded')

    def test_excludes_lanthanide(self):
        viable = set(get_viable_elements())
        for el in SKIP_LANTHANIDE:
            self.assertNotIn(el, viable, f'{el} (lanthanide) should be excluded')

    def test_excludes_orthorhombic_rhombohedral(self):
        viable = set(get_viable_elements())
        # Ga/Sb/Bi have complex struct that ASE bulk can't handle
        for el in ['Ga', 'Sb', 'Bi']:
            self.assertNotIn(el, viable,
                             f'{el} (complex struct) should be excluded')

    def test_all_viable_have_tm(self):
        for el in get_viable_elements():
            self.assertIsInstance(ELEMENT_PHASE_DATA[el].get('Tm'),
                                  (int, float),
                                  f'{el} should have numeric Tm')

    def test_all_viable_have_rt_structure(self):
        for el in get_viable_elements():
            rt = ELEMENT_PHASE_DATA[el].get('rt_structure')
            self.assertIn(rt, SUPPORTED_STRUCTURES,
                          f'{el} rt_structure={rt} not supported')

    def test_all_viable_have_mass(self):
        viable = get_viable_elements()
        mass_list = _get_mass_map(viable)
        self.assertEqual(len(mass_list), len(viable))

    def test_return_type(self):
        viable = get_viable_elements()
        self.assertIsInstance(viable, list)
        for el in viable:
            self.assertIsInstance(el, str)

    def test_sorted_order(self):
        viable = get_viable_elements()
        for i in range(len(viable) - 1):
            self.assertLessEqual(
                ELEMENT_PHASE_DATA[viable[i]]['z'],
                ELEMENT_PHASE_DATA[viable[i+1]]['z'],
                'Elements should be sorted by atomic number')


class TestPhaseSegments(unittest.TestCase):
    def test_ti_three_phases(self):
        segs = get_phase_segments('Ti')
        self.assertEqual(len(segs), 3)
        labels = [s['label'] for s in segs]
        self.assertEqual(labels[0], 'Ti-HCP')
        self.assertEqual(labels[1], 'Ti-BCC')
        self.assertTrue(segs[2]['label'].endswith('-LIQ'))

    def test_li_two_phases_alpha_discarded(self):
        segs = get_phase_segments('Li')
        self.assertEqual(len(segs), 2)
        # alpha (0-80K) should be discarded by drop_below_T=200
        self.assertEqual(segs[0]['label'], 'Li-BCC')

    def test_sn_two_phases_alpha_discarded(self):
        segs = get_phase_segments('Sn')
        self.assertEqual(len(segs), 2)
        # alpha diamond (0-286K) should be discarded by T_range intersection
        self.assertEqual(segs[0]['label'], 'Sn-BCT')

    def test_be_three_phases(self):
        segs = get_phase_segments('Be')
        self.assertEqual(len(segs), 3)
        self.assertEqual(segs[0]['label'], 'Be-HCP')
        self.assertEqual(segs[1]['label'], 'Be-BCC')
        self.assertEqual(segs[2]['label'], 'Be-LIQ')
        # BCC is narrow (1527-1560K)
        core_low, core_high = segs[1]['T_core']
        self.assertAlmostEqual(core_low, 1527, delta=1)
        self.assertAlmostEqual(core_high, 1560, delta=1)

    def test_al_two_phases(self):
        segs = get_phase_segments('Al')
        self.assertEqual(len(segs), 2)
        self.assertEqual(segs[0]['label'], 'Al-FCC')
        self.assertTrue(segs[1]['label'].endswith('-LIQ'))

    def test_w_two_phases(self):
        segs = get_phase_segments('W')
        self.assertEqual(len(segs), 2)
        self.assertEqual(segs[0]['label'], 'W-BCC')
        self.assertTrue(segs[1]['label'].endswith('-LIQ'))
        # W has high Tm
        self.assertGreater(segs[0]['T_core'][1], 3000)
        self.assertGreater(segs[1]['T_explore'][1], 6000)

    def test_no_duplicate_labels(self):
        for el in get_viable_elements():
            segs = get_phase_segments(el)
            labels = [s['label'] for s in segs]
            self.assertEqual(len(labels), len(set(labels)),
                             f'{el} has duplicate labels: {labels}')

    def test_explore_bounds_in_range(self):
        for el in get_viable_elements():
            segs = get_phase_segments(el)
            Tm = ELEMENT_PHASE_DATA[el]['Tm']
            T_max = 2 * Tm
            for s in segs:
                lo, hi = s['T_explore']
                self.assertGreaterEqual(lo, 300,
                    f'{el} {s["label"]} T_explore lower {lo} < 300')
                self.assertLessEqual(hi, T_max,
                    f'{el} {s["label"]} T_explore upper {hi} > 2*Tm')

    def test_liquid_last_seg(self):
        for el in get_viable_elements():
            segs = get_phase_segments(el)
            self.assertTrue(segs[-1]['label'].endswith('-LIQ'),
                           f'{el} last segment should be LIQ')

    def test_explore_overlaps_core(self):
        for el in get_viable_elements():
            segs = get_phase_segments(el)
            for s in segs:
                # explore should extend BEYOND core (due to overlap)
                lo, hi = s['T_explore']
                core_lo, core_hi = s['T_core']
                self.assertLessEqual(lo, core_lo,
                    f'{el} {s["label"]} explore lower ({lo}) > core lower ({core_lo})')
                self.assertGreaterEqual(hi, core_hi,
                    f'{el} {s["label"]} explore upper ({hi}) < core upper ({core_hi})')

    def test_error_no_tm(self):
        with self.assertRaises(ValueError):
            get_phase_segments('C')  # C has Tm=None

    def test_error_unknown_element(self):
        with self.assertRaises(ValueError):
            get_phase_segments('Xx')

    def test_ti_hcp_core_range(self):
        segs = get_phase_segments('Ti')
        self.assertEqual(segs[0]['T_core'], (300, 1155))

    def test_ti_bcc_core_range(self):
        segs = get_phase_segments('Ti')
        self.assertEqual(segs[1]['T_core'], (1155, 1941))


class TestComputeOverlap(unittest.TestCase):
    def test_auto_below_500(self):
        # max(0.2*100, 100) = 100
        self.assertEqual(_compute_overlap(100, False, 'auto'), 100)

    def test_auto_above_500(self):
        # max(0.2*1941, 100) = max(388, 100) = 388
        self.assertAlmostEqual(_compute_overlap(1941, True, 'auto'), 388.2)

    def test_auto_melting_uses_Tm(self):
        # same formula regardless of is_melting with 'auto'
        self.assertAlmostEqual(_compute_overlap(1155, True, 'auto'),
                               max(0.2*1155, 100))
        self.assertAlmostEqual(_compute_overlap(1155, False, 'auto'),
                               max(0.2*1155, 100))

    def test_callable_rule(self):
        result = _compute_overlap(1000, False, lambda T, m: T * 0.1)
        self.assertEqual(result, 100.0)

    def test_float_rule(self):
        result = _compute_overlap(1000, False, 200)
        self.assertEqual(result, 200.0)

    def test_low_temperature_uses_min_100(self):
        # max(0.2*300, 100) = max(60, 100) = 100
        self.assertEqual(_compute_overlap(300, False, 'auto'), 100)


class TestDeduplicate(unittest.TestCase):
    def test_no_collision(self):
        segs = [{'label': 'Ti-HCP'}, {'label': 'Ti-BCC'}]
        _deduplicate_labels(segs)
        self.assertEqual(segs[0]['label'], 'Ti-HCP')
        self.assertEqual(segs[1]['label'], 'Ti-BCC')

    def test_single_collision(self):
        segs = [{'label': 'Li-BCC'}, {'label': 'Li-BCC'}]
        _deduplicate_labels(segs)
        self.assertEqual(segs[0]['label'], 'Li-BCC')
        self.assertEqual(segs[1]['label'], 'Li-BCC-2')

    def test_triple_collision(self):
        segs = [{'label': 'A'}, {'label': 'A'}, {'label': 'A'}]
        _deduplicate_labels(segs)
        self.assertEqual(segs[0]['label'], 'A')
        self.assertEqual(segs[1]['label'], 'A-2')
        self.assertEqual(segs[2]['label'], 'A-3')

    def test_mixed_collision(self):
        segs = [{'label': 'A'}, {'label': 'B'}, {'label': 'A'}]
        _deduplicate_labels(segs)
        self.assertEqual(segs[0]['label'], 'A')
        self.assertEqual(segs[1]['label'], 'B')
        self.assertEqual(segs[2]['label'], 'A-2')


if __name__ == '__main__':
    unittest.main()
