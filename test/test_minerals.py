import unittest
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from extrempy.lazy.minerals import (
    parse_formula, formula_elements, load_mineral_db,
    list_compounds, get_compound, get_compound_phases,
    make_extreme_segments, lookup_earth_layers, validate_db,
    EARTH_LAYERS, DEEP_EARTH_ELEMENTS,
    GLOBAL_P_RANGE_GPa, GLOBAL_T_RANGE_K,
)


class TestParseFormula(unittest.TestCase):
    def test_single_element(self):
        self.assertEqual(parse_formula('Mg'), {'Mg': 1})

    def test_binary(self):
        self.assertEqual(parse_formula('MgO'), {'Mg': 1, 'O': 1})

    def test_with_counts(self):
        self.assertEqual(parse_formula('MgSiO3'),
                         {'Mg': 1, 'Si': 1, 'O': 3})

    def test_multi_digit_count(self):
        self.assertEqual(parse_formula('Mg2SiO4'),
                         {'Mg': 2, 'Si': 1, 'O': 4})

    def test_complex(self):
        self.assertEqual(parse_formula('NaAlSi3O8'),
                         {'Na': 1, 'Al': 1, 'Si': 3, 'O': 8})

    def test_formula_elements_sorted(self):
        self.assertEqual(formula_elements('Mg2SiO4'),
                         ['Mg', 'O', 'Si'])

    def test_invalid_char(self):
        with self.assertRaises(ValueError):
            parse_formula('Mg#O')

    def test_empty(self):
        with self.assertRaises(ValueError):
            parse_formula('')


class TestDatabaseLoading(unittest.TestCase):
    def test_load_default(self):
        db = load_mineral_db()
        self.assertIn('compounds', db)
        self.assertIn('MgO', db['compounds'])

    def test_has_16_elements(self):
        db = load_mineral_db()
        for el in DEEP_EARTH_ELEMENTS:
            self.assertIn(el, db['compounds'],
                          f'{el} missing from mineral_phases.json')

    def test_has_mg_sio_compounds(self):
        db = load_mineral_db()
        for c in ['MgO', 'SiO2', 'MgSiO3', 'Mg2SiO4']:
            self.assertIn(c, db['compounds'])

    def test_list_compounds(self):
        compounds = list_compounds()
        self.assertIsInstance(compounds, list)
        self.assertGreaterEqual(len(compounds), 20)
        self.assertIn('MgO', compounds)


class TestQuery(unittest.TestCase):
    def test_get_compound(self):
        c = get_compound('MgO')
        self.assertEqual(c['Tm_at_1bar_K'], 3098)
        self.assertGreaterEqual(len(c['phases']), 2)

    def test_get_compound_missing(self):
        with self.assertRaises(KeyError):
            get_compound('XxYz')

    def test_get_phases(self):
        phases = get_compound_phases('MgSiO3')
        labels = [p['label'] for p in phases]
        self.assertIn('MgSiO3-bridgmanite', labels)
        self.assertIn('MgSiO3-postperovskite', labels)

    def test_phase_has_required_fields(self):
        for formula in list_compounds():
            for ph in get_compound_phases(formula):
                for field in ('label', 'P_range_GPa', 'T_range_K',
                              'structure_source', 'refs',
                              'spacegroup', 'spg_intl'):
                    self.assertIn(field, ph,
                                  f'{formula}/{ph.get("label")}: missing {field}')

    def test_non_liq_mc3d_phases_have_spacegroup(self):
        """All non-LIQ phases with structure_source='mc3d' must have spacegroup."""
        for formula in list_compounds():
            for ph in get_compound_phases(formula):
                if ph['label'].endswith('-LIQ'):
                    self.assertIsNone(ph.get('spacegroup'),
                        f'{formula}/{ph["label"]}: LIQ should have spacegroup=null')
                    continue
                if ph.get('structure_source') == 'mc3d':
                    self.assertIsNotNone(ph.get('spacegroup'),
                        f'{formula}/{ph["label"]}: mc3d phase missing spacegroup')

    def test_liq_phases_have_null_sg(self):
        """LIQ phases should have spacegroup=null and spg_intl=null."""
        for formula in list_compounds():
            for ph in get_compound_phases(formula):
                if ph['label'].endswith('-LIQ'):
                    self.assertIsNone(ph['spacegroup'],
                        f'{formula}/{ph["label"]}: LIQ spacegroup should be null')
                    self.assertIsNone(ph['spg_intl'],
                        f'{formula}/{ph["label"]}: LIQ spg_intl should be null')


class TestEarthLayers(unittest.TestCase):
    def test_upper_mantle(self):
        layers = lookup_earth_layers((0, 10))
        self.assertIn('upper_mantle', layers)

    def test_lower_mantle(self):
        layers = lookup_earth_layers((30, 100))
        self.assertIn('lower_mantle', layers)

    def test_inner_core(self):
        layers = lookup_earth_layers((340, 400))
        self.assertIn('inner_core', layers)

    def test_multi_layer(self):
        layers = lookup_earth_layers((10, 30))
        self.assertIn('upper_mantle', layers)
        self.assertIn('transition_zone', layers)
        self.assertIn('lower_mantle', layers)

    def test_no_match(self):
        layers = lookup_earth_layers((-1, -0.5))
        self.assertEqual(layers, [])


class TestMakeSegments(unittest.TestCase):
    def test_mgo_segments(self):
        segs = make_extreme_segments('MgO')
        self.assertGreaterEqual(len(segs), 2)
        for s in segs:
            self.assertIn('label', s)
            self.assertIn('P_core', s)
            self.assertIn('T_core', s)
            self.assertIn('P_explore', s)
            self.assertIn('T_explore', s)
            self.assertIn('earth_layers', s)
            self.assertIn('refs', s)
            # 2D segs should carry spacegroup + spg_intl from the database
            self.assertIn('spacegroup', s)
            self.assertIn('spg_intl', s)

    def test_explore_wider_than_core_P(self):
        segs = make_extreme_segments('MgO')
        for s in segs:
            self.assertLessEqual(s['P_explore'][0], s['P_core'][0],
                f'{s["label"]}: P_explore lower > P_core lower')
            self.assertGreaterEqual(s['P_explore'][1], s['P_core'][1],
                f'{s["label"]}: P_explore upper < P_core upper')

    def test_explore_wider_than_core_T(self):
        segs = make_extreme_segments('MgSiO3')
        for s in segs:
            self.assertLessEqual(s['T_explore'][0], s['T_core'][0],
                f'{s["label"]}: T_explore lower > T_core lower')
            self.assertGreaterEqual(s['T_explore'][1], s['T_core'][1],
                f'{s["label"]}: T_explore upper < T_core upper')

    def test_liq_label(self):
        segs = make_extreme_segments('Mg2SiO4')
        liq_segs = [s for s in segs if s['label'].endswith('-LIQ')]
        self.assertGreater(len(liq_segs), 0,
                          'Should have a LIQ phase')

    def test_custom_margin(self):
        segs_default = make_extreme_segments('MgO', p_margin=0.15)
        segs_wide = make_extreme_segments('MgO', p_margin=0.30)
        # wider margin → wider explore
        self.assertGreater(segs_wide[0]['P_explore'][1] - segs_wide[0]['P_explore'][0],
                          segs_default[0]['P_explore'][1] - segs_default[0]['P_explore'][0])

    def test_earth_layers_populated(self):
        segs = make_extreme_segments('Fe')
        for s in segs:
            # Fe-hcp-epsilon spans 13-400 GPa → multiple layers
            if s['label'] == 'Fe-hcp-epsilon':
                self.assertGreater(len(s['earth_layers']), 1)

    def test_bridgmanite_in_lower_mantle(self):
        segs = make_extreme_segments('MgSiO3')
        bv = [s for s in segs if s['label'] == 'MgSiO3-bridgmanite'][0]
        self.assertIn('lower_mantle', bv['earth_layers'])


class TestValidate(unittest.TestCase):
    def test_validate_returns_list(self):
        warnings = validate_db()
        self.assertIsInstance(warnings, list)

    def test_no_missing_required_fields(self):
        """No warning should mention missing P_range, T_range, or spacegroup."""
        warnings = validate_db()
        field_warnings = [w for w in warnings
                          if 'missing P_range' in w
                          or 'missing T_range' in w
                          or 'spacegroup is null' in w]
        self.assertEqual(field_warnings, [],
                         f'Fields missing: {field_warnings}')

    def test_mc3d_uuid_null_warns(self):
        """mc3d_uuid=null with structure_source='mc3d' should warn (before backfill)."""
        warnings = validate_db()
        uuid_warnings = [w for w in warnings if 'mc3d_uuid is null' in w]
        self.assertGreater(len(uuid_warnings), 0,
                          'Should warn about null mc3d_uuid before backfill')


class TestGlobalRanges(unittest.TestCase):
    def test_global_p_range(self):
        self.assertEqual(GLOBAL_P_RANGE_GPa, (0, 400))

    def test_global_t_range(self):
        self.assertEqual(GLOBAL_T_RANGE_K, (300, 10000))

    def test_16_elements(self):
        self.assertEqual(len(DEEP_EARTH_ELEMENTS), 16)


if __name__ == '__main__':
    unittest.main()
