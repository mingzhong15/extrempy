"""Tests for extrempy.campaign.chunk.

Loaded directly from file path to avoid the heavy extrempy/__init__.py
import chain.  numpy is the only real dependency needed.
"""

import os
import sys
import importlib.util
import tempfile
import unittest

# Load chunk.py directly from its file path (bypasses package __init__).
_CHUNK_PATH = os.path.join(os.path.dirname(__file__), '..',
                           'extrempy', 'campaign', 'chunk.py')
_spec = importlib.util.spec_from_file_location('extrempy.campaign.chunk', _CHUNK_PATH)
chunk = importlib.util.module_from_spec(_spec)
sys.modules['extrempy.campaign.chunk'] = chunk
_spec.loader.exec_module(chunk)

from extrempy.campaign.chunk import (
    read_chunk_profile,
    diagnose_block,
    diagnose_rdf,
    cross_validate_with_rdf,
    diagnose_case,
    Q4_SOLID, Q6_SOLID, DELTA_Q4, DELTA_RHO,
)


def _make_block(q4_lower, q4_upper, q6_lower=0.4, q6_upper=0.0,
                rho_lower=2.7, rho_upper=2.5, n=10):
    """Build a synthetic block with n chunks split into halves."""
    half = n // 2
    chunks = []
    for i in range(half):
        chunks.append({'id': i + 1, 'coord': 0.05 + 0.1 * i,
                       'q4': q4_lower, 'q6': q6_lower,
                       'density': rho_lower})
    for i in range(half, n):
        chunks.append({'id': i + 1, 'coord': 0.05 + 0.1 * i,
                       'q4': q4_upper, 'q6': q6_upper,
                       'density': rho_upper})
    return {'timestep': 1000, 'chunks': chunks}


def _write_chunk_file(path, blocks, with_virial=True):
    """Write a synthetic chunk.profile with the LAMMPS header format."""
    header_cols = ['Chunk', 'Coord1', 'Ncount', 'density/mass', 'temp']
    if with_virial:
        header_cols.append('v_virial_atom')
    header_cols += ['c_Q[1]', 'c_Q[2]']

    with open(path, 'w') as f:
        for blk in blocks:
            f.write(f"{blk['timestep']}\n")
            f.write('# ' + ' '.join(header_cols) + '\n')
            for c in blk['chunks']:
                row = [str(c['id']), f"{c['coord']:.4f}", '100',
                       f"{c.get('density', 2.7):.4f}", f"{c.get('temp', 300):.1f}"]
                if with_virial:
                    row.append(f"{c.get('virial', 0.0):.4f}")
                row.append(f"{c.get('q4', 0.0):.4f}")
                row.append(f"{c.get('q6', 0.0):.4f}")
                f.write(' '.join(row) + '\n')


class TestReadChunkProfile(unittest.TestCase):
    """Header-driven parsing (C1): column reordering must not break it."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmpdir)

    def test_parses_standard_format(self):
        path = os.path.join(self.tmpdir, 'chunk.profile')
        blocks = [_make_block(0.18, 0.02)]
        _write_chunk_file(path, blocks)
        result = read_chunk_profile(path)
        self.assertEqual(len(result), 1)
        self.assertEqual(len(result[0]['chunks']), 10)
        c0 = result[0]['chunks'][0]
        self.assertIn('q4', c0)
        self.assertIn('q6', c0)
        self.assertIn('density', c0)
        self.assertAlmostEqual(c0['q4'], 0.18, places=3)

    def test_column_reorder_safe(self):
        """C1: even if columns are reordered, canonical keys are correct."""
        path = os.path.join(self.tmpdir, 'chunk.profile')
        # Custom header with reordered columns
        with open(path, 'w') as f:
            f.write('1000\n')
            f.write('# Chunk Coord1 Ncount c_Q[1] c_Q[2] density/mass temp\n')
            f.write('1 0.0500 100 0.1800 0.4000 2.7000 300.0\n')
            f.write('2 0.1500 100 0.1800 0.4000 2.7000 300.0\n')
        result = read_chunk_profile(path)
        c0 = result[0]['chunks'][0]
        self.assertAlmostEqual(c0['q4'], 0.18, places=3)
        self.assertAlmostEqual(c0['q6'], 0.40, places=3)
        self.assertAlmostEqual(c0['density'], 2.70, places=3)

    def test_missing_file_returns_empty(self):
        self.assertEqual(read_chunk_profile('/nonexistent/path'), [])

    def test_multiple_blocks(self):
        path = os.path.join(self.tmpdir, 'chunk.profile')
        blocks = [_make_block(0.18, 0.02, ),
                  _make_block(0.18, 0.02)]
        blocks[1]['timestep'] = 2000
        _write_chunk_file(path, blocks)
        result = read_chunk_profile(path)
        self.assertEqual(len(result), 2)
        self.assertEqual(result[0]['timestep'], 1000)
        self.assertEqual(result[1]['timestep'], 2000)

    def test_header_with_coordination_first_word(self):
        """B: newer LAMMPS may emit 'Coordination' instead of 'Chunk' as
        the leading header word, but the data still has a Chunk ID column.
        The relaxed recogniser matches on any of chunk/coord1/coord/ncount.
        """
        path = os.path.join(self.tmpdir, 'chunk.profile')
        with open(path, 'w') as f:
            f.write('1000\n')
            f.write('# Coordination Coord1 Ncount density/mass temp c_Q[1] c_Q[2]\n')
            f.write('1 0.0500 100 2.7000 300.0 0.1800 0.4000\n')
            f.write('2 0.1500 100 2.7000 300.0 0.1800 0.4000\n')
        result = read_chunk_profile(path)
        self.assertEqual(len(result), 1)
        c0 = result[0]['chunks'][0]
        self.assertIn('q4', c0)
        self.assertAlmostEqual(c0['q4'], 0.18, places=3)
        self.assertAlmostEqual(c0['q6'], 0.40, places=3)


class TestDiagnoseBlock(unittest.TestCase):
    """Per-block verdict logic."""

    def test_both_solid(self):
        blk = _make_block(q4_lower=0.18, q4_upper=0.19,
                          q6_lower=0.5, q6_upper=0.5,
                          rho_lower=2.7, rho_upper=2.7)
        d = diagnose_block(blk)
        self.assertEqual(d['verdict'], 'solid')

    def test_both_liquid(self):
        blk = _make_block(q4_lower=0.02, q4_upper=0.01,
                          q6_lower=0.01, q6_upper=0.01,
                          rho_lower=2.5, rho_upper=2.5)
        d = diagnose_block(blk)
        self.assertEqual(d['verdict'], 'liquid')

    def test_coexist(self):
        # Lower solid, upper liquid, with stratification in q4 and rho.
        blk = _make_block(q4_lower=0.18, q4_upper=0.02,
                          q6_lower=0.5, q6_upper=0.01,
                          rho_lower=2.7, rho_upper=2.5)
        d = diagnose_block(blk)
        self.assertEqual(d['verdict'], 'coexist')
        self.assertGreater(d['delta_rho'], DELTA_RHO)

    def test_unknown_when_ambiguous(self):
        # Lower solid, upper also above Q4_SOLID but q6 too low → not solid,
        # not liquid, not coexist.
        blk = _make_block(q4_lower=0.18, q4_upper=0.12,
                          q6_lower=0.5, q6_upper=0.05,
                          rho_lower=2.7, rho_upper=2.65)
        d = diagnose_block(blk)
        self.assertEqual(d['verdict'], 'unknown')

    def test_empty_block(self):
        d = diagnose_block({'timestep': 0, 'chunks': []})
        self.assertEqual(d['verdict'], 'unknown')

    def test_n_half_adaptive(self):
        # 6 chunks, n_half defaults to 3
        blk = {'timestep': 0, 'chunks': [
            {'q4': 0.2, 'q6': 0.5, 'density': 2.7}] * 3 + [
            {'q4': 0.02, 'q6': 0.01, 'density': 2.5}] * 3}
        d = diagnose_block(blk)
        self.assertEqual(d['verdict'], 'coexist')


class TestDiagnoseRdf(unittest.TestCase):

    def test_solid_multiple_peaks(self):
        import numpy as np
        # Need >= 3 points with g>1.5 in r∈[2.5, 6.0]
        r = np.array([1.0, 2.0, 2.7, 3.5, 4.3, 5.0, 5.8])
        g = np.array([0.0, 3.0, 1.8, 2.0, 1.7, 1.9, 0.3])
        self.assertEqual(diagnose_rdf(r, g), 'solid')

    def test_liquid_single_peak(self):
        import numpy as np
        r = np.array([1.0, 2.0, 2.7, 3.5, 4.3, 5.0, 5.8])
        g = np.array([0.0, 3.0, 0.5, 0.8, 0.4, 0.6, 0.3])
        self.assertEqual(diagnose_rdf(r, g), 'liquid')

    def test_none_input(self):
        self.assertEqual(diagnose_rdf(None, None), 'unknown')

    def test_custom_thresholds(self):
        import numpy as np
        r = np.array([3.0, 4.0, 5.0])
        g = np.array([1.6, 1.7, 1.8])
        # With default min_peaks=3 → solid; with min_peaks=5 → liquid
        self.assertEqual(diagnose_rdf(r, g, min_peaks=3), 'solid')
        self.assertEqual(diagnose_rdf(r, g, min_peaks=5), 'liquid')


class TestCrossValidateRdf(unittest.TestCase):
    """C3: RDF must NOT mutate the chunk verdict, only annotate confidence."""

    def test_coexist_plus_solid_rdf_low_confidence(self):
        result = cross_validate_with_rdf('coexist', 'solid')
        self.assertEqual(result['confidence'], 'low')
        self.assertEqual(result['rdf_verdict'], 'solid')

    def test_coexist_plus_liquid_rdf_high_confidence(self):
        result = cross_validate_with_rdf('coexist', 'liquid')
        self.assertEqual(result['confidence'], 'high')

    def test_solid_chunk_medium_confidence(self):
        result = cross_validate_with_rdf('solid', 'liquid')
        self.assertEqual(result['confidence'], 'medium')

    def test_liquid_chunk_medium_confidence(self):
        result = cross_validate_with_rdf('liquid', 'solid')
        self.assertEqual(result['confidence'], 'medium')


class TestDiagnoseCase(unittest.TestCase):
    """End-to-end: synthetic chunk.profile + optional RDF."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmpdir)

    def test_coexist_detection(self):
        path = os.path.join(self.tmpdir, 'chunk.profile')
        # 20 blocks, all coexist
        blocks = [_make_block(0.18, 0.02, q6_lower=0.5, q6_upper=0.01,
                              rho_lower=2.7, rho_upper=2.5)
                  for _ in range(20)]
        for i, b in enumerate(blocks):
            b['timestep'] = 1000 * (i + 1)
        _write_chunk_file(path, blocks)

        result = diagnose_case(path)
        self.assertEqual(result['verdict'], 'coexist')
        self.assertGreater(result['n_blocks_voted'], 0)
        self.assertIsNone(result['confidence'])  # no RDF given

    def test_solid_detection(self):
        path = os.path.join(self.tmpdir, 'chunk.profile')
        blocks = [_make_block(0.18, 0.19, q6_lower=0.5, q6_upper=0.5,
                              rho_lower=2.7, rho_upper=2.7)
                  for _ in range(20)]
        for i, b in enumerate(blocks):
            b['timestep'] = 1000 * (i + 1)
        _write_chunk_file(path, blocks)
        result = diagnose_case(path)
        self.assertEqual(result['verdict'], 'solid')

    def test_vote_window_adaptive_short_sequence(self):
        """C5: with only 4 blocks, adaptive_window = min(5, 1) = 1."""
        path = os.path.join(self.tmpdir, 'chunk.profile')
        blocks = [_make_block(0.18, 0.02, q6_lower=0.5, q6_upper=0.01,
                              rho_lower=2.7, rho_upper=2.5)
                  for _ in range(4)]
        for i, b in enumerate(blocks):
            b['timestep'] = 1000 * (i + 1)
        _write_chunk_file(path, blocks)
        result = diagnose_case(path, vote_window=5)
        self.assertEqual(result['n_blocks_voted'], 1)
        self.assertEqual(result['verdict'], 'coexist')

    def test_missing_chunk_file(self):
        result = diagnose_case('/nonexistent/chunk.profile')
        self.assertEqual(result['verdict'], 'unknown')
        self.assertEqual(result['n_blocks_voted'], 0)

    def test_rdf_annotation(self):
        """RDF adds confidence without changing verdict."""
        chunk_path = os.path.join(self.tmpdir, 'chunk.profile')
        rdf_path = os.path.join(self.tmpdir, 'rdf_top.txt')

        blocks = [_make_block(0.18, 0.02, q6_lower=0.5, q6_upper=0.01,
                              rho_lower=2.7, rho_upper=2.5)
                  for _ in range(20)]
        for i, b in enumerate(blocks):
            b['timestep'] = 1000 * (i + 1)
        _write_chunk_file(chunk_path, blocks)

        # Write a minimal RDF file (LAMMPS fix ave/time format)
        with open(rdf_path, 'w') as f:
            f.write('# Time-averaged data for fix RDF\n')
            f.write('100000 200\n')
            for i in range(200):
                f.write(f'{100001 + i} {0.05 * (i + 1)} {3.0 * (1.0 if i < 5 else 0.5)}\n')

        result = diagnose_case(chunk_path, rdf_top_path=rdf_path)
        # verdict stays coexist; confidence is annotated
        self.assertEqual(result['verdict'], 'coexist')
        self.assertIsNotNone(result['rdf_top_verdict'])
        self.assertIn(result['confidence'], ('high', 'low', 'medium'))


if __name__ == '__main__':
    unittest.main()
