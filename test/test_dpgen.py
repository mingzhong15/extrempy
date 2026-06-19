import unittest
import sys, os, json
from unittest.mock import MagicMock

# Setup shared mocks for heavy deps
from test_helpers import setup_mocks
setup_mocks()
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from extrempy.lazy.dpgen import (_generate_temp_list, _generate_pres_list,
                                  split_temperature_range)


class TestGenerateTempList(unittest.TestCase):
    def test_single_point(self):
        # xmin == xmax returns single element list
        self.assertEqual(_generate_temp_list(300, 300), [300])

    def test_range_increasing(self):
        tlist = _generate_temp_list(300, 2000)
        self.assertGreater(len(tlist), 1)
        self.assertEqual(tlist[0], 300)
        self.assertGreater(tlist[-1], 2000)

    def test_strictly_increasing(self):
        tlist = _generate_temp_list(300, 4000)
        for i in range(len(tlist) - 1):
            self.assertGreater(tlist[i+1], tlist[i],
                               'T list must be strictly increasing')

    def test_last_exceeds_xmax(self):
        tlist = _generate_temp_list(300, 1000)
        self.assertGreaterEqual(tlist[-1], 1000,
                                'last element should >= xmax')

    def test_no_gaps(self):
        tlist = _generate_temp_list(300, 3000)
        for i in range(len(tlist) - 1):
            step = tlist[i+1] - tlist[i]
            self.assertGreaterEqual(step, 50)  # min delta_x

    def test_ti_hcp_explore(self):
        # T_explore for Ti-HCP is (300, 1386)
        tlist = _generate_temp_list(300, 1386)
        self.assertGreater(len(tlist), 3)

    def test_ti_liq_explore(self):
        # T_explore for Ti-LIQ is (1553, 3882)
        tlist = _generate_temp_list(1553, 3882)
        self.assertGreater(len(tlist), 5)


class TestSplitTemperatureRange(unittest.TestCase):
    def test_empty_input(self):
        self.assertEqual(split_temperature_range([], 2000), [])

    def test_single_temperature(self):
        result = split_temperature_range([500], 2000)
        self.assertEqual(result, [[500]])

    def test_split_into_multiple_ranges(self):
        tlist = [300, 500, 700, 3000, 3200, 3400]
        ranges = split_temperature_range(tlist, delta_T=2000)
        self.assertGreaterEqual(len(ranges), 2)

    def test_all_in_one_range(self):
        tlist = [300, 310, 320]
        ranges = split_temperature_range(tlist, delta_T=50)
        self.assertEqual(len(ranges), 1)

    def test_preserves_order(self):
        tlist = [100, 300, 200]
        ranges = split_temperature_range(tlist, 500)
        self.assertEqual(ranges[0], [100, 200, 300])


class TestGeneratePresList(unittest.TestCase):
    def test_single_point(self):
        self.assertEqual(_generate_pres_list(1, 1), [1])

    def test_two_arguments(self):
        plist = _generate_pres_list(1, 100)
        self.assertGreater(len(plist), 1)

    def test_log_scale_large_ratio(self):
        plist = _generate_pres_list(1, 10000)
        self.assertGreater(len(plist), 3)


class TestDpgenJobsLogic(unittest.TestCase):
    """Test _set_model_devi_jobs_from_segments logic using a mock object."""

    def setUp(self):
        # Build segments similar to get_phase_segments('Ti')
        self.segs = [
            {'label': 'Ti-HCP', 'structure': 'hcp',
             'T_core': (300, 1155), 'T_explore': (300, 1386)},
            {'label': 'Ti-BCC', 'structure': 'bcc',
             'T_core': (1155, 1941), 'T_explore': (924, 2329)},
            {'label': 'Ti-LIQ', 'structure': 'bcc',
             'T_core': (1941, 3882), 'T_explore': (1553, 3882)},
        ]
        # Build a minimal object with jparam to act as the method's self
        class FakeGen:
            jparam = {"model_devi_jobs": []}
        self.gen = FakeGen()

    def _run(self, **kw):
        from extrempy.lazy.dpgen import DPGENGenerator
        method = DPGENGenerator._set_model_devi_jobs_from_segments
        defaults = dict(nsteps_per_phase=5,
                        init_steps=[1000, 2000, 4000, 8000, 16000],
                        press_grid=[1, 10, 100, 1000, 10000],
                        trj_freq=20, numb_frame_per_iter_per_PT=5,
                        ensemble='npt')
        defaults.update(kw)
        method(self.gen, self.segs, **defaults)
        return self.gen.jparam['model_devi_jobs']

    def test_job_count_15(self):
        jobs = self._run()
        self.assertEqual(len(jobs), 15)

    def test_sys_idx_sequence(self):
        jobs = self._run()
        self.assertEqual(jobs[0]['sys_idx'], [0])    # HCP
        self.assertEqual(jobs[4]['sys_idx'], [0])    # HCP last
        self.assertEqual(jobs[5]['sys_idx'], [1])    # BCC first
        self.assertEqual(jobs[9]['sys_idx'], [1])    # BCC last
        self.assertEqual(jobs[10]['sys_idx'], [2])   # LIQ first
        self.assertEqual(jobs[14]['sys_idx'], [2])   # LIQ last

    def test_nsteps_restart_per_phase(self):
        jobs = self._run()
        # HCP phase
        self.assertEqual(jobs[0]['nsteps'], 1000)
        self.assertEqual(jobs[4]['nsteps'], 16000)
        # BCC phase (restarts at 1000)
        self.assertEqual(jobs[5]['nsteps'], 1000)
        self.assertEqual(jobs[9]['nsteps'], 16000)
        # LIQ phase (restarts)
        self.assertEqual(jobs[10]['nsteps'], 1000)
        self.assertEqual(jobs[14]['nsteps'], 16000)

    def test_press_grid(self):
        jobs = self._run()
        for j in jobs:
            self.assertEqual(j['press'],
                             [1, 10, 100, 1000, 10000])

    def test_trj_freq(self):
        jobs = self._run(trj_freq=50)
        for j in jobs:
            self.assertEqual(j['trj_freq'], 50)

    def test_fp_task_max(self):
        jobs = self._run()
        # Ti-HCP has the most T points (10 for (300,1386))
        from extrempy.lazy.dpgen import _generate_temp_list
        n_t_hcp = len(_generate_temp_list(300, 1386))
        n_t_bcc = len(_generate_temp_list(924, 2329))
        n_t_liq = len(_generate_temp_list(1553, 3882))
        expected_max = max(n_t_hcp, n_t_bcc, n_t_liq) * 5 * 5
        self.assertEqual(self.gen.jparam['fp_task_max'],
                         expected_max)

    def test_fp_task_min(self):
        jobs = self._run()
        from extrempy.lazy.dpgen import _generate_temp_list
        n_t_hcp = len(_generate_temp_list(300, 1386))
        n_t_bcc = len(_generate_temp_list(924, 2329))
        n_t_liq = len(_generate_temp_list(1553, 3882))
        expected_min = min(n_t_hcp, n_t_bcc, n_t_liq) * 5 * 1
        self.assertEqual(self.gen.jparam['fp_task_min'],
                         expected_min)

    def test_nsteps_custom_sequence(self):
        jobs = self._run(init_steps=[2000, 4000, 8000, 16000, 32000])
        self.assertEqual(jobs[0]['nsteps'], 2000)
        self.assertEqual(jobs[4]['nsteps'], 32000)
        self.assertEqual(jobs[5]['nsteps'], 2000)  # restart

    def test_press_grid_custom(self):
        jobs = self._run(press_grid=[1, 100, 10000])
        for j in jobs:
            self.assertEqual(j['press'], [1, 100, 10000])

    def test_index_sequence(self):
        jobs = self._run()
        for i in range(15):
            self.assertEqual(jobs[i]['_idx'], i)

    def test_ensemble_default(self):
        jobs = self._run()
        for j in jobs:
            self.assertEqual(j['ensemble'], 'npt')

    def test_ensemble_custom(self):
        jobs = self._run(ensemble='nvt')
        for j in jobs:
            self.assertEqual(j['ensemble'], 'nvt')


if __name__ == '__main__':
    unittest.main()
