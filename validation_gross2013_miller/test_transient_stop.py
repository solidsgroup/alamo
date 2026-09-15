#!/usr/bin/env python3
"""Check that automatic stopping requires a confirmed, adequately sampled rate."""
from datetime import datetime,timedelta,timezone
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
import assess_steady
import analyze_runs
import monitor_transient_sweep as monitor


class StopTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.case=Path(self.tmp.name)
        (self.case/'case.json').write_text(json.dumps(dict(packing=dict(bins={'20':{}}))))

    def rows(self,end,accelerating=False):
        return [dict(time_s=float(t),mean_solid_height_m=float(.001-.02*t-(5*t*t if accelerating else 0)),
                     temperature_max_K=3000.,gas_heat_release_W_m2=1e6)
                for t in np.arange(0,end+1e-10,25e-6)]

    def assess(self,rows):
        with patch.object(assess_steady,'history',return_value=rows):
            return assess_steady.assessment(self.case,completed_only=True)

    def test_requires_particle_sampling_and_later_confirmation(self):
        short=self.assess(self.rows(.00075));self.assertFalse(short['accepted'])
        first=self.assess(self.rows(.002));self.assertEqual(first['status'],'steady_candidate')
        same=self.assess(self.rows(.002));self.assertFalse(same['accepted'])
        accepted=self.assess(self.rows(.0025));self.assertTrue(accepted['accepted'])
        start=datetime(2026,9,15,tzinfo=timezone.utc)
        receipt=dict(status='running',started_utc=start.isoformat())
        self.assertFalse(monitor.record_acceptance(self.case,first,receipt,start))
        self.assertFalse((self.case/'STOP').exists())
        self.assertTrue(monitor.record_acceptance(self.case,accepted,receipt,start+timedelta(seconds=17)))
        self.assertTrue((self.case/'STOP').exists())
        record=json.loads((self.case/'analysis/steady_stop.json').read_text())
        self.assertEqual(record['time_to_steady_wall_s'],17)
        self.assertAlmostEqual(record['time_to_steady_simulation_s'],.0025)
        self.assertAlmostEqual(record['rate_cm_s'],2.)
        self.assertFalse(monitor.record_acceptance(self.case,accepted,receipt,start+timedelta(seconds=20)))

    def test_accelerating_rate_is_not_steady(self):
        self.assertFalse(self.assess(self.rows(.004,accelerating=True))['accepted'])
        self.assertFalse(self.assess(self.rows(.006,accelerating=True))['accepted'])

    def test_incomplete_plot_is_excluded(self):
        output=self.case/'output';output.mkdir()
        for name in ('00000cell','00001cell'):
            p=output/name;p.mkdir();(p/'Header').touch();(p/'Checkpoint').touch()
        (output/'celloutput.visit').write_text('00000cell/Header\n00001cell/Header')
        fake=dict(plotfile='00000cell',time_s=0.,snapshot_schema=3,gas_heat_release_W_m2=1.)
        with patch.object(analyze_runs,'snapshot',return_value=fake) as snapshot:
            rows=analyze_runs.history(self.case,completed_only=True)
        self.assertEqual(len(rows),1)
        snapshot.assert_called_once_with(output/'00000cell')


if __name__=='__main__':unittest.main()
