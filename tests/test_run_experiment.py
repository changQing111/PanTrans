"""Runner integration checks with tiny fake inputs; no Slurm jobs submitted."""
import importlib.util
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

spec = importlib.util.spec_from_file_location('runner', Path(__file__).resolve().parents[1] / 'scripts/run_experiment.py')
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)


class RunnerTest(unittest.TestCase):
    def make_run(self, root):
        run = Path(root) / 'example'; run.mkdir(); (run / 'logs').mkdir()
        (run / 'snapshot').mkdir(); (run / 'snapshot/code.py').write_text('x=1\n')
        source = Path(root) / 'input.fa'; source.write_text('>g\nACGT\n')
        d = dict(run_id='example', status='prepared', inputs={'input': runner.identity(source)},
                 snapshot_files={'code.py': runner.digest(run / 'snapshot/code.py')}, argv=['construct'],
                 output_dir=str(run / 'results'))
        runner.save(run / 'run.json', d)
        return run

    def test_execution_and_index_capture(self):
        with tempfile.TemporaryDirectory() as root:
            run = self.make_run(root)
            def command(argv, **kwargs):
                self.assertEqual(kwargs['env']['PYTHONPATH'], str(run / 'snapshot/src'))
                (run / 'logs/time.txt').write_text('Maximum resident set size (kbytes): 1234\n')
                return 0
            with patch.object(runner, 'capture', return_value={'exit_code': 0, 'output': 'test'}), patch.object(runner.subprocess, 'call', side_effect=command):
                self.assertEqual(runner.execute(run), 0)
            d = json.loads((run / 'run.json').read_text())
            self.assertEqual(d['status'], 'succeeded')
            self.assertEqual(d['gnu_time_max_rss_kib'], 1234)
            self.assertIn('sha256', d['inputs']['input'])
            self.assertIn('succeeded', (Path(root) / 'runs.tsv').read_text())

    def test_changed_snapshot_fails_before_analysis(self):
        with tempfile.TemporaryDirectory() as root:
            run = self.make_run(root)
            (run / 'snapshot/code.py').write_text('x=2\n')
            with patch.object(runner.subprocess, 'call') as call:
                self.assertEqual(runner.execute(run), 1)
                call.assert_not_called()
            self.assertIn('Code snapshot changed', json.loads((run / 'run.json').read_text())['error'])

    def test_oom_accounting_finalizes_unfinished_run(self):
        with tempfile.TemporaryDirectory() as root:
            run = self.make_run(root)
            d = json.loads((run / 'run.json').read_text()); d.update(job_id='42', status='running')
            runner.save(run / 'run.json', d)
            result = {'exit_code': 0, 'output': 'JobID|State|MaxRSS|ExitCode\n42|OUT_OF_MEMORY||0:9\n42.batch|OUT_OF_MEMORY|1000K|0:9\n'}
            with patch.object(runner, 'capture', return_value=result): runner.collect(run)
            d = json.loads((run / 'run.json').read_text())
            self.assertEqual(d['status'], 'out_of_memory')
            self.assertEqual(len(d['slurm_accounting']), 2)

    def test_submit_holds_until_collector_recorded(self):
        with tempfile.TemporaryDirectory() as root:
            run = self.make_run(root)
            with patch.object(runner.subprocess, 'check_output', side_effect=['42\n', '43\n']) as submit, patch.object(runner.subprocess, 'run') as release:
                runner.submit(run)
                self.assertIn('--hold', submit.call_args_list[0].args[0])
                self.assertIn('--dependency=afterany:42', submit.call_args_list[1].args[0])
                release.assert_called_once_with(['scontrol', 'release', '42'], check=True)
            self.assertEqual(json.loads((run / 'run.json').read_text())['collector_job_id'], '43')


if __name__ == '__main__':
    unittest.main()
