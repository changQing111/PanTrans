"""Official CDS dispatch and real command-line integration."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from pantrans.main import read_parameters,main
from test_cds_append import CDSAppendTests

ROOT=Path(__file__).resolve().parents[1]

class CDSParserTests(unittest.TestCase):
    def test_construct_defaults(self):
        a=read_parameters(['construct','--prepared','prepared','-o','results'])
        self.assertEqual(a.prepared,Path('prepared'))
        self.assertEqual((a.threads,a.workers),(32,8))
        self.assertFalse(a.annotate_pre)
        self.assertFalse(a.no_evidence_cache)
        self.assertFalse(hasattr(a,'cdna'))

    def test_append_contract(self):
        a=read_parameters(['append','--history-results','old','--history-prepared','prepared',
                          '--query-prepared','new','-o','out'])
        self.assertEqual(a.history_results,Path('old'))
        self.assertEqual(a.query_prepared,Path('new'))

    def test_rejects_invalid_and_mixed_modes(self):
        base=['construct','--prepared','p','-o','o']
        for extra in (['-t','1'],['--workers','0'],['--reuse-graph','r','--genes','g'],
                      ['--no-evidence-cache','--history-evidence-cache','c'],['--cdna','c']):
            with self.subTest(extra=extra),self.assertRaises(SystemExit):read_parameters(base+extra)
        with self.assertRaises(SystemExit):read_parameters(['append','--prepared','p','-o','o'])
        with self.assertRaises(SystemExit):read_parameters(['construct','-c','c','-g','g','-b','b','-r','CS','-p','P','-o','out'])

    def test_dispatch_and_failure_status(self):
        with patch.object(sys,'argv',['pantrans','construct','--prepared','p','-o','o']),patch('pantrans.cds_construct.run') as run:
            self.assertEqual(main(),0);run.assert_called_once()
        with patch.object(sys,'argv',['pantrans','construct','--prepared','p','-o','o']),patch('pantrans.cds_construct.run',side_effect=ValueError('bad prepared')):
            self.assertEqual(main(),1)


@unittest.skipUnless(shutil.which('minimap2') and shutil.which('samtools'),'alignment executables required')
class CDSCLIIntegrationTests(unittest.TestCase):
    setUp=CDSAppendTests.setUp
    prepare=CDSAppendTests.prepare
    compare=CDSAppendTests.compare

    def command(self,*args):
        env=dict(os.environ,PYTHONPATH=str(ROOT/'src'))
        result=subprocess.run([sys.executable,'-m','pantrans',*map(str,args)],cwd=self.root,env=env,text=True,capture_output=True)
        self.assertEqual(result.returncode,0,result.stdout+'\n'+result.stderr)
        return result

    def test_official_construct_append_and_repeated_append(self):
        history=self.root/'history';extended=self.root/'extended';third=self.root/'third';direct=self.root/'direct'
        self.command('construct','--prepared',self.a,'-o',history,'-t',2,'--workers',2)
        self.assertFalse((history/'pre.gff3').exists())
        self.assertTrue((history/'last.evidence.sqlite').exists())
        self.assertFalse(json.loads((history/'summary.json').read_text())['experimental'])
        self.command('append','--history-results',history,'--history-prepared',self.a,
                     '--query-prepared',self.b,'-o',extended,'-t',2,'--workers',2)
        self.command('append','--history-results',extended,'--history-prepared',extended/'merged_prepared',
                     '--query-prepared',self.c,'-o',third,'-t',2,'--workers',2)
        self.command('construct','--prepared',third/'merged_prepared','-o',direct,'-t',2,'--workers',2)
        self.compare(third,direct)
        self.assertGreater(json.loads((third/'last.cache_stats.json').read_text())['hits'],0)

    def test_official_and_compatibility_wrapper_match(self):
        official=self.root/'official';wrapper=self.root/'wrapper'
        self.command('construct','--prepared',self.a,'-o',official,'-t',2,'--workers',1)
        result=subprocess.run([sys.executable,str(ROOT/'scripts/experimental_cds_construct.py'),
            '--prepared',str(self.a),'-o',str(wrapper),'-t','2','--workers','1'],cwd=self.root,
            env=dict(os.environ,PYTHONPATH=str(ROOT/'src')),capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)
        self.compare(official,wrapper)
