import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import threading
from concurrent.futures import ThreadPoolExecutor
from unittest.mock import patch
import unittest
import pysam
from pantrans.annotation_runtime import FastaHandles, ordered_results

ROOT=Path(__file__).resolve().parents[1]
from pantrans import cds_construct as p

class AnnotationOptimizationTests(unittest.TestCase):
    def test_thread_private_handles_reused_and_closed(self):
        with tempfile.TemporaryDirectory() as td:
            fasta=Path(td)/'ref.fa';fasta.write_text('>g\nACGT\n');pysam.faidx(str(fasta))
            barrier=threading.Barrier(2)
            with FastaHandles() as handles:
                def read(i):
                    a=handles.get(fasta);barrier.wait()
                    self.assertIs(a,handles.get(fasta));self.assertEqual(a.fetch('g'),'ACGT')
                    return a
                with ThreadPoolExecutor(2) as pool:result=list(pool.map(read,range(2)))
                self.assertIsNot(result[0],result[1])
            self.assertTrue(all(h.closed for h in result))

    def test_ordered_bounded_results_and_errors(self):
        with ThreadPoolExecutor(2) as pool:
            self.assertEqual(list(ordered_results(pool,lambda x:x*x,range(10),3)),[x*x for x in range(10)])
            with self.assertRaises(ZeroDivisionError):
                list(ordered_results(pool,lambda x:1/x,[1,0,2],2))

    def test_default_last_only_and_optional_pre(self):
        for extra,expected in [([],['last']),(['--annotate-pre'],['pre','last'])]:
            with tempfile.TemporaryDirectory() as td:
                out=Path(td)/'out'
                with patch.object(sys,'argv',['pipeline','--prepared',td,'--output',str(out)]+extra),patch.object(p,'load_bundle',return_value=({}, {}, {})),patch.object(p,'graph_stage',return_value=([['a']],[['b']])),patch.object(p,'annotation_stage',return_value={}) as annotate,patch('pantrans.cds_package.write_package'):
                    p.main()
                self.assertEqual([c.args[0] for c in annotate.call_args_list],expected)

    def test_reuse_rejects_incompatible_partition(self):
        with tempfile.TemporaryDirectory() as td:
            source=Path(td)
            (source/'pre.cluster').write_text('a\ta\n')
            with self.assertRaisesRegex(ValueError,'partition'):
                p.reuse_graph(source,source,{'a':{}})

    def test_reuse_main_skips_mapping_and_preserves_inputs(self):
        with tempfile.TemporaryDirectory() as td:
            source=Path(td)/'old';source.mkdir();out=Path(td)/'new'
            for name in ('pre.cluster','last.cluster'):(source/name).write_text('a\tb\n')
            for name in ('cds.raw.bam','cds.filtered.bam','graph.edges.tsv','graph.summary.json'):
                (source/name).write_text('fixture')
            before={f.name:f.read_bytes() for f in source.iterdir()}
            with patch.object(sys,'argv',['pipeline','--prepared',td,'--reuse-graph',str(source),'--output',str(out)]),patch.object(p,'load_bundle',return_value=({'a':{},'b':{}},{},{})),patch.object(p,'graph_stage') as graph,patch.object(p,'annotation_stage',return_value={}) as annotate,patch('pantrans.cds_package.write_package'):
                p.main()
            graph.assert_not_called()
            self.assertEqual(annotate.call_args.args[:2],('last',[['a','b']]))
            self.assertEqual(before,{f.name:f.read_bytes() for f in source.iterdir()})
            self.assertTrue((out/'reused_graph.json').is_file())
