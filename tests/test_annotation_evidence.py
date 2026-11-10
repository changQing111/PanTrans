"""Evidence-cache equivalence, negative evidence, and content invalidation."""
import copy
import json
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch

from test_cds_append import CDSAppendTests,construct
from pantrans.annotation_evidence import EvidenceCache

@unittest.skipUnless(shutil.which('minimap2') and shutil.which('samtools'),'alignment executables required')
class EvidenceIntegrationTests(unittest.TestCase):
    setUp=CDSAppendTests.setUp
    prepare=CDSAppendTests.prepare

    def stage(self,prepared,models,tids,clusters,name,sources=None,local_fasta=False):
        out=self.root/name;out.mkdir()
        construct.annotation_stage('last',clusters,prepared,models,tids,out,2,cache_sources=sources,local_fasta=local_fasta)
        return out

    def assert_outputs(self,a,b):
        for suffix in ('gff3','cdna.fa','model_sources.tsv','alignment_status.tsv','summary.json'):
            self.assertEqual((a/('last.'+suffix)).read_bytes(),(b/('last.'+suffix)).read_bytes(),suffix)
        import pysam
        with pysam.AlignmentFile(a/'last.annotation.bam') as x,pysam.AlignmentFile(b/'last.annotation.bam') as y:
            self.assertEqual([r.to_string() for r in x],[r.to_string() for r in y])

    def test_partial_hits_sequence_model_changes_and_negative_evidence(self):
        from pantrans.cds_merge import merge_prepared
        from pantrans.prepared_bundle import load_prepared_bundle,prepared_file
        ab=self.root/'ab';merge_prepared(self.a,self.b,ab)
        allp=self.root/'all';merge_prepared(ab,self.c,allp)
        bundle=self.root/'bundle';bundle.mkdir()
        _,models,tids=load_prepared_bundle(allp,bundle)
        # A pre-existing transfer reused when another cluster member is added.
        cold=self.stage(allp,models,tids,[['CS_gene','A_gene']],'cold',[])
        cluster=[['CS_gene','A_gene','Z_gene']]
        baseline=self.stage(allp,models,tids,cluster,'baseline')
        warm=self.stage(allp,models,tids,cluster,'warm',[cold/'last.evidence.sqlite'])
        self.assert_outputs(baseline,warm)
        counts=json.loads((warm/'last.cache_stats.json').read_text())
        self.assertEqual((counts['hits'],counts['misses']),(1,1))
        hot=self.stage(allp,models,tids,cluster,'hot',[warm/'last.evidence.sqlite'])
        self.assert_outputs(baseline,hot)
        self.assertEqual(json.loads((hot/'last.cache_stats.json').read_text())['computed_pairs'],0)
        with patch.object(construct,'evidence_policy',return_value='changed-algorithm-version'):
            invalidated=self.stage(allp,models,tids,cluster,'policy_changed',[warm/'last.evidence.sqlite'])
        self.assert_outputs(baseline,invalidated)
        self.assertEqual(json.loads((invalidated/'last.cache_stats.json').read_text())['misses'],2)
        local=self.stage(allp,models,tids,cluster,'local',[warm/'last.evidence.sqlite'],local_fasta=True)
        self.assert_outputs(baseline,local)
        self.assertTrue((local/'last.local_fasta.json').is_file())
        # Changed CDS coordinate metadata must invalidate old evidence.
        altered=copy.deepcopy(models);altered['A_gene.1']['qcds']=[(41,940)]
        changed=self.stage(allp,altered,tids,cluster,'changed',[warm/'last.evidence.sqlite'])
        fresh=self.stage(allp,altered,tids,cluster,'changed_fresh')
        self.assert_outputs(fresh,changed)
        self.assertEqual(json.loads((changed/'last.cache_stats.json').read_text())['misses'],1)
        self.assertIn('no_accepted_model',(changed/'last.alignment_status.tsv').read_text())
        negative=self.stage(allp,altered,tids,cluster,'negative',[changed/'last.evidence.sqlite'])
        self.assert_outputs(fresh,negative)
        self.assertEqual(json.loads((negative/'last.cache_stats.json').read_text())['computed_pairs'],0)
        # Same IDs and coordinates but changed target sequence must invalidate all pairs.
        path=prepared_file(allp,'gdna');path.write_text(path.read_text().replace('ATG','ACG'))
        changedseq=self.stage(allp,models,tids,cluster,'changedseq',[warm/'last.evidence.sqlite'])
        freshseq=self.stage(allp,models,tids,cluster,'freshseq')
        self.assert_outputs(freshseq,changedseq)
        self.assertEqual(json.loads((changedseq/'last.cache_stats.json').read_text())['misses'],2)

class CacheLifecycleTests(unittest.TestCase):
    def test_failed_stage_never_publishes_cache(self):
        with tempfile.TemporaryDirectory() as td:
            path=Path(td)/'cache.sqlite'
            with self.assertRaises(RuntimeError):
                with EvidenceCache(path):raise RuntimeError('failed annotation')
            self.assertFalse(path.exists())

    def test_missing_parent_never_created(self):
        with tempfile.TemporaryDirectory() as td:
            path=Path(td)/'missing.sqlite'
            with self.assertRaises(Exception):
                with EvidenceCache(Path(td)/'out.sqlite',[path]):pass
            self.assertFalse(path.exists())

class MergeIOTests(unittest.TestCase):
    def test_links_are_opt_in_and_content_matches(self):
        from test_preprocess_cds import PreprocessCDSTests
        from pantrans.prepared_bundle import load_prepared_bundle
        fixture=PreprocessCDSTests();fixture.setUp()
        try:
            fixture.prepare()
            out=fixture.root/'linked';out.mkdir()
            load_prepared_bundle(fixture.root/'prepared',out,graph_fasta_mode='link')
            self.assertTrue((out/'body.fa').is_symlink())
            self.assertEqual((out/'body.fa').read_bytes(),(fixture.root/'prepared/custom_gene.fa').read_bytes())
        finally:fixture.doCleanups()
