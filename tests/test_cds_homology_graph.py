import sys,tempfile,unittest
from pathlib import Path
import networkx as nx
import pysam
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from pantrans import cds_construct as p
from pantrans.cds_package import POLICY

class HomologyGraphTests(unittest.TestCase):
    def read(self,cigar,nm=0):
        r=pysam.AlignedSegment();r.reference_start=0;r.cigarstring=cigar;r.set_tag('NM',nm);return r
    def test_intron_retention_example(self):
        r=self.read('948M',2);blocks=[(0,896),(938,948)]
        self.assertIsNone(p.coding_support(r,[(0,948)],blocks))
        hit=p.graph_coding_support(r,948,blocks)
        self.assertAlmostEqual(hit['pid'],946/948)
        self.assertAlmostEqual(hit['query_cds_coverage'],906/948)
        self.assertEqual(hit['target_cds_coverage'],1)
    def test_novel_splice_and_indels_allowed_for_graph_only(self):
        for cigar,qlen,blocks in [('100M30N100M',200,[(0,230)]),('100M1I100M',201,[(0,200)]),('100M1D100M',200,[(0,201)])]:
            r=self.read(cigar,1)
            self.assertIsNotNone(p.graph_coding_support(r,qlen,blocks))
            self.assertIsNone(p.coding_support(r,[(0,qlen)],blocks))
    def test_bam_thresholds_include_exact_90_exclude_below(self):
        with tempfile.TemporaryDirectory() as td:
            raw=Path(td)/'raw.bam';filtered=Path(td)/'filtered.bam'
            with pysam.AlignmentFile(raw,'wb',header={'SQ':[{'SN':'target','LN':200}]}) as bam:
                for tid,cigar,nm in [('pass_pid','100M',10),('fail_pid','100M',11),('pass_cov','90M10S',0),('fail_cov','89M11S',0)]:
                    r=self.read(cigar,nm);r.reference_id=0;r.query_name=tid;bam.write(r)
            models={tid:{'gene':tid,'cds_len':100} for tid in ['pass_pid','fail_pid','pass_cov','fail_cov']}
            models['target.1']={'valid':True,'cds':[(0,100)]}
            g=nx.DiGraph()
            p.filter_cds_graph(raw,filtered,{'target':{'body_start0':0}},models,{'target':['target.1']},g)
            self.assertEqual(set(g.edges()),{('pass_pid','target'),('pass_cov','target')})
    def test_invalid_records_still_rejected(self):
        for flag in (4,16,2048):
            r=self.read('100M');r.flag=flag
            self.assertIsNone(p.graph_coding_support(r,100,[(0,100)]))
        r=self.read('100M');r.set_tag('NM',None)
        self.assertIsNone(p.graph_coding_support(r,100,[(0,100)]))
    def test_versioned_policy(self):
        self.assertEqual(POLICY['filter'],'cds-homology-v2')
        self.assertFalse(POLICY['structure_checks'])
        self.assertEqual([POLICY[k] for k in ('pid','query_cds_coverage','target_cds_coverage')],[.9]*3)
