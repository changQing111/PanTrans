"""Real minimap2 integration and input-integrity tests for CDS append."""
import argparse
import importlib.util
import json
from pathlib import Path
import random
import shutil
import sys
import tempfile
import unittest

import pysam
from pantrans.preprocess import preprocess
from pantrans.cds_package import write_package, load_package
from pantrans.cds_merge import merge_prepared, check_inputs

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from pantrans import cds_construct as construct
from pantrans import cds_append as append


@unittest.skipUnless(shutil.which('minimap2') and shutil.which('samtools'),'alignment executables required')
class CDSAppendTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        rng=random.Random(8723)
        codons=['AAA','AAC','AAG','ACA','ACC','ACG','ACT','AGC','AGG','AGT','ATA','ATC','ATT','CAA','CAC','CAG','CCA','CCC','CCG','CCT','CGA','CGC','CGG','CGT','CTA','CTC','CTG','CTT','GAA','GAC','GAG','GAT','GCA','GCC','GCG','GCT','GGA','GGC','GGG','GGT','GTA','GTC','GTG','GTT','TAC','TAT','TCA','TCC','TCG','TCT','TGC','TGG','TGT','TTA','TTC','TTG','TTT']
        self.cds='ATG'+''.join(rng.choices(codons,k=298))+'TAA'
        self.a=self.prepare('A','A_gene',self.cds)
        self.b=self.prepare('CS','CS_gene',self.cds)
        self.c=self.prepare('Z','Z_gene',self.cds)

    def prepare(self,variety,gene,cds):
        genome=self.root/(variety+'.fa');annotation=self.root/(variety+'.gff3')
        seq='C'*40+cds+'G'*40
        genome.write_text('>Chr1A\n'+seq+'\n')
        annotation.write_text('##gff-version 3\n'+ '\n'.join([
            f'Chr1A\tx\tgene\t1\t{len(seq)}\t.\t+\t.\tID={gene}',
            f'Chr1A\tx\tmRNA\t1\t{len(seq)}\t.\t+\t.\tID={gene}.1;Parent={gene}',
            f'Chr1A\tx\texon\t1\t{len(seq)}\t.\t+\t.\tParent={gene}.1',
            f'Chr1A\tx\tCDS\t41\t{40+len(cds)}\t.\t+\t0\tParent={gene}.1'])+'\n')
        directory=self.root/(variety+'_prepared')
        preprocess([genome],[annotation],[variety],directory,prefix='sample',flank=3000)
        return directory

    def build(self,prepared,name):
        out=self.root/name;out.mkdir()
        r,m,t=construct.load_bundle(prepared,[],out)
        _,last=construct.graph_stage(out,r,m,t,2)
        write_package(out,prepared,kind='construct')
        construct.annotation_stage('last',last,prepared,m,t,out,2,cache_sources=[])
        return out

    def extend(self,history,hp,query,name,legacy=False):
        out=self.root/name
        append.run(argparse.Namespace(history_results=history,history_prepared=hp,
            query_prepared=query,output=out,threads=2,workers=2,import_legacy_history=legacy))
        return out

    def compare(self,left,right):
        for name in ('graph.edges.tsv','pre.cluster','last.cluster','last.gff3','last.cdna.fa','last.model_sources.tsv','last.alignment_status.tsv','last.summary.json'):
            self.assertEqual((left/name).read_bytes(),(right/name).read_bytes(),name)
        with pysam.AlignmentFile(left/'last.annotation.bam') as a,pysam.AlignmentFile(right/'last.annotation.bam') as b:
            self.assertEqual(a.header.to_dict(),b.header.to_dict())
            self.assertEqual([r.to_string() for r in a],[r.to_string() for r in b])

    def test_append_construct_equivalence_hub_change_and_chaining(self):
        history=self.build(self.a,'history')
        result=self.extend(history,self.a,self.b,'append')
        direct=self.build(result/'merged_prepared','direct')
        self.compare(result,direct)
        self.assertEqual((history/'last.cluster').read_text().split()[0],'A_gene')
        self.assertEqual((result/'last.cluster').read_text().split()[0],'CS_gene')
        self.assertFalse((result/'pre.gff3').exists())
        third=self.extend(result,result/'merged_prepared',self.c,'third')
        direct3=self.build(third/'merged_prepared','direct3')
        self.compare(third,direct3)
        self.assertGreater(json.loads((third/'last.cache_stats.json').read_text())['hits'],0)

    def test_explicit_legacy_import_and_sequence_rejection(self):
        history=self.build(self.a,'history');(history/'cds.graph.json').unlink()
        with self.assertRaisesRegex(ValueError,'Missing cds.graph'):
            self.extend(history,self.a,self.b,'missing')
        self.assertFalse((self.root/'missing').exists())
        result=self.extend(history,self.a,self.b,'legacy',legacy=True)
        direct=self.build(result/'merged_prepared','direct')
        self.compare(result,direct)
        (history/'cds.fa').write_text('bad')
        with self.assertRaisesRegex(ValueError,'mismatch'):
            self.extend(history,self.a,self.b,'bad_legacy',legacy=True)

    def test_duplicate_variety_and_input_tampering(self):
        with self.assertRaisesRegex(ValueError,'varieties overlap'):check_inputs(self.a,self.a)
        history=self.build(self.a,'history')
        (self.a/'sample_mrna.fa').write_text('changed')
        with self.assertRaisesRegex(ValueError,'inputs changed'):load_package(history,self.a)

    def test_no_valid_new_cds_keeps_gene_and_annotation(self):
        no_cds=self.prepare('N','N_gene','AAA'+self.cds[3:])
        history=self.build(self.a,'history')
        result=self.extend(history,self.a,no_cds,'no_cds_append')
        direct=self.build(result/'merged_prepared','direct')
        self.compare(result,direct)
        summary=json.loads((result/'summary.json').read_text())
        self.assertTrue(summary['cross_alignments']['new_to_all']['skipped_empty_cds'])
        self.assertIn('N_gene',(result/'last.gff3').read_text())

    def test_graph_edge_tampering_rejected(self):
        history=self.build(self.a,'history')
        with (history/'graph.edges.tsv').open('a') as f:f.write('A_gene\tmissing\n')
        with self.assertRaisesRegex(ValueError,'edges changed'):load_package(history,self.a)

    def test_streamed_merge_hashes_match_disk(self):
        from pantrans.cds_package import digest
        from pantrans.prepared_bundle import prepared_file
        hashes={};out=self.root/'hash_merge'
        merge_prepared(self.a,self.b,out,hashes=hashes)
        for key,value in hashes.items():
            self.assertEqual(value,digest(prepared_file(out,key)),key)

    def test_policy_mismatch_and_gene_collision(self):
        history=self.build(self.a,'history')
        path=history/'cds.graph.json';data=json.loads(path.read_text())
        data['policy']['pid']=0.80;path.write_text(json.dumps(data))
        with self.assertRaisesRegex(ValueError,'policy differs'):
            self.extend(history,self.a,self.b,'bad_policy')
        self.assertFalse((self.root/'bad_policy').exists())
        collision=self.prepare('Collision','A_gene',self.cds)
        with self.assertRaisesRegex(ValueError,'IDs overlap'):
            check_inputs(self.a,collision)

    def test_missing_cdna_record_rejected_before_mapping(self):
        (self.b/'sample_mrna.fa').write_text('')
        with self.assertRaisesRegex(ValueError,'FASTA identifiers'):
            merge_prepared(self.a,self.b,self.root/'merged')

if __name__=='__main__':unittest.main()
