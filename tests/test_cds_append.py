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
from unittest.mock import patch

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
        out=self.root/name
        construct.main(['--prepared',str(prepared),'-o',str(out),'-t','2','--workers','2'])
        return out

    def extend(self,history,hp,query,name,legacy=False):
        out=self.root/name
        append.run(argparse.Namespace(history_results=history,history_prepared=hp,
            query_prepared=query,output=out,threads=2,workers=2,import_legacy_history=legacy))
        return out

    def compare(self,left,right):
        for name in ('graph.edges.tsv','pre.cluster','last.cluster','last.gff3','last.cdna.fa','last.gdna.fa','last.complete.cdna.fa','last.model_sources.tsv','last.alignment_status.tsv','last.summary.json'):
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
        with self.assertRaisesRegex(ValueError,'does not import'):
            self.extend(history,self.a,self.b,'legacy',legacy=True)
        self.assertFalse((self.root/'legacy').exists())

    def test_duplicate_variety_and_input_tampering(self):
        with self.assertRaisesRegex(ValueError,'varieties overlap'):check_inputs(self.a,self.a)
        history=self.build(self.a,'history')
        (self.a/'sample_mrna.fa').write_text('changed')
        with self.assertRaisesRegex(ValueError,'inputs changed'):load_package(history,self.a,longest=True)

    def test_no_valid_new_cds_keeps_gene_and_annotation(self):
        no_cds=self.prepare('N','N_gene','AAA'+self.cds[3:])
        history=self.build(self.a,'history')
        result=self.extend(history,self.a,no_cds,'no_cds_append')
        direct=self.build(result/'merged_prepared','direct')
        self.compare(result,direct)
        summary=json.loads((result/'summary.json').read_text())
        self.assertTrue(summary['cross_alignments']['A']['skipped_empty_cds'])
        self.assertFalse(summary['cross_alignments']['N']['skipped_empty_cds'])
        self.assertIn('N_gene',(result/'last.gff3').read_text())

    def test_graph_edge_tampering_rejected(self):
        history=self.build(self.a,'history')
        with (history/'graph.edges.tsv').open('a') as f:f.write('A_gene\tmissing\n')
        with self.assertRaisesRegex(ValueError,'edges changed'):load_package(history,self.a,longest=True)

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

    def test_multiple_old_and_new_varieties_map_only_missing_pairs(self):
        from pantrans.cds_package import digest
        from pantrans.prepared_bundle import prepared_file
        old=self.root/'old_prepared';merge_prepared(self.a,self.b,old)
        d=self.prepare('D','D_gene',self.cds)
        new=self.root/'new_prepared';merge_prepared(self.c,d,new)
        history=self.build(old,'history_multi')
        before=digest(history/'graph.edges.tsv')
        extended=self.extend(history,old,new,'extended_multi')
        mapping=json.loads((extended/'graph.mapping.json').read_text())
        self.assertEqual(mapping['reused_variety_pairs'],4)
        self.assertEqual(mapping['ordered_variety_pairs'],12)
        self.assertEqual(mapping['total_variety_pairs'],16)
        for task in mapping['targets']:
            directory=extended/task['directory']
            command=json.loads((directory/'command.json').read_text())['minimap2']
            self.assertEqual(command[command.index('-f')+1],'20000')
            self.assertEqual(command[command.index('-N')+1],'2000')
            with pysam.AlignmentFile(directory/'raw.bam') as bam:
                records=[r for r in bam if not r.is_secondary and not r.is_supplementary]
                queries={r.query_name for r in records}
                self.assertEqual(len(bam.references),1)
                self.assertEqual(queries,{'Z_gene.1','D_gene.1'} if task['variety'] in ('A','CS')
                                 else {'A_gene.1','CS_gene.1','Z_gene.1','D_gene.1'})
            summary=json.loads((directory/'summary.json').read_text())
            self.assertEqual(summary['counts']['edges'],2 if task['query_scope']=='new' else 4)
        package,_=load_package(extended,extended/'merged_prepared',longest=True)
        self.assertEqual(package['mapping'],load_package(history,old,longest=True)[0]['mapping'])
        self.assertEqual(before,digest(history/'graph.edges.tsv'))
        direct=self.build(extended/'merged_prepared','direct_multi')
        self.compare(extended,direct)
        restarted=self.root/'restarted_append'
        construct.main(['--prepared',str(extended/'merged_prepared'),'--reuse-graph',str(extended),
                        '-o',str(restarted),'-t','2','--workers','2'])
        self.assertEqual((extended/'last.gff3').read_bytes(),(restarted/'last.gff3').read_bytes())

    def test_n100_history_rejected_without_output(self):
        from pantrans.cds_package import mapping_policy
        history=self.build(self.a,'old_policy')
        path=history/'cds.graph.json';data=json.loads(path.read_text())
        data['mapping']=mapping_policy();path.write_text(json.dumps(data))
        with self.assertRaisesRegex(ValueError,'policy differs'):
            self.extend(history,self.a,self.b,'rejected')
        self.assertFalse((self.root/'rejected').exists())

    def test_mapping_failure_blocks_final_package_and_annotation(self):
        history=self.build(self.a,'history_failed')
        with patch('pantrans.cds_target_graph.minimap2_map',side_effect=RuntimeError('mapping failed')):
            with self.assertRaisesRegex(RuntimeError,'mapping failed'):
                self.extend(history,self.a,self.b,'failed')
        failed=self.root/'failed'
        self.assertEqual(json.loads((failed/'graph.mapping.json').read_text())['status'],'failed')
        for name in ['summary.json','cds.graph.json','last.gff3']:
            self.assertFalse((failed/name).exists())

    def test_no_transcript_new_gene_is_retained(self):
        self.prepare('N','N_gene','AAA'+self.cds[3:])
        annotation=self.root/'N.gff3'
        annotation.write_text('Chr1A\tx\tgene\t1\t20\t.\t+\t.\tID=N_gene\n')
        gene_only=self.root/'gene_only'
        preprocess([self.root/'N.fa'],[annotation],['N'],gene_only)
        history=self.build(self.a,'gene_only_history')
        result=self.extend(history,self.a,gene_only,'gene_only_append')
        self.assertIn('N_gene',(result/'last.cluster').read_text())
        self.assertIn('N_gene',(result/'last.unannotated_representatives.tsv').read_text())
        self.assertEqual(json.loads((result/'validation.json').read_text())['unannotated_representatives'],1)

if __name__=='__main__':unittest.main()
