"""Exercise real per-variety mapping, provenance and complete transcript export."""
import json
from pathlib import Path
import shutil
import subprocess
import unittest
from unittest.mock import patch

from Bio.SeqIO.FastaIO import SimpleFastaParser
import networkx as nx
import pysam

from pantrans import cds_construct as construct
from pantrans.cds_merge import merge_prepared
from pantrans.cds_package import load_package
from pantrans.cds_target_graph import filter_target
from pantrans.preprocess import preprocess
from test_cds_append import CDSAppendTests


@unittest.skipUnless(shutil.which('minimap2') and shutil.which('samtools'), 'alignment tools required')
class TargetGraphTests(unittest.TestCase):
    setUp = CDSAppendTests.setUp
    prepare = CDSAppendTests.prepare

    def run_construct(self, prepared, name, *extra):
        out = self.root/name
        construct.main(['--prepared', str(prepared), '-o', str(out), '-t', '2', '--workers', '2', *map(str,extra)])
        return out

    def test_real_targets_queries_edges_and_all_isoforms_for_annotation(self):
        annotation = self.root/'A.gff3'
        text = annotation.read_text()
        with annotation.open('a') as f:
            f.write('\n'.join(line.replace('A_gene.1', 'A_gene.2') for line in text.splitlines()
                               if not line.startswith('#') and '\tgene\t' not in line)+'\n')
        a = self.root/'A_multi'
        preprocess([self.root/'A.fa'], [annotation], ['A'], a)
        merged = self.root/'all'
        merge_prepared(a, self.b, merged)
        out = self.run_construct(merged, 'new')
        mapping = json.loads((out/'graph.mapping.json').read_text())
        self.assertEqual(mapping['query_count'], 2)
        self.assertEqual(mapping['ordered_variety_pairs'], 4)
        self.assertEqual(mapping['status'], 'completed')
        expected = {'A_gene.1','CS_gene.1'}
        for target in mapping['targets']:
            d = out/target['directory']
            cmd = json.loads((d/'command.json').read_text())['minimap2']
            self.assertEqual(cmd[cmd.index('-N')+1], '2000')
            self.assertEqual(cmd[cmd.index('-f')+1], '20000')
            with pysam.AlignmentFile(d/'raw.bam') as bam:
                self.assertEqual(set(bam.references), {'A_gene'} if target['variety']=='A' else {'CS_gene'})
                self.assertEqual({r.query_name for r in bam}, expected)
        self.assertEqual(set(map(tuple,(line.split() for line in (out/'graph.edges.tsv').read_text().splitlines()))),
                         {(a,b) for a in ['A_gene','CS_gene'] for b in ['A_gene','CS_gene']})
        self.assertIn('A_gene.2', (out/'last.model_sources.tsv').read_text())
        package,_ = load_package(out, merged, longest=True)
        self.assertIn('cds_longest', package['prepared']['files'])
        with self.assertRaisesRegex(ValueError, 'policy differs'):load_package(out,merged)
        audit = json.loads((out/'validation.json').read_text())
        self.assertTrue(audit['all_cdna_reconstructed_and_verified'])
        self.assertEqual(audit['unannotated_representatives'], 0)
        with (out/'last.complete.cdna.fa').open() as f:
            self.assertEqual(len(list(SimpleFastaParser(f))), audit['complete_mRNAs'])
        restarted = self.run_construct(merged, 'restart', '--reuse-graph', out)
        again = self.run_construct(merged, 'restart_again', '--reuse-graph', restarted)
        for dest in [restarted,again]:
            for name in ('last.gff3','last.cdna.fa','last.gdna.fa','last.complete.cdna.fa'):
                self.assertEqual((out/name).read_bytes(), (dest/name).read_bytes())
            self.assertFalse((dest/'targets').exists())

    def test_incomplete_or_old_prepared_rejected_before_output(self):
        path = next(self.a.glob('*_preprocess.json'))
        data = json.loads(path.read_text())
        data['files'].pop('cds_longest')
        path.write_text(json.dumps(data))
        from pantrans.prepared_bundle import manifest
        manifest.cache_clear()
        with self.assertRaisesRegex(ValueError, 'rerun pantrans preprocess'):
            self.run_construct(self.a, 'bad')
        self.assertFalse((self.root/'bad').exists())

    def test_tampered_query_and_table_rejected_before_mapping(self):
        from pantrans.prepared_bundle import prepared_file
        table = prepared_file(self.a, 'cds_longest_selection')
        original = table.read_bytes()
        table.write_text(table.read_text().replace('A_gene.1','wrong.1'))
        with patch('pantrans.cds_target_graph.minimap2_map') as mapper:
            with self.assertRaisesRegex(ValueError, 'disagrees with models'):
                self.run_construct(self.a, 'bad_table')
            mapper.assert_not_called()
        table.write_bytes(original)
        query = prepared_file(self.a, 'cds_longest')
        query.write_text(query.read_text().replace(self.cds[:15], 'C'+self.cds[1:15]))
        with patch('pantrans.cds_target_graph.minimap2_map') as mapper:
            with self.assertRaisesRegex(ValueError, 'sequence mismatch'):
                self.run_construct(self.a, 'bad_sequence')
            mapper.assert_not_called()

    def test_mapping_failure_does_not_publish_success(self):
        with patch('pantrans.cds_target_graph.minimap2_map', side_effect=subprocess.CalledProcessError(1, ['minimap2'])):
            with self.assertRaises(subprocess.CalledProcessError):self.run_construct(self.a, 'failed')
        out = self.root/'failed'
        self.assertEqual(json.loads((out/'graph.mapping.json').read_text())['status'], 'failed')
        self.assertFalse((out/'summary.json').exists())
        self.assertFalse((out/'cds.graph.json').exists())

    def test_subset_and_no_transcript_genes_remain_nodes(self):
        annotation = self.root/'A.gff3'
        with annotation.open('a') as f:
            f.write('Chr1A\tx\tgene\t1\t20\t.\t+\t.\tID=unmodeled\n')
        prepared = self.root/'with_unmodeled'
        preprocess([self.root/'A.fa'], [annotation], ['A'], prepared)
        genes = self.root/'genes.txt'
        genes.write_text('A_gene\nunmodeled\n')
        out = self.run_construct(prepared, 'subset', '--genes', genes)
        self.assertFalse((out/'cds.graph.json').exists())
        self.assertEqual(set((out/'last.cluster').read_text().split()), {'A_gene','unmodeled'})
        audit = json.loads((out/'validation.json').read_text())
        self.assertEqual(audit['unannotated_representatives'],1)
        self.assertIn('unmodeled', (out/'last.unannotated_representatives.tsv').read_text())

    def test_target_manifest_tampering_rejected(self):
        out = self.run_construct(self.a,'attested')
        path = out/'graph.mapping.json'
        path.write_text(path.read_text()+'\n')
        with self.assertRaisesRegex(ValueError, 'manifest changed'):load_package(out,self.a,longest=True)


class TargetAuditTests(unittest.TestCase):
    def test_saturation_counts_raw_secondary_not_filtered_hits(self):
        import tempfile
        with tempfile.TemporaryDirectory() as temp:
            d = Path(temp)
            with pysam.AlignmentFile(d/'raw.bam','wb',header={'SQ':[{'SN':'g','LN':100}]}) as bam:
                for i in range(2001):
                    read=pysam.AlignedSegment();read.query_name='g.1';read.reference_id=0
                    read.reference_start=0;read.cigarstring='100M';read.flag=256 if i else 0
                    read.set_tag('NM', 20 if i else 0);bam.write(read)
            regions={'g':dict(variety='CS',body_start0='0',body_end0='100')}
            models={'g.1':dict(gene='g',valid=True,cds_len=100,cds=[(0,100)])}
            graph=nx.DiGraph();graph.add_node('g')
            result,sat=filter_target(d,'CS',{'g.1':'g'},regions,models,{'g':['g.1']},graph)
            self.assertEqual(sat, {'g'})
            self.assertEqual(result['counts']['max_secondary'],2000)
            self.assertEqual(result['counts']['retained'],1)
            self.assertEqual(result['counts']['edges'],1)


if __name__ == '__main__':unittest.main()
