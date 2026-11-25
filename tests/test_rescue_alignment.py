import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import Bio.SeqIO
import pysam

from pantrans import pipeline
from pantrans.gff3 import make_row, write_gff3


class PooledRescueTests(unittest.TestCase):
    def test_mixed_clusters_share_alignment_and_reject_cross_cluster_models(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            cdna, gdna, bed = (root / name for name in ('c.fa', 'g.fa', 'genes.bed'))
            genes = ['existing', 'pre', 'single', 'unresolved', 'hub', 'member']
            cdna.write_text(''.join(f'>{g}.1\n' + 'A' * 100 + '\n' for g in genes))
            gdna.write_text(''.join(f'>{g}\n' + 'A' * 200 + '\n' for g in genes))
            pysam.faidx(str(gdna))
            bed.write_text(''.join(f'Chr1\t0\t200\t{g}\tVariety\t+\n' for g in genes))
            final, pre = root / 'final.gff3', root / 'pre.gff3'

            def model(g):
                return [make_row(g, 'gene', 1, 100, '+', g),
                        make_row(g, 'mRNA', 1, 100, '+', g + '.1', g),
                        make_row(g, 'exon', 1, 100, '+', g + '.1.exon1', g + '.1')]

            write_gff3(final, model('existing'))
            write_gff3(pre, model('pre'))
            observed = {}

            def align(query, target, threads, destination):
                with open(query) as handle:
                    observed['queries'] = {r.id for r in Bio.SeqIO.parse(handle, 'fasta')}
                with open(target) as handle:
                    targets = [r.id for r in Bio.SeqIO.parse(handle, 'fasta')]
                observed['targets'] = set(targets)
                self.assertEqual(threads, 4)
                header = {'HD': {'VN': '1.6'}, 'SQ': [{'SN': g, 'LN': 200} for g in targets]}
                with pysam.AlignmentFile(destination, 'wb', header=header) as bam:
                    for q, ref in [('single.1', 'single'), ('member.1', 'hub'),
                                   ('single.1', 'unresolved')]:
                        read = pysam.AlignedSegment(bam.header)
                        read.query_name = q
                        read.query_sequence = 'A' * 100
                        read.flag = 0
                        read.reference_id = targets.index(ref)
                        read.reference_start = 0
                        read.mapping_quality = 60
                        read.cigarstring = '100M'
                        read.set_tag('NM', 0)
                        bam.write(read)

            clusters = {g: [g] for g in genes if g != 'member'}
            clusters['hub'].append('member')
            with patch.object(pipeline, 'minimap2_map_rescue', side_effect=align) as mapper:
                count = pipeline._rescue_missing_cluster_genes(
                    clusters, str(final), str(cdna), str(gdna), str(bed),
                    None, directory, 'test', '', 4, pre_gff3_path=str(pre))
            mapper.assert_called_once()
            self.assertEqual(observed['targets'], {'single', 'unresolved', 'hub', 'member'})
            self.assertEqual(observed['queries'], {g + '.1' for g in observed['targets']})
            self.assertEqual(count, 3)
            self.assertEqual(pipeline._get_gff3_hub_genes(final), {'existing', 'pre', 'single', 'hub'})
            self.assertTrue((root / 'test_rescue_cdna_align_gdna.filtered.bam').exists())
            self.assertFalse((root / 'test_rescue_work').exists())

    def test_no_missing_genes_does_not_align(self):
        with tempfile.TemporaryDirectory() as directory:
            final = Path(directory) / 'final.gff3'
            write_gff3(final, [])
            with patch.object(pipeline, 'minimap2_map_rescue') as mapper:
                self.assertEqual(pipeline._rescue_missing_cluster_genes(
                    {}, str(final), None, None, None, None, directory, 'test', '', 4), 0)
            mapper.assert_not_called()


if __name__ == '__main__':
    unittest.main()
