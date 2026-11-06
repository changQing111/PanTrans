import csv
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pysam

from pantrans.exon_filter import ExonPolicy, LOGIC_ID
from pantrans.align_filter import filter_bam, minimap2_map
from pantrans.graph_package import validate_edge_provenance
from pantrans import pipeline


class ExonFilterTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.gene = self.root / "genes.tsv"
        self.exon = self.root / "exons.tsv"
        self.gene.write_text("gene\tgdna_length\tbody_start0\tbody_end0\nT\t1000\t200\t800\n")
        self.exon.write_text("gene\ttranscript\tgdna_start0\tgdna_end0\tgdna_strand\n"
                             "T\tT.1\t200\t300\t+\nT\tT.1\t700\t800\t+\n")
        self.policy = ExonPolicy([self.gene], [self.exon])
        self.header = pysam.AlignmentHeader.from_dict({"SQ": [{"SN": "T", "LN": 1000}]})

    def read(self, name="Q.1", start=200, cigar="100M", flag=256, nm=0):
        r = pysam.AlignedSegment(self.header)
        r.query_name, r.reference_id, r.reference_start = name, 0, start
        r.cigarstring, r.flag = cigar, flag
        if nm is not None:
            r.set_tag("NM", nm)
        return r

    def test_side_flank_intron_and_reverse_rejected(self):
        self.assertTrue(self.policy.accepts(self.read(), 100))
        self.assertFalse(self.policy.accepts(self.read(start=10), 100))
        self.assertFalse(self.policy.accepts(self.read(start=400), 100))
        self.assertFalse(self.policy.accepts(self.read(flag=272), 100))

    def test_splice_gaps_excluded_from_overlap(self):
        self.assertTrue(self.policy.accepts(self.read(cigar="100M400N100M"), 200))
        self.assertFalse(self.policy.accepts(self.read(cigar="100M400N100M"), 300))

    def test_same_isoform_required(self):
        self.exon.write_text("gene\ttranscript\tgdna_start0\tgdna_end0\tgdna_strand\n"
                             "T\tT.1\t200\t300\t+\nT\tT.2\t700\t800\t+\n")
        p = ExonPolicy([self.gene], [self.exon])
        self.assertFalse(p.accepts(self.read(cigar="100M400N100M"), 200))

    def test_boundary_threshold(self):
        self.assertTrue(self.policy.accepts(self.read(start=190), 100))
        self.assertFalse(self.policy.accepts(self.read(start=189), 100))

    def test_indels_do_not_count_as_exon_matches(self):
        self.assertFalse(self.policy.accepts(self.read(cigar="30M70I"), 100))
        self.assertTrue(self.policy.accepts(self.read(cigar="40M10D50M"), 90))

    def test_mask_keeps_coordinates_lengths_and_intron(self):
        src, dst = self.root / "src.fa", self.root / "masked.fa"
        src.write_text(">T description\n" + "A" * 1000 + "\n")
        self.policy.mask(src, dst)
        import Bio.SeqIO
        with dst.open() as handle:
            record = next(Bio.SeqIO.parse(handle, "fasta"))
        self.assertEqual(record.id, "T")
        self.assertEqual(str(record.seq), "N" * 200 + "A" * 600 + "N" * 200)

    def test_filter_secondary_real_length_self_and_edges(self):
        source, output = self.root / "raw.bam", self.root / "filtered.bam"
        records = [self.read(), self.read(name="T.1", start=0),
                   self.read(name="shortcoverage.1", cigar="100M100S"),
                   self.read(name="badpid.1", nm=20), self.read(name="noNM.1", nm=None)]
        with pysam.AlignmentFile(source, "wb", header=self.header) as bam:
            for record in records:
                bam.write(record)
        edges, _ = filter_bam(source, output, {"T": ["c", 0, 600]}, exon_policy=self.policy,
                              query_lengths={"Q.1": 100, "T.1": 100, "shortcoverage.1": 200,
                                             "badpid.1": 100, "noNM.1": 100})
        self.assertEqual(edges, [("Q", "T")])
        with pysam.AlignmentFile(output) as bam:
            self.assertEqual(len(list(bam)), 1)
            self.assertIn("pantrans-filtered=" + LOGIC_ID, bam.header.to_dict()["CO"])
        self.assertEqual(json.loads(Path(str(output) + ".filter.json").read_text())["counts"]["retained_records"], 1)

    def test_invalid_tables_targets_and_lengths(self):
        with self.assertRaises(ValueError):
            self.policy.validate_targets({"T": 999})
        with self.assertRaises(ValueError):
            self.policy.validate_targets({"unknown": 1000})
        with self.assertRaises(ValueError):
            ExonPolicy([self.gene], [self.exon], overlap=0)

    def test_annotation_and_mode_provenance(self):
        p = dict(self.policy.provenance(), filter_thresholds_assumed=False)
        validate_edge_provenance({"edge_generations": [p]}, p)
        for changed in [dict(p, target_mask=None), dict(p, exon_thresholds={}),
                        dict(p, gene_region_digests={"T": "changed"}), dict(p, filter_logic_id="legacy")]:
            with self.assertRaisesRegex(ValueError, "incompatible"):
                validate_edge_provenance({"edge_generations": [p]}, changed)

    def test_pipeline_wrappers_forward_policy_for_rescue(self):
        with patch.object(pipeline, "minimap2_map_rescue") as mapper:
            pipeline._map_exons("c", "g", 2, "b", self.policy, rescue=True)
            self.assertIs(mapper.call_args.kwargs["exon_policy"], self.policy)
        with patch.object(pipeline, "filter_bam") as filt, patch.object(pipeline, "get_fasta_len", return_value={"Q.1": 100}):
            pipeline._filter_exons("b", "o", {}, self.policy, "c", collect_edges=False)
            self.assertEqual(filt.call_args.kwargs["query_lengths"], {"Q.1": 100})


if __name__ == "__main__":
    unittest.main()
