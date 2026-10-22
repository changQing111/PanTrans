import tempfile
import unittest
from pathlib import Path

import pysam
from Bio import SeqIO
from Bio.Seq import Seq

from pantrans.transcript_processor import transcript_dedup, get_cdna_from_gff3
from pantrans.gff3 import load_gff3_transcript_models, write_gff3, make_row


class TranscriptStrandTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def run_model(self, dna, bed_strand, alignments, seed=None):
        fasta = self.root / "gene.fa"
        fasta.write_text(">g\n" + dna + "\n")
        bam = self.root / "reads.bam"
        with pysam.AlignmentFile(bam, "wb", header={"HD": {"VN": "1.6"}, "SQ": [{"SN": "g", "LN": len(dna)}]}) as out:
            for tid, flag in alignments:
                read = pysam.AlignedSegment()
                read.query_name = tid
                read.query_sequence = "A" * 12
                read.flag = flag
                read.reference_id = 0
                read.reference_start = 0
                read.cigartuples = [(0, 6), (3, 6), (0, 6)]
                out.write(read)
        annotation = self.root / "out.gff3"
        transcript_dedup(str(bam), {"g": ["g"]}, {"g.1": 12, "g.2": 12},
                         {"g": len(dna)}, {"g": bed_strand}, None, annotation, seed_gff3_path=seed)
        cdna = self.root / "cdna.fa"
        get_cdna_from_gff3(str(fasta), annotation, cdna)
        return load_gff3_transcript_models(annotation), [str(r.seq) for r in SeqIO.parse(cdna, "fasta")]

    def test_oriented_negative_chromosome_gene_is_not_reversed_twice(self):
        models, cdna = self.run_model("ATGAAANNNNNNCCCTAA", "-", [("g.1", 256)])
        self.assertEqual(models["transcript_context"]["g.1"], ("g", "+"))
        self.assertEqual(cdna, ["ATGAAACCCTAA"])

    def test_reverse_alignment_extracts_reverse_complement_regardless_of_bed(self):
        dna = str(Seq("ATGAAANNNNNNCCCTAA").reverse_complement())
        models, cdna = self.run_model(dna, "+", [("g.1", 16)])
        self.assertEqual(models["transcript_context"]["g.1"], ("g", "-"))
        self.assertEqual(cdna, ["ATGAAACCCTAA"])

    def test_opposite_strands_are_not_deduplicated(self):
        models, cdna = self.run_model("ATGAAANNNNNNCCCTAA", "-", [("g.1", 0), ("g.2", 16)])
        self.assertEqual(len(models["transcript_gene"]), 2)
        self.assertCountEqual(cdna, ["ATGAAACCCTAA", "TTAGGGTTTCAT"])
        self.assertEqual(models["gene_rows"]["g"][6], ".")

    def test_strand_follows_selected_alignment(self):
        models, cdna = self.run_model("ATGAAANNNNNNCCCTAA", "-", [("g.1", 272), ("g.1", 0)])
        self.assertEqual(models["transcript_context"]["g.1"][1], "+")
        self.assertEqual(cdna, ["ATGAAACCCTAA"])

    def test_historical_reference_relative_strand_is_preserved(self):
        seed = self.root / "seed.gff3"
        write_gff3(seed, [make_row("g", "gene", 1, 18, "-", "g"),
                         make_row("g", "mRNA", 1, 18, "-", "g.1", "g"),
                         make_row("g", "exon", 1, 6, "-", "g.1.e1", "g.1"),
                         make_row("g", "exon", 13, 18, "-", "g.1.e2", "g.1")])
        models, cdna = self.run_model("TTAGGGNNNNNNTTTCAT", "+", [], str(seed))
        self.assertEqual(models["transcript_context"]["g.1"][1], "-")
        self.assertEqual(cdna, ["ATGAAACCCTAA"])

    def test_unknown_strand_fails_instead_of_assuming_forward(self):
        seed = self.root / "unknown.gff3"
        write_gff3(seed, [make_row("g", "gene", 1, 6, ".", "g"),
                         make_row("g", "mRNA", 1, 6, ".", "g.1", "g"),
                         make_row("g", "exon", 1, 6, ".", "g.1.e1", "g.1")])
        fasta = self.root / "gene.fa"
        fasta.write_text(">g\nATGTAA\n")
        with self.assertRaisesRegex(ValueError, "strand must be"):
            get_cdna_from_gff3(str(fasta), seed, self.root / "out.fa")
