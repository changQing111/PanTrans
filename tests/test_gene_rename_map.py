import os
import tempfile
import unittest
import importlib.util
from pathlib import Path

from pantrans.pipeline import _build_gene_rename_map
from pantrans.gff3 import rename_gff3_ids, make_row, write_gff3, read_gff3


class GeneRenameMapTests(unittest.TestCase):
    def test_preserves_bed_chromosomes_and_transcript_suffixes(self):
        with tempfile.TemporaryDirectory() as directory:
            bed = os.path.join(directory, "genes.bed")
            chromosomes = ["Chr1A", "Chr1B", "A1", "1A", "scaffold_42"]
            genes = ["V.g%d" % i for i in range(len(chromosomes))]
            with open(bed, "w") as handle:
                for chrom, gene in zip(chromosomes, genes):
                    handle.write(f"{chrom}\t0\t100\t{gene}\tV\t+\n")
            mapping = _build_gene_rename_map(
                {gene: [gene] for gene in reversed(genes)}, bed, ["V"], "Pan",
                gene_variety_dic={gene: "V" for gene in genes},
            )
            self.assertEqual(mapping, {gene: f"Pan{chrom}000001" for gene, chrom in zip(genes, chromosomes)})
            source = os.path.join(directory, "source.gff3")
            output = os.path.join(directory, "official.gff3")
            write_gff3(source, [make_row("V.g0", "gene", 1, 100, "+", "V.g0"),
                                make_row("V.g0", "mRNA", 1, 100, "+", "V.g0.12", "V.g0"),
                                make_row("V.g0", "exon", 1, 100, "+", "V.g0.12.exon1", "V.g0.12")])
            rename_gff3_ids(source, output, mapping)
            with open(output) as handle:
                row = handle.read()
            self.assertTrue(row.startswith("##gff-version 3\nV.g0\t"))
            self.assertIn('ID=PanChr1A000001', row)
            self.assertIn('ID=PanChr1A000001.12;Parent=PanChr1A000001', row)

    def test_missing_representative_has_no_silent_contig_fallback(self):
        with tempfile.TemporaryDirectory() as directory:
            bed = os.path.join(directory, "genes.bed")
            with open(bed, "w") as handle:
                handle.write("Chr1A\t0\t100\tV.g1\tV\t+\n")
            with self.assertRaisesRegex(ValueError, "absent from BED"):
                _build_gene_rename_map({"V.missing": []}, bed, ["V"], "Pan")

    def test_repair_preserves_sequences_and_backs_up_previous_outputs(self):
        spec = importlib.util.spec_from_file_location(
            "repair_official_ids", Path(__file__).resolve().parents[1] / "scripts/repair_official_ids.py"
        )
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "genes.bed").write_text("Chr1A\t0\t4\tV.g1\tV\t+\n")
            (root / "test_last.tmp.cluster").write_text("V.g1\n")
            write_gff3(root / "test_unrenamed.gff3", [
                make_row("V.g1", "gene", 1, 4, "+", "V.g1"),
                make_row("V.g1", "mRNA", 1, 4, "+", "V.g1.2", "V.g1"),
                make_row("V.g1", "exon", 1, 4, "+", "V.g1.2.exon1", "V.g1.2")])
            source_fasta = b">V.g1.2 cDNA\nAcGT\n"
            (root / "test_unrenamed_cdna.refer.fasta").write_bytes(source_fasta)
            (root / "test.gff3").write_text("old GFF3")
            (root / "test_cdna.refer.fasta").write_text("old FASTA")
            report = module.repair(root / "genes.bed", root, "test")
            self.assertEqual(report["transcripts"], 1)
            self.assertEqual((Path(report["backup"]) / "test.gff3").read_text(), "old GFF3")
            self.assertEqual((Path(report["backup"]) / "test_cdna.refer.fasta").read_text(), "old FASTA")
            self.assertEqual((root / "test_cdna.refer.fasta").read_bytes(), b">PanChr1A000001.2 cDNA\nAcGT\n")
            self.assertEqual((root / "test_unrenamed_cdna.refer.fasta").read_bytes(), source_fasta)


if __name__ == "__main__":
    unittest.main()
