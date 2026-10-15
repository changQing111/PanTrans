import os
import tempfile
import unittest
import importlib.util
from pathlib import Path

from pantrans.pipeline import _build_gene_rename_map
from pantrans.transcript_processor import rename_gtf_ids


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
            source = os.path.join(directory, "source.gtf")
            output = os.path.join(directory, "official.gtf")
            with open(source, "w") as handle:
                handle.write('V.g0\tPan\texon\t1\t100\t.\t+\t.\tgene_id "V.g0"; transcript_id "V.g0.12";\n')
            rename_gtf_ids(source, output, mapping)
            with open(output) as handle:
                row = handle.read()
            self.assertTrue(row.startswith("V.g0\t"))
            self.assertIn('gene_id "PanChr1A000001"', row)
            self.assertIn('transcript_id "PanChr1A000001.12"', row)

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
            (root / "test_unrenamed.gtf").write_text('V.g1\tPan\texon\t1\t4\t.\t+\t.\tgene_id "V.g1"; transcript_id "V.g1.2";\n')
            source_fasta = b">V.g1.2 cDNA\nAcGT\n"
            (root / "test_unrenamed_cdna.refer.fasta").write_bytes(source_fasta)
            (root / "test.gtf").write_text("old GTF")
            (root / "test_cdna.refer.fasta").write_text("old FASTA")
            report = module.repair(root / "genes.bed", root, "test")
            self.assertEqual(report["transcripts"], 1)
            self.assertEqual((Path(report["backup"]) / "test.gtf").read_text(), "old GTF")
            self.assertEqual((Path(report["backup"]) / "test_cdna.refer.fasta").read_text(), "old FASTA")
            self.assertEqual((root / "test_cdna.refer.fasta").read_bytes(), b">PanChr1A000001.2 cDNA\nAcGT\n")
            self.assertEqual((root / "test_unrenamed_cdna.refer.fasta").read_bytes(), source_fasta)


if __name__ == "__main__":
    unittest.main()
