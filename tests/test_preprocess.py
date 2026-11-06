import csv
import json
import tempfile
import unittest
import subprocess
import sys
from pathlib import Path

from Bio.Seq import Seq
from Bio import SeqIO

from pantrans.preprocess import preprocess, read_list
from pantrans.main import read_parameters
from pantrans.common_func import get_bed, get_bed_varieties
from pantrans.align_filter import transcript_to_gene_id


class PreprocessTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.seq = "ACGTTGCAAGTCCTAGGCTAACGT"
        self.fa = self.root / "genome.fa"
        self.fa.write_text(">contig1\n" + self.seq + "\n")
        self.gff = self.root / "models.gff3"
        # Exons intentionally precede transcripts and genes; one shared exon.
        self.gff.write_text("\n".join([
            "##gff-version 3",
            "contig1\tx\texon\t2\t4\t.\t+\t.\tParent=isoform,gplus.2",
            "contig1\tx\texon\t7\t8\t.\t+\t.\tParent=isoform",
            "contig1\tx\texon\t18\t19\t.\t-\t.\tParent=negative",
            "contig1\tx\texon\t22\t23\t.\t-\t.\tParent=negative",
            "contig1\tx\tmRNA\t2\t8\t.\t+\t.\tID=isoform;Parent=gplus",
            "contig1\tx\tmRNA\t2\t4\t.\t+\t.\tID=gplus.2;Parent=gplus",
            "contig1\tx\tmRNA\t18\t23\t.\t-\t.\tID=negative;Parent=gminus",
            "contig1\tx\tgene\t2\t8\t.\t+\t.\tID=gplus",
            "contig1\tx\tgene\t18\t23\t.\t-\t.\tID=gminus",
        ]) + "\n")

    def run_preprocess(self, **kwargs):
        return preprocess([self.fa], [self.gff], ["V1"], self.root / "out", **kwargs)

    def read_table(self, path):
        with open(path) as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def fasta(self, path):
        with open(path) as handle:
            return {r.id: str(r.seq) for r in SeqIO.parse(handle, "fasta")}

    def test_clipped_flanks_splicing_and_local_coordinates(self):
        paths = self.run_preprocess(flank=3)
        gdna, cdna = self.fasta(paths["gdna"]), self.fasta(paths["cdna"])
        self.assertEqual(gdna["gplus"], self.seq[:11])
        self.assertEqual(gdna["gminus"], str(Seq(self.seq[14:]).reverse_complement()))
        self.assertEqual(cdna["gplus.1"], self.seq[1:4] + self.seq[6:8])
        self.assertEqual(cdna["gplus.2"], self.seq[1:4])
        self.assertEqual(cdna["gminus.1"], str(Seq(self.seq[17:19] + self.seq[21:23]).reverse_complement()))
        gene_rows = {r["gene"]: r for r in self.read_table(paths["genes"])}
        for gene in ("gplus", "gminus"):
            self.assertEqual(gene_rows[gene]["actual_upstream"], "1")
            self.assertEqual(gene_rows[gene]["actual_downstream"], "3")
        blocks = {}
        for row in self.read_table(paths["exons"]):
            self.assertEqual(row["gdna_strand"], "+")
            blocks.setdefault(row["transcript"], []).append(
                gdna[row["gene"]][int(row["gdna_start0"]):int(row["gdna_end0"])])
        self.assertEqual({t: "".join(e) for t, e in blocks.items()}, cdna)
        bed, strand = get_bed(paths["bed"])
        self.assertEqual(bed["gminus"], ["contig1", 17, 23])
        self.assertEqual(strand["gminus"], "-")
        self.assertEqual(get_bed_varieties(paths["bed"], required=True), {"gplus": "V1", "gminus": "V1"})
        self.assertEqual({transcript_to_gene_id(t) for t in cdna}, set(gdna))
        self.assertFalse(Path(str(self.fa) + ".fai").exists())
        self.assertTrue(Path(paths["gdna"] + ".fai").exists())

    def test_contig_shorter_than_flanks(self):
        paths = self.run_preprocess()
        self.assertTrue(all(len(s) == len(self.seq) for s in self.fasta(paths["gdna"]).values()))
        self.assertEqual(json.loads(Path(paths["manifest"]).read_text())["varieties"][0]["clipped_genes"], 2)

    def test_zero_flank(self):
        paths = self.run_preprocess(flank=0)
        self.assertEqual(self.fasta(paths["gdna"])["gplus"], self.seq[1:8])
        self.assertTrue(all(r["body_start0"] == "0" for r in self.read_table(paths["genes"])))

    def test_duplicate_ids_require_namespace_and_no_partial_outputs(self):
        with self.assertRaisesRegex(ValueError, "Duplicate gene"):
            preprocess([self.fa, self.fa], [self.gff, self.gff], ["V1", "V2"], self.root / "out")
        self.assertEqual(list((self.root / "out").iterdir()), [])
        paths = preprocess([self.fa, self.fa], [self.gff, self.gff], ["V1", "V2"], self.root / "out", namespace_ids=True)
        self.assertEqual(set(self.fasta(paths["gdna"])), {"V1__gplus", "V1__gminus", "V2__gplus", "V2__gminus"})

    def test_refuses_overwrite(self):
        paths = self.run_preprocess()
        original = Path(paths["cdna"]).read_bytes()
        with self.assertRaises(FileExistsError):
            self.run_preprocess()
        self.assertEqual(Path(paths["cdna"]).read_bytes(), original)

    def test_invalid_gene_boundary(self):
        self.gff.write_text(self.gff.read_text().replace("gene\t18\t23", "gene\t18\t40"))
        with self.assertRaisesRegex(ValueError, "out-of-bounds"):
            self.run_preprocess()
        self.assertEqual(list((self.root / "out").iterdir()), [])

    def test_missing_contig(self):
        self.gff.write_text(self.gff.read_text().replace("contig1", "unknown"))
        with self.assertRaisesRegex(ValueError, "missing contig"):
            self.run_preprocess()

    def test_overlapping_exons_rejected(self):
        with self.gff.open("a") as handle:
            handle.write("contig1\tx\texon\t3\t7\t.\t+\t.\tParent=isoform\n")
        with self.assertRaisesRegex(ValueError, "Overlapping exons"):
            self.run_preprocess()

    def test_orphan_exons_rejected(self):
        self.gff.write_text(self.gff.read_text().replace("Parent=negative", "Parent=missing"))
        with self.assertRaises(ValueError):
            self.run_preprocess()

    def test_duplicate_fasta_contig_rejected(self):
        with self.fa.open("a") as handle:
            handle.write(">contig1\nACGT\n")
        with self.assertRaisesRegex(ValueError, "duplicate FASTA"):
            self.run_preprocess()

    def test_cli_lists_relative_to_list_files(self):
        for filename, text in [("genomes.txt", "# genomes\ngenome.fa\n"),
                               ("annotations.txt", "models.gff3\n"), ("varieties.txt", "V1\n")]:
            (self.root / filename).write_text(text)
        args = read_parameters(["preprocess", "--genome-list", str(self.root / "genomes.txt"),
                                "--annotation-list", str(self.root / "annotations.txt"),
                                "--variety-list", str(self.root / "varieties.txt"), "-o", str(self.root / "out")])
        self.assertEqual(args.genomes, [str(self.fa)])
        self.assertEqual(args.annotations, [str(self.gff)])
        self.assertEqual(args.varieties, ["V1"])
        self.assertEqual(args.flank, 3000)

    def test_cli_rejects_list_length_mismatch(self):
        with self.assertRaises(SystemExit):
            read_parameters(["preprocess", "--genomes", "a", "b", "--annotations", "a", "--varieties", "V1", "-o", "out"])

    def test_cli_executes_and_reports_invalid_annotation(self):
        command = [sys.executable, "-m", "pantrans", "preprocess", "--genomes", str(self.fa),
                   "--annotations", str(self.gff), "--varieties", "V1", "-o", str(self.root / "cli")]
        done = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(done.returncode, 0, done.stderr)
        self.assertTrue((self.root / "cli/Pan_preprocess.json").is_file())
        self.gff.write_text(self.gff.read_text().replace("Parent=gminus", "Parent=missing"))
        command[-1] = str(self.root / "bad-cli")
        failed = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(failed.returncode, 1)
        self.assertIn("missing gene", failed.stderr)
        self.assertNotIn("Traceback", failed.stderr)


if __name__ == "__main__":
    unittest.main()
