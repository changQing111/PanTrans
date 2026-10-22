import tempfile
import unittest
from pathlib import Path
from unittest import mock

from pantrans.coding_annotation import annotate_gff3, run_transdecoder, resolve_transdecoder
from pantrans.main import read_parameters
from pantrans.pipeline import _write_cluster_outputs
from pantrans.gff3 import load_gff3_transcript_models, rename_gff3_ids, make_row, write_gff3, read_gff3


class CodingAnnotationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.gff3 = self.root / "models.gff3"
        self.fasta = self.root / "models.fa"
        self.gff = self.root / "orfs.gff3"

    def fixture(self, strand="+", sequence="CCATGAAACCCTAAG", start=3, end=14,
                kind="complete", phase=0, exons=((101, 107), (201, 208))):
        rows = [make_row("g", "gene", exons[0][0], exons[-1][1], strand, "g"),
                make_row("g", "mRNA", exons[0][0], exons[-1][1], strand, "g.1", "g")]
        rows += [make_row("g", "exon", a, b, strand, f"g.1.exon{i}", "g.1")
                 for i, (a, b) in enumerate(exons, 1)]
        write_gff3(self.gff3, rows)
        self.fasta.write_text(f">g.1\n{sequence}\n")
        self.gff.write_text(
            f"g.1\tTD\tmRNA\t1\t{len(sequence)}\t.\t+\t.\tID=g.1.p1;Name=ORF%20type:{kind}%20len:4\n"
            f"g.1\tTD\tCDS\t{start}\t{end}\t.\t+\t{phase}\tParent=g.1.p1\n")

    def features(self):
        result = {}
        rows = read_gff3(self.gff3)
        rows.sort(key=lambda r: int(r[3]), reverse=rows[0][6] == "-")
        for f in rows:
            result.setdefault(f[2], []).append((int(f[3]), int(f[4]), f[7]))
        return result

    def test_positive_spliced_cds_and_utrs(self):
        self.fixture(end=14, sequence="CCATGAAACCCTAAG")
        # ORF is ATG AAA CCC TAA at cDNA positions 3..14.
        annotate_gff3(self.gff3, self.fasta, self.gff)
        features = self.features()
        self.assertEqual(features["CDS"], [(103, 107, "0"), (201, 207, "1")])
        self.assertEqual(features["five_prime_UTR"], [(101, 102, ".")])
        self.assertEqual(features["stop_codon"], [(205, 207, ".")])
        self.assertEqual(features["three_prime_UTR"], [(208, 208, ".")])
        original = self.gff3.read_bytes()
        annotate_gff3(self.gff3, self.fasta, self.gff)
        self.assertEqual(self.gff3.read_bytes(), original)
        renamed = self.root / "official.gff3"
        rename_gff3_ids(self.gff3, renamed, {"g": "PanChr1000001"})
        original_rows, renamed_rows = read_gff3(self.gff3), read_gff3(renamed)
        self.assertEqual([r[:8] for r in original_rows], [r[:8] for r in renamed_rows])
        self.assertTrue(all(r[8]["ID"].startswith("PanChr1000001") for r in renamed_rows))
        self.assertTrue(all(p.startswith("PanChr1000001") for r in renamed_rows for p in r[8].get("Parent", [])))
        self.assertEqual(load_gff3_transcript_models(self.gff3)["exon_coords"]["g.1"], [(101, 107), (201, 208)])

    def test_negative_strand(self):
        self.fixture(strand="-")
        annotate_gff3(self.gff3, self.fasta, self.gff)
        f = self.features()
        self.assertEqual(f["CDS"], [(201, 206, "0"), (102, 107, "0")])
        self.assertEqual(f["five_prime_UTR"], [(207, 208, ".")])
        self.assertEqual(f["three_prime_UTR"], [(101, 101, ".")])
        self.assertEqual(f["stop_codon"], [(102, 104, ".")])

    def test_split_stop_codon(self):
        self.fixture(exons=((101, 113), (201, 202)))
        annotate_gff3(self.gff3, self.fasta, self.gff)
        self.assertEqual(self.features()["stop_codon"], [(112, 113, "."), (201, 201, ".")])
        self.assertEqual(self.features()["three_prime_UTR"], [(202, 202, ".")])

    def test_partials_only_have_supported_utr(self):
        for kind, expected, absent in [("5prime_partial", "three_prime_UTR", "five_prime_UTR"),
                                       ("3prime_partial", "five_prime_UTR", "three_prime_UTR")]:
            with self.subTest(kind=kind):
                self.fixture(kind=kind)
                annotate_gff3(self.gff3, self.fasta, self.gff)
                self.assertIn(expected, self.features())
                self.assertNotIn(absent, self.features())

    def test_internal_partial_phase(self):
        self.fixture(kind="internal", phase=2)
        annotate_gff3(self.gff3, self.fasta, self.gff)
        self.assertEqual(self.features()["CDS"], [(103, 107, "2"), (201, 207, "0")])
        self.assertNotIn("five_prime_UTR", self.features())
        self.assertNotIn("three_prime_UTR", self.features())

    def test_empty_predictions_preserve_transcripts(self):
        self.fixture()
        original = self.gff3.read_bytes()
        self.gff.write_text("##gff-version 3\n")
        self.assertEqual(annotate_gff3(self.gff3, self.fasta, self.gff), 0)
        self.assertEqual(self.gff3.read_bytes(), original)

    def test_invalid_prediction_does_not_modify_gff3(self):
        for edit in (lambda s: s.replace("\t14\t", "\t90\t"),
                     lambda s: s.replace("g.1", "unknown"),
                     lambda s: s + s,
                     lambda s: s.replace("\t+\t0\t", "\t-\t0\t")):
            self.fixture()
            original = self.gff3.read_bytes()
            self.gff.write_text(edit(self.gff.read_text()))
            with self.assertRaises(ValueError):
                annotate_gff3(self.gff3, self.fasta, self.gff)
            self.assertEqual(self.gff3.read_bytes(), original)

    def test_runner_uses_fresh_directory_and_sense_only(self):
        self.fixture()
        for executables in (["long", "predict"], ["wrapper"]):
            def fake_run(command, cwd, **kwargs):
                (Path(cwd) / "transcripts.fasta.transdecoder.gff3").write_text(self.gff.read_text())
            with mock.patch("pantrans.coding_annotation.resolve_transdecoder", return_value=executables), \
                 mock.patch("pantrans.coding_annotation.subprocess.run", side_effect=fake_run) as run:
                first = run_transdecoder(self.fasta, self.root / "runs")
                second = run_transdecoder(self.fasta, self.root / "runs")
                self.assertNotEqual(first.parent, second.parent)
                self.assertIn("-S", run.call_args_list[0].args[0])
                command = run.call_args_list[0].args[0]
                self.assertEqual(command[command.index("-m") + 1], "50")
                self.assertIn("--single_best_only", run.call_args.args[0])

    def test_missing_tool_and_missing_output_fail(self):
        with mock.patch("pantrans.coding_annotation.shutil.which", return_value=None):
            with self.assertRaises(FileNotFoundError):
                resolve_transdecoder()
        self.fixture()
        with mock.patch("pantrans.coding_annotation.resolve_transdecoder", return_value=["wrapper"]), \
             mock.patch("pantrans.coding_annotation.subprocess.run"):
            with self.assertRaises(RuntimeError):
                run_transdecoder(self.fasta, self.root / "runs")

    def test_pipeline_annotates_before_renaming(self):
        self.fixture()
        def dedup(*args, **kwargs):
            Path(kwargs["gff3_path"]).write_text(self.gff3.read_text())
        def cdna(gdna, gff3, output):
            Path(output).write_text(self.fasta.read_text())
        with mock.patch("pantrans.pipeline.transcript_dedup", side_effect=dedup), \
             mock.patch("pantrans.pipeline.get_cdna_from_gff3", side_effect=cdna), \
             mock.patch("pantrans.pipeline.extract_fasta_subset_by_names"), \
             mock.patch("pantrans.pipeline.get_subset_bed"), \
             mock.patch("pantrans.pipeline.run_transdecoder", return_value=self.gff) as run:
            _write_cluster_outputs({"g": ["g"]}, "gdna", "bed", "bam", {}, {}, {},
                                   {"g": "PanChr1000001"}, str(self.root), "Pan", "", annotate_cds=True)
        run.assert_called_once()
        self.assertEqual(run.call_args.kwargs["threads"], 1)
        official = (self.root / "Pan.gff3").read_text()
        unrenamed = (self.root / "Pan_unrenamed.gff3").read_text()
        self.assertIn("\tCDS\t", official)
        self.assertIn("\tCDS\t", unrenamed)
        self.assertIn('ID=PanChr1000001.1;Parent=PanChr1000001', official)
        self.assertIn('ID=g.1;Parent=g', unrenamed)

    def test_pipeline_parallel_uses_thread_budget(self):
        self.fixture()
        def dedup(*args, **kwargs):
            Path(kwargs["gff3_path"]).write_text(self.gff3.read_text())
        def cdna(gdna, gff3, output):
            Path(output).write_text(self.fasta.read_text())
        with mock.patch("pantrans.pipeline.transcript_dedup", side_effect=dedup), \
             mock.patch("pantrans.pipeline.get_cdna_from_gff3", side_effect=cdna), \
             mock.patch("pantrans.pipeline.extract_fasta_subset_by_names"), \
             mock.patch("pantrans.pipeline.get_subset_bed"), \
             mock.patch("pantrans.pipeline.run_transdecoder", return_value=self.gff) as run:
            _write_cluster_outputs({"g": ["g"]}, "gdna", "bed", "bam", {}, {}, {},
                                   {"g": "PanChr1000001"}, str(self.root), "Pan", "", annotate_cds=True, threads=7, transdecoder_parallel=True)
        run.assert_called_once()
        self.assertEqual(run.call_args.kwargs["threads"], 7)
        official = (self.root / "Pan.gff3").read_text()
        unrenamed = (self.root / "Pan_unrenamed.gff3").read_text()
        self.assertIn("\tCDS\t", official)
        self.assertIn("\tCDS\t", unrenamed)
        self.assertIn('ID=PanChr1000001.1;Parent=PanChr1000001', official)
        self.assertIn('ID=g.1;Parent=g', unrenamed)

    def test_cli(self):
        base = ["construct", "-c", "c", "-g", "g", "-b", "b", "-r", "r", "-p", "p", "-o", "o"]
        self.assertFalse(read_parameters(base).annotate_cds)
        self.assertFalse(read_parameters(base).transdecoder_parallel)
        args = read_parameters(base + ["--annotate-cds", "--transdecoder-parallel", "-t", "4"])
        self.assertTrue(args.transdecoder_parallel)
        self.assertEqual(args.threads, 4)
        self.assertEqual(read_parameters(base).min_orf_length, 50)
        self.assertEqual(read_parameters(base + ["--min-orf-length", "75"]).min_orf_length, 75)
        self.assertTrue(read_parameters(base + ["--transdecoder-dir", "/tools"]).annotate_cds)
        with self.assertRaises(SystemExit):
            read_parameters(base + ["--min-orf-length", "0"])


if __name__ == "__main__":
    unittest.main()
