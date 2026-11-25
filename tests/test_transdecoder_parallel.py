import io
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest import mock

from pantrans.transdecoder_parallel import HEADER, records, score_parallel
from pantrans.coding_annotation import _parallel_predict


class ParallelScoringTest(unittest.TestCase):
    def test_order_header_and_coverage(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fasta = root / "input.fa"
            fasta.write_bytes(b">z desc\nATG\nAAA\n>a\nATGNNNTAA\n>b\nGC")
            model = root / "model"
            observed = []
            def scorer(command, stdout, stderr, check):
                self.assertEqual(command[-1], str(model))
                observed.append(command[1])
                stdout.write(HEADER)
                for record in records(command[1]):
                    stdout.write(record[0][1:].split()[0] + b"\t5\t9\t1\t2\t3\t4\t5\t6\n")
            output = io.BytesIO()
            with mock.patch("pantrans.transdecoder_parallel.subprocess.run", side_effect=scorer):
                score_parallel(fasta, model, "scorer", 8, output)
            rows = output.getvalue().splitlines()
            self.assertEqual(rows[0] + b"\n", HEADER)
            self.assertEqual([row.split(b"\t")[0] for row in rows[1:]], [b"z", b"a", b"b"])
            self.assertLessEqual(len(observed), 8)

    def test_failed_or_truncated_worker_never_publishes_scores(self):
        for failure in (True, False):
            with self.subTest(failure=failure), tempfile.TemporaryDirectory() as tmp:
                fasta = Path(tmp) / "input.fa"
                fasta.write_bytes(b">one\nATG\n>two\nTAA\n")
                def scorer(command, stdout, **kwargs):
                    if failure:
                        raise subprocess.CalledProcessError(1, command)
                    stdout.write(HEADER)
                output = io.BytesIO()
                with mock.patch("pantrans.transdecoder_parallel.subprocess.run", side_effect=scorer):
                    with self.assertRaises((ValueError, subprocess.CalledProcessError)):
                        score_parallel(fasta, "model", "scorer", 2, output)
                self.assertEqual(output.getvalue(), b"")

    def test_empty_input(self):
        with tempfile.TemporaryDirectory() as tmp:
            fasta = Path(tmp) / "empty.fa"
            fasta.touch()
            output = io.BytesIO()
            score_parallel(fasta, "model", "scorer", 2, output)
            self.assertEqual(output.getvalue(), HEADER)

    def test_unknown_predict_falls_back_without_changes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            predict = root / "Predict"
            predict.write_text("unknown version")
            self.assertIsNone(_parallel_predict(["LongOrfs", str(predict)], root, 4))
            self.assertEqual(predict.read_text(), "unknown version")
            self.assertIsNone(_parallel_predict(["wrapper"], root, 4))

    def test_both_modes_use_identical_tie_breaking(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            tools = root / "tools"
            (tools / "util").mkdir(parents=True)
            predict = tools / "Predict"
            predict.write_text('my $VERSION = "5.7.1";\n'
                               '$FindBin::RealBin $FindBin::RealBin\n'
                               '$UTIL_DIR/select_best_ORFs_per_transcript.pl\n'
                               '$UTIL_DIR/score_CDS_likelihood_all_6_frames.pl $cds_file $hexamer_scores_file > $cds_scores_file')
            selector = tools / "util/select_best_ORFs_per_transcript.pl"
            original = ('$FindBin::Bin/../PerlLib\n$b->{length} <=> $a->{length}\n'
                        'if ($SINGLE_BEST_ORF_FLAG) {\n$b->{length} <=> $a->{length}\n}')
            selector.write_text(original)
            (tools / "util/score_CDS_likelihood_all_6_frames.pl").touch()
            copies = []
            for threads in (1, 4):
                work = root / str(threads)
                work.mkdir()
                adapter, env = _parallel_predict(["long", str(predict)], work, threads)
                copies.append((work / "select_best_ORFs_per_transcript.deterministic.pl").read_text())
                self.assertIn("PANTRANS_SELECT_COMMAND", Path(adapter).read_text())
                self.assertEqual("PANTRANS_SCORE_COMMAND" in Path(adapter).read_text(), threads > 1)
            self.assertEqual(copies[0], copies[1])
            self.assertEqual(copies[0].count('Model_feat_name} cmp'), 2)
            self.assertIn('$b->{cds_scores}->[0] <=> $a->{cds_scores}->[0]', copies[0])
            self.assertEqual(selector.read_text(), original)
