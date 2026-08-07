import importlib.util
import json
import os
import tempfile
import unittest

import pysam


ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SCRIPT_PATH = os.path.join(ROOT, "scripts", "compare_append_to_construct.py")
SLURM_PATH = os.path.join(
    ROOT, "validation", "run_incremental_graph_append_JM22.slurm"
)
SPEC = importlib.util.spec_from_file_location("compare_append_to_construct", SCRIPT_PATH)
COMPARE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(COMPARE)


class CompareAppendToConstructTest(unittest.TestCase):
    def _write_bam(self, path):
        header = {
            "HD": {"VN": "1.6"},
            "SQ": [{"SN": name, "LN": 100} for name in ("A", "B", "C")],
        }
        with pysam.AlignmentFile(path, "wb", header=header) as output:
            for query, target in (("A.1", "A"), ("A.1", "B"), ("B.1", "A")):
                read = pysam.AlignedSegment()
                read.query_name = query
                read.query_sequence = "A" * 20
                read.flag = 0
                read.reference_id = ("A", "B", "C").index(target)
                read.reference_start = 0
                read.mapping_quality = 60
                read.cigar = [(0, 20)]
                read.query_qualities = pysam.qualitystring_to_array("I" * 20)
                output.write(read)

    def test_reports_representative_ids_independently_of_member_digest(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            construct_bam = os.path.join(temp_dir, "construct.bam")
            incremental_edges = os.path.join(temp_dir, "incremental.edges.tsv")
            incremental_pre = os.path.join(temp_dir, "incremental.pre.cluster")
            construct_pre = os.path.join(temp_dir, "construct.pre.cluster")
            incremental_last = os.path.join(temp_dir, "incremental.last.cluster")
            construct_last = os.path.join(temp_dir, "construct.last.cluster")
            output_dir = os.path.join(temp_dir, "comparison")
            self._write_bam(construct_bam)
            with open(incremental_edges, "w", encoding="utf-8") as output:
                output.write("A\tA\nA\tB\nB\tA\nB\tC\n")
            for path, content in (
                (incremental_pre, "A\tA\tB\nC\tC\n"),
                (construct_pre, "A\tB\tA\nD\tD\n"),
                (incremental_last, "A\tA\n"),
                (construct_last, "A\tB\n"),
            ):
                with open(path, "w", encoding="utf-8") as output:
                    output.write(content)

            COMPARE.main(
                [
                    "--incremental-edges", incremental_edges,
                    "--construct-bam", construct_bam,
                    "--incremental-pre", incremental_pre,
                    "--construct-pre", construct_pre,
                    "--incremental-last", incremental_last,
                    "--construct-last", construct_last,
                    "--output-dir", output_dir,
                ]
            )

            with open(
                os.path.join(output_dir, "comparison_summary.json"),
                encoding="utf-8",
            ) as handle:
                summary = json.load(handle)
            self.assertEqual(summary["directed_edges"]["common"], 3)
            self.assertEqual(summary["reciprocal_pairs"]["common"], 1)
            self.assertEqual(summary["clusters"]["pre"]["common_representatives"], 1)
            self.assertEqual(
                summary["clusters"]["pre"]["exact_members_for_common_representative"],
                1,
            )
            self.assertEqual(summary["clusters"]["last"]["common_representatives"], 1)
            self.assertEqual(
                summary["clusters"]["last"]["exact_members_for_common_representative"],
                0,
            )
            with open(
                os.path.join(output_dir, "last.incremental_only_representatives.txt"),
                encoding="utf-8",
            ) as handle:
                self.assertEqual(handle.read(), "")


class IncrementalAppendSlurmTest(unittest.TestCase):
    def test_uses_historical_unrenamed_transcriptome_and_query_cdna(self):
        with open(SLURM_PATH, encoding="utf-8") as handle:
            script = handle.read()

        self.assertIn("#SBATCH --cpus-per-task=32", script)
        self.assertIn("#SBATCH --mem=240G", script)
        self.assertIn("#SBATCH --time=3-00:00:00", script)
        self.assertIn("--coverage-min 0.80", script)
        self.assertIn("--identity-min 0.90", script)
        self.assertIn("--soft-clip-max 0.10", script)
        self.assertIn("history_construct_dir=", script)
        self.assertIn("-m pantrans.main construct", script)
        self.assertIn(
            "test_pantrans_cdna_align_gdna.filtered.bam", script
        )
        self.assertIn("history_manifest=", script)
        self.assertIn("history_gtf=", script)
        self.assertIn("_unrenamed.gtf", script)
        self.assertIn("history_cdna=", script)
        self.assertIn("_unrenamed_cdna.refer.fasta", script)
        self.assertIn('--cdna "${history_cdna}"', script)
        self.assertIn('--history-gtf "${history_gtf}"', script)
        self.assertIn("--query-cdna /data/changq/PanTrans/test/JM22_cdna.fasta", script)
        self.assertIn('--history-graph "${history_manifest}"', script)
        self.assertIn("test_pantrans_pre.refer_append_JM22.bed", script)
        self.assertIn("--incremental-edges", script)
        self.assertIn("--incremental-pre", script)
        self.assertIn("--construct-pre", script)
        self.assertIn("--incremental-last", script)
        self.assertIn("--construct-last", script)

    def test_feature_version_is_consistent(self):
        with open(os.path.join(ROOT, "pyproject.toml"), encoding="utf-8") as handle:
            pyproject = handle.read()
        with open(
            os.path.join(ROOT, "src", "pantrans", "version.py"),
            encoding="utf-8",
        ) as handle:
            version_module = handle.read()

        self.assertIn('version = "0.3.0"', pyproject)
        self.assertIn('__version__ = "0.3.0"', version_module)


if __name__ == "__main__":
    unittest.main()
