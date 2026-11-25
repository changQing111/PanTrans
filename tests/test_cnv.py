import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from pantrans.cnv import (
    GenePosition, detect_cnv, first_proximity_group, load_cnv_bed,
    select_cnv_representative,
)


class CnvTests(unittest.TestCase):
    def test_inclusive_gap_and_transitive_chain(self):
        positions = {
            "a": GenePosition("Chr1A", 0, 100, "CS"),
            "b": GenePosition("Chr1A", 1_000_100, 1_000_200, "CS"),
            "c": GenePosition("Chr1A", 2_000_200, 2_000_300, "CS"),
        }
        self.assertEqual(first_proximity_group(["c", "a", "b"], positions), ["a", "b", "c"])
        self.assertEqual(first_proximity_group(["a", "c"], positions), [])
        self.assertEqual(first_proximity_group(["a", "b"], positions, 999_999), [])

    def test_nested_overlapping_and_touching_intervals(self):
        positions = {
            "outer": GenePosition("Chr1A", 0, 5000, "CS"),
            "inner": GenePosition("Chr1A", 10, 20, "CS"),
            "touching": GenePosition("Chr1A", 5000, 6000, "CS"),
        }
        self.assertEqual(first_proximity_group(list(positions), positions, 0), list(positions))

    def test_requires_same_variety_and_chromosome(self):
        positions = {
            "a": GenePosition("Chr1A", 0, 100, "CS"),
            "b": GenePosition("Chr1A", 200, 300, "Q"),
            "c": GenePosition("Chr1B", 200, 300, "CS"),
        }
        self.assertEqual(first_proximity_group(list(positions), positions), [])

    def test_first_hit_only_and_skips_isolated_first_interval(self):
        positions = {
            "a": GenePosition("Chr1A", 0, 100, "Q"),
            "b": GenePosition("Chr1A", 3_000_000, 3_000_100, "Q"),
            "c": GenePosition("Chr1A", 3_000_150, 3_000_200, "Q"),
            "d": GenePosition("Chr2A", 0, 100, "CS"),
            "e": GenePosition("Chr2A", 200, 300, "CS"),
        }
        self.assertEqual(first_proximity_group(list(positions), positions), ["b", "c"])

    def test_reference_then_bed_length_then_id_across_entire_cluster(self):
        positions = {
            "short_ref": GenePosition("Chr1A", 0, 10, "CS"),
            "z_ref": GenePosition("Chr2A", 0, 100, "CS"),
            "a_ref": GenePosition("Chr3A", 0, 100, "CS"),
            "long_q": GenePosition("Chr1A", 0, 1000, "Q"),
            "short_q": GenePosition("Chr1A", 1500, 1600, "Q"),
        }
        self.assertEqual(select_cnv_representative(list(positions), positions, "CS"), "a_ref")
        self.assertEqual(select_cnv_representative(["short_q", "long_q"], positions, "CS"), "long_q")

    def test_output_once_per_cluster_and_preserved_on_invalid_input(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bed, clusters, output = root / "genes.bed", root / "pre.cluster", root / "result.txt"
            bed.write_text(
                "Chr1A\t0\t100\tq1\tQ\t+\n"
                "Chr1A\t200\t400\tq2\tQ\t-\n"
                "Chr2A\t0\t10\tr1\tCS\t+\n"
                "Chr2A\t20\t30\tr2\tCS\t+\n"
                "Chr3A\t0\t20\tsingle\tCS\t+\n"
            )
            clusters.write_text("# pre clusters\nq1 q2 r1 r2\n\nsingle\n")
            self.assertEqual(detect_cnv(clusters, bed, "CS", output), 1)
            self.assertEqual(output.read_text(), "r1\n")
            clusters.write_text("q1 q2\nunknown\n")
            with self.assertRaisesRegex(ValueError, "absent from BED"):
                detect_cnv(clusters, bed, "CS", output)
            self.assertEqual(output.read_text(), "r1\n")
            self.assertFalse(list(root.glob(".pantrans-cnv-*")))
            clusters.write_text("q1 q1\n")
            with self.assertRaisesRegex(ValueError, "duplicate gene"):
                detect_cnv(clusters, bed, "CS", output)
            clusters.write_text("single\n")
            self.assertEqual(detect_cnv(clusters, bed, "CS", output), 0)
            self.assertEqual(output.read_text(), "")
            with self.assertRaisesRegex(ValueError, "Reference variety"):
                detect_cnv(clusters, bed, "typo", output)
            with self.assertRaisesRegex(ValueError, "overwrite an input"):
                detect_cnv(clusters, bed, "CS", bed)

    def test_invalid_bed_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "genes.bed"
            for content in (
                "Chr1A\t0\t10\ta\t.\n",
                "Chr1A\t10\t0\ta\tCS\n",
                "Chr1A\t-1\t10\ta\tCS\n",
                "Chr1A\t0\t10\ta\tCS\nChr1B\t0\t10\ta\tCS\n",
            ):
                with self.subTest(content=content):
                    path.write_text(content)
                    with self.assertRaises(ValueError):
                        load_cnv_bed(path)

    def test_cli_success_and_error_exit_codes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bed, clusters, output = root / "genes.bed", root / "pre.cluster", root / "result.txt"
            bed.write_text("Chr1A\t0\t100\ta\tCS\nChr1A\t200\t400\tb\tCS\n")
            clusters.write_text("a b\n")
            env = dict(os.environ, PYTHONPATH=str(Path(__file__).resolve().parents[1] / "src"))
            command = [sys.executable, "-m", "pantrans", "cnv", "--pre-cluster", str(clusters),
                       "-b", str(bed), "-r", "CS", "-o", str(output)]
            result = subprocess.run(command, env=env, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(output.read_text(), "b\n")
            result = subprocess.run(command + ["--max-gap", "-1"], env=env, capture_output=True, text=True)
            self.assertEqual(result.returncode, 2)
            clusters.write_text("missing\n")
            result = subprocess.run(command, env=env, capture_output=True, text=True)
            self.assertEqual(result.returncode, 1)
            self.assertEqual(output.read_text(), "b\n")


if __name__ == "__main__":
    unittest.main()
