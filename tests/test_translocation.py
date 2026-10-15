import csv
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

from pantrans.translocation import detect_translocations


class TranslocationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.bed = self.root / "genes.bed"
        self.pre = self.root / "pre.cluster"
        self.edges = self.root / "graph.edges.tsv"
        self.main = self.root / "main_chroms"
        self.prefix = self.root / "result"
        self.main.write_text("Chr1A\nChr1B\nChr2A\nChr3B\n", encoding="utf-8")

    def fixture(self, query_starts=(0, 1000, 2000), ref_starts=(0, 1000, 2000), same_cluster=False):
        records, clusters, edges = [], [], []
        for index, (ref_start, query_start) in enumerate(zip(ref_starts, query_starts)):
            ref, query = f"r{index}", f"q{index}"
            records.extend([
                f"Chr1A\t{ref_start}\t{ref_start + 100}\t{ref}\tCS\t+",
                f"Chr3B\t{query_start}\t{query_start + 100}\t{query}\tQ\t-",
            ])
            clusters.append(f"{ref} {query}")
            edges.extend([f"{ref}\t{query}", f"{query}\t{ref}"])
        self.bed.write_text("\n".join(records) + "\n", encoding="utf-8")
        self.pre.write_text((" " if same_cluster else "\n").join(clusters) + "\n", encoding="utf-8")
        self.edges.write_text("\n".join(edges) + "\n", encoding="utf-8")

    def read_report(self, suffix):
        with open(str(self.prefix) + suffix, encoding="utf-8", newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def run_screen(self, **kwargs):
        summary = detect_translocations(self.pre, self.bed, self.edges, self.main, "CS", self.prefix, **kwargs)
        return summary, self.read_report(".translocation.genes.tsv"), self.read_report(".translocation.blocks.tsv")

    def test_three_independent_anchors_support_forward_block(self):
        self.fixture()
        summary, genes, blocks = self.run_screen()
        self.assertEqual(summary["supported_blocks"], 1)
        self.assertEqual(len(genes), 3)
        self.assertEqual({gene["status"] for gene in genes}, {"supported_block"})
        block = blocks[0]
        self.assertEqual(block["orientation"], "+")
        self.assertEqual(block["anchor_count"], "3")
        self.assertEqual(block["pre_cluster_count"], "3")
        self.assertEqual(block["reference_genes"], "r0;r1;r2")
        self.assertEqual(block["query_genes"], "q0;q1;q2")
        self.assertEqual((block["reference_start"], block["reference_end"]), ("0", "2100"))
        self.assertEqual({gene["block_id"] for gene in genes}, {block["block_id"]})

    def test_reverse_requires_explicit_option(self):
        self.fixture(query_starts=(2000, 1000, 0))
        summary, genes, blocks = self.run_screen()
        self.assertEqual(summary["supported_blocks"], 0)
        self.assertEqual(blocks, [])
        self.assertEqual({gene["status"] for gene in genes}, {"single_gene_candidate"})
        summary, genes, blocks = self.run_screen(allow_reverse=True)
        self.assertEqual(summary["supported_blocks"], 1)
        self.assertEqual(blocks[0]["orientation"], "-")
        self.assertEqual(blocks[0]["query_genes"], "q0;q1;q2")
        self.assertEqual((blocks[0]["query_start"], blocks[0]["query_end"]), ("0", "2100"))

    def test_gap_is_interval_based_inclusive_and_required_on_both_sides(self):
        starts = (0, 1_000_100, 2_000_200)
        self.fixture(query_starts=starts, ref_starts=starts)
        self.assertEqual(self.run_screen()[0]["supported_blocks"], 1)
        self.assertEqual(self.run_screen(max_gap=999_999)[2], [])
        self.fixture(query_starts=(0, 1_000_101, 2_000_202))
        self.assertEqual(self.run_screen()[2], [])
        self.fixture(ref_starts=(0, 1_000_101, 2_000_202))
        self.assertEqual(self.run_screen()[2], [])

    def test_two_independent_anchors_only_supported_when_threshold_lowered(self):
        self.fixture(query_starts=(0, 1000), ref_starts=(0, 1000))
        summary, _, blocks = self.run_screen()
        self.assertEqual(summary["supported_blocks"], 0)
        self.assertEqual(blocks[0]["status"], "low_support_block")
        self.assertEqual(self.run_screen(min_anchor_genes=2)[0]["supported_blocks"], 1)

    def test_mixed_order_cannot_pool_forward_and_reverse_support(self):
        self.fixture(query_starts=(0, 2000, 1000))
        summary, genes, blocks = self.run_screen(allow_reverse=True, max_skipped_genes=0)
        self.assertEqual(summary["supported_blocks"], 0)
        used = [gene for block in blocks for gene in block["reference_genes"].split(";")]
        self.assertEqual(len(set(used)), len(used))
        self.assertEqual({row["status"] for row in genes}, {"single_gene_candidate"})

    def test_bounded_chaining_can_skip_an_out_of_order_anchor(self):
        self.fixture(query_starts=(0, 2000, 1000, 3000), ref_starts=(0, 1000, 2000, 3000))
        summary, genes, blocks = self.run_screen(max_skipped_genes=1)
        self.assertEqual(summary["supported_blocks"], 1)
        self.assertEqual(blocks[0]["anchor_count"], "3")
        self.assertEqual(sum(row["status"] == "single_gene_candidate" for row in genes), 1)

    def test_intervening_genes_from_complete_bed_limit_chain(self):
        self.fixture()
        with self.bed.open("a") as handle:
            handle.write("Chr3B\t200\t300\textra\tQ\t+\n")
        summary, _, blocks = self.run_screen(max_skipped_genes=0)
        self.assertEqual(summary["supported_blocks"], 0)
        self.assertEqual(blocks[0]["query_genes"], "q1;q2")
        self.assertEqual(blocks[0]["status"], "low_support_block")
        self.assertEqual(self.run_screen(max_skipped_genes=1)[0]["supported_blocks"], 1)

    def test_one_large_family_cannot_supply_independent_block_support(self):
        self.fixture(same_cluster=True)
        summary, genes, blocks = self.run_screen()
        self.assertEqual(summary["supported_blocks"], 0)
        self.assertEqual(blocks[0]["anchor_count"], "3")
        self.assertEqual(blocks[0]["pre_cluster_count"], "1")
        self.assertEqual(blocks[0]["status"], "low_support_block")
        self.assertEqual({gene["status"] for gene in genes}, {"single_gene_candidate"})

    def test_same_chromosome_alternative_marks_cross_chromosome_copy_ambiguous(self):
        self.fixture()
        with self.bed.open("a") as handle:
            handle.write("Chr1A\t10\t110\tq_same\tQ\t+\n")
        self.pre.write_text("r0 q0 q_same\nr1 q1\nr2 q2\n")
        with self.edges.open("a") as handle:
            handle.write("r0\tq_same\nq_same\tr0\n")
        summary, genes, _ = self.run_screen()
        self.assertEqual(summary["same_chromosome_pairs"], 1)
        hit = next(row for row in genes if row["reference_gene"] == "r0")
        self.assertEqual(hit["status"], "ambiguous_copy")
        self.assertEqual(hit["reference_has_same_chrom_match"], "1")
        self.assertIn("same_chromosome_alternative", hit["reason"])
        self.assertEqual(hit["block_id"], ".")
        self.assertEqual(summary["supported_blocks"], 0)

    def test_query_with_multiple_reference_matches_is_ambiguous(self):
        self.fixture()
        with self.bed.open("a") as handle:
            handle.write("Chr2A\t0\t100\tr_other\tCS\t+\n")
        self.pre.write_text("r0 q0 r_other\nr1 q1\nr2 q2\n")
        with self.edges.open("a") as handle:
            handle.write("r_other\tq0\nq0\tr_other\n")
        _, genes, _ = self.run_screen()
        copies = [row for row in genes if row["query_gene"] == "q0"]
        self.assertEqual(len(copies), 2)
        self.assertTrue(all(row["status"] == "ambiguous_copy" for row in copies))
        self.assertTrue(all(row["query_match_count"] == "2" for row in copies))

    def test_normal_homeologs_do_not_form_translocation_block(self):
        self.bed.write_text("Chr1A\t0\t100\trA\tCS\nChr1B\t0\t100\trB\tCS\n"
                            "Chr1A\t0\t100\tqA\tQ\nChr1B\t0\t100\tqB\tQ\n")
        self.pre.write_text("rA rB qA qB\n")
        self.edges.write_text("".join(f"{r}\t{q}\n{q}\t{r}\n"
                                       for r in ("rA", "rB") for q in ("qA", "qB")))
        summary, genes, blocks = self.run_screen()
        self.assertEqual(summary["same_chromosome_pairs"], 2)
        self.assertEqual(len(genes), 2)
        self.assertTrue(all(row["status"] == "ambiguous_copy" for row in genes))
        self.assertEqual(blocks, [])

    def test_unplaced_either_side_is_excluded_unless_requested(self):
        self.fixture()
        self.bed.write_text(self.bed.read_text().replace("Chr3B\t0\t100\tq0", "contig\t0\t100\tq0")
                            .replace("Chr1A\t1000\t1100\tr1", "Ctg42\t1000\t1100\tr1"))
        summary, genes, blocks = self.run_screen()
        self.assertEqual(summary["excluded_unplaced_pairs"], 2)
        self.assertEqual(summary["reciprocal_pairs"], 1)
        self.assertNotIn("unplaced", summary["pair_statuses"])
        self.assertEqual([row["query_gene"] for row in genes], ["q2"])
        self.assertEqual(blocks, [])
        summary, genes, blocks = self.run_screen(include_unplaced=True)
        self.assertEqual(summary["excluded_unplaced_pairs"], 0)
        self.assertEqual(summary["pair_statuses"]["unplaced"], 2)
        self.assertEqual(blocks, [])
        self.assertTrue(all(row["block_id"] == "." for row in genes))

    def test_excluded_ctg_partners_do_not_create_copy_ambiguity(self):
        self.fixture()
        with self.bed.open("a") as handle:
            handle.write("contig\t0\t100\tquery_ctg\tQ\nCtg42\t0\t100\tref_ctg\tCS\n")
        self.pre.write_text("r0 q0 query_ctg\nr1 q1 ref_ctg\nr2 q2\n")
        with self.edges.open("a") as handle:
            handle.write("r0\tquery_ctg\nquery_ctg\tr0\nref_ctg\tq1\nq1\tref_ctg\n")
        summary, genes, _ = self.run_screen()
        self.assertEqual(summary["supported_blocks"], 1)
        self.assertEqual(summary["excluded_unplaced_pairs"], 2)
        self.assertTrue(all(row["reference_match_count"] == row["query_match_count"] == "1" for row in genes))
        summary, genes, _ = self.run_screen(include_unplaced=True)
        self.assertEqual(summary["supported_blocks"], 0)
        self.assertEqual(summary["pair_statuses"]["unplaced"], 2)
        self.assertEqual(summary["pair_statuses"]["ambiguous_copy"], 2)

    def test_one_way_indirect_cross_cluster_and_self_edges_are_not_anchors(self):
        self.fixture()
        self.edges.write_text("r0\tq0\nr1\tq1\nq1\tr1\nr2\tq2\nr0\tr0\n"
                              "q2\tr0\nr0\tq2\nq2\tq0\n")
        summary, genes, blocks = self.run_screen()
        self.assertEqual(summary["reciprocal_pairs"], 1)
        self.assertEqual([row["query_gene"] for row in genes], ["q1"])
        self.assertEqual(blocks, [])

    def test_other_varieties_are_compared_independently(self):
        self.fixture()
        with self.bed.open("a") as bed, self.edges.open("a") as edges:
            for index in range(3):
                bed.write(f"Chr2A\t{index * 1000}\t{index * 1000 + 100}\tz{index}\tZ\n")
                edges.write(f"r{index}\tz{index}\nz{index}\tr{index}\n")
        self.pre.write_text("\n".join(f"r{i} q{i} z{i}" for i in range(3)) + "\n")
        summary, genes, blocks = self.run_screen()
        self.assertEqual(summary["supported_blocks"], 2)
        self.assertEqual({row["query_variety"] for row in blocks}, {"Q", "Z"})
        self.assertTrue(all(row["reference_match_count"] == "1" for row in genes))

    def test_empty_results_have_headers_and_missing_reference_clusters_are_counted(self):
        self.fixture()
        self.pre.write_text("r0\nq0\nr1\nq1\nq2\n")
        summary, genes, blocks = self.run_screen()
        self.assertEqual(summary["clusters_without_reference"], 3)
        self.assertEqual(genes, [])
        self.assertEqual(blocks, [])
        self.assertTrue(Path(str(self.prefix) + ".translocation.genes.tsv").read_text().startswith("pre_cluster\t"))

    def test_input_order_and_duplicate_edges_do_not_change_reports(self):
        self.fixture()
        self.run_screen()
        before = [Path(str(self.prefix) + suffix).read_bytes() for suffix in
                  (".translocation.genes.tsv", ".translocation.blocks.tsv")]
        edges = self.edges.read_text().splitlines()
        self.edges.write_text("\n".join(reversed(edges + edges)) + "\n")
        self.bed.write_text("\n".join(reversed(self.bed.read_text().splitlines())) + "\n")
        self.pre.write_text("q0 r0\nq1 r1\nq2 r2\n")
        self.run_screen()
        after = [Path(str(self.prefix) + suffix).read_bytes() for suffix in
                 (".translocation.genes.tsv", ".translocation.blocks.tsv")]
        self.assertEqual(before, after)

    def test_overlaps_allowed_but_equal_start_positions_do_not_show_order(self):
        self.fixture(query_starts=(0, 50, 100), ref_starts=(0, 50, 100))
        self.assertEqual(self.run_screen(max_gap=0)[0]["supported_blocks"], 1)
        self.fixture(query_starts=(0, 0, 0))
        self.assertEqual(self.run_screen()[2], [])

    def test_validation_preserves_existing_reports(self):
        self.fixture()
        self.run_screen()
        gene_path = Path(str(self.prefix) + ".translocation.genes.tsv")
        before = gene_path.read_bytes()
        for content in ("r0 q0\nr0 q1\n", "r0 unknown\n", "r0 r0\n"):
            with self.subTest(content=content):
                self.pre.write_text(content)
                with self.assertRaises(ValueError):
                    self.run_screen()
                self.assertEqual(gene_path.read_bytes(), before)
        self.fixture()
        for content in ("r0\n", "r0\tunknown\n"):
            self.edges.write_text(content)
            with self.assertRaises(ValueError):
                self.run_screen()
        self.fixture()
        self.main.write_text("unknown_chromosome\n")
        with self.assertRaisesRegex(ValueError, "match --main-chroms"):
            self.run_screen()
        self.assertEqual(gene_path.read_bytes(), before)

    def test_publishing_failure_restores_both_previous_reports(self):
        self.fixture()
        self.run_screen()
        gene_path = Path(str(self.prefix) + ".translocation.genes.tsv")
        block_path = Path(str(self.prefix) + ".translocation.blocks.tsv")
        before = (gene_path.read_bytes(), block_path.read_bytes())
        self.edges.write_text("")
        real_replace = os.replace

        def fail_second(source, destination):
            if str(destination) == str(block_path):
                raise OSError("simulated failure")
            real_replace(source, destination)

        with mock.patch("pantrans.translocation.os.replace", side_effect=fail_second):
            with self.assertRaisesRegex(OSError, "simulated"):
                self.run_screen()
        self.assertEqual((gene_path.read_bytes(), block_path.read_bytes()), before)
        self.assertFalse(list(self.root.glob(".pantrans-translocation-*")))

    def test_cli_and_invalid_thresholds(self):
        self.fixture()
        env = dict(os.environ, PYTHONPATH=str(Path(__file__).resolve().parents[1] / "src"))
        command = [sys.executable, "-m", "pantrans", "translocation", "--pre-cluster", str(self.pre),
                   "-b", str(self.bed), "--edges", str(self.edges), "--main-chroms", str(self.main),
                   "-r", "CS", "-o", str(self.prefix)]
        result = subprocess.run(command, env=env, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.read_report(".translocation.blocks.tsv")[0]["status"], "supported_block")
        self.bed.write_text(self.bed.read_text().replace("Chr3B\t0\t100\tq0", "contig\t0\t100\tq0"))
        result = subprocess.run(command, env=env, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertFalse(any(row["status"] == "unplaced" for row in self.read_report(".translocation.genes.tsv")))
        result = subprocess.run(command + ["--include-unplaced"], env=env, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(sum(row["status"] == "unplaced" for row in self.read_report(".translocation.genes.tsv")), 1)
        for option, value in (("--min-anchor-genes", "1"), ("--max-gap", "-1"), ("--max-skipped-genes", "-1")):
            result = subprocess.run(command + [option, value], env=env, capture_output=True, text=True)
            self.assertEqual(result.returncode, 2, result.stderr)
        self.pre.write_text("missing\n")
        result = subprocess.run(command, env=env, capture_output=True, text=True)
        self.assertEqual(result.returncode, 1)

    def test_input_cannot_be_overwritten_by_report(self):
        self.fixture()
        self.edges = Path(str(self.prefix) + ".translocation.genes.tsv")
        self.edges.write_text("r0\tq0\nq0\tr0\n")
        before = self.edges.read_bytes()
        with self.assertRaisesRegex(ValueError, "overwrite an input"):
            self.run_screen()
        self.assertEqual(self.edges.read_bytes(), before)


if __name__ == "__main__":
    unittest.main()
