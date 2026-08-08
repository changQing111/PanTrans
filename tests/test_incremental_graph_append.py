import json
import hashlib
import os
import tempfile
import unittest
from unittest import mock
import inspect

import pysam

from pantrans import pipeline
from pantrans.align_filter import (
    bam_alignment_provenance,
    filter_bam,
    iter_filtered_bam_gene_edges,
    merge_bams_with_full_header,
    validate_resume_bam,
)
from pantrans.graph_package import (
    edge_provenance_generations,
    iter_graph_edges,
    load_graph_package,
    merge_history_and_query_bed,
    validate_edge_provenance,
    write_graph_package,
)


class GraphPackageTest(unittest.TestCase):
    def _write_file(self, path, content):
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(content)

    def test_round_trip_serializes_unique_edges_and_metadata(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            cdna_path = os.path.join(temp_dir, "history.cdna.fa")
            gdna_path = os.path.join(temp_dir, "history.gdna.fa")
            bed_path = os.path.join(temp_dir, "history.bed")
            bam_path = os.path.join(temp_dir, "history.filtered.bam")
            manifest_path = os.path.join(temp_dir, "history.graph.json")
            self._write_file(cdna_path, ">Ref.g1.1\nAAAA\n")
            self._write_file(gdna_path, ">Ref.g1\nAAAA\n>Ref.g2\nAAA\n")
            self._write_file(
                bed_path,
                "chr1\t0\t4\tRef.g1\t0\t+\n"
                "chr1\t10\t13\tRef.g2\t0\t-\n",
            )
            with pysam.AlignmentFile(
                bam_path,
                "wb",
                header={
                    "HD": {"VN": "1.6"},
                    "SQ": [
                        {"SN": "Ref.g1", "LN": 4},
                        {"SN": "Ref.g2", "LN": 3},
                    ],
                },
            ):
                pass

            write_graph_package(
                manifest_path=manifest_path,
                edge_iter=[
                    ("Ref.g2", "Ref.g1"),
                    ("Ref.g1", "Ref.g2"),
                    ("Ref.g1", "Ref.g2"),
                ],
                gene_len_dic={"Ref.g1": 4, "Ref.g2": 3},
                bed_path=bed_path,
                filtered_bam_path=bam_path,
                cdna_paths=[cdna_path],
                gdna_paths=[gdna_path],
                variety_names=["Ref"],
                reference_name="Ref",
                main_chroms=["chr1"],
            )

            package = load_graph_package(manifest_path)
            self.assertEqual(
                list(iter_graph_edges(package["edges_path"])),
                [("Ref.g1", "Ref.g2"), ("Ref.g2", "Ref.g1")],
            )
            self.assertEqual(package["variety_names"], ["Ref"])
            self.assertEqual(package["reference_name"], "Ref")
            self.assertEqual(package["main_chroms"], ["chr1"])
            self.assertEqual(package["gene_len_dic"], {"Ref.g1": 4, "Ref.g2": 3})
            self.assertEqual(package["bed_path"], bed_path)
            self.assertEqual(package["filtered_bam_path"], bam_path)

    def test_load_rejects_changed_source_file(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            cdna_path = os.path.join(temp_dir, "history.cdna.fa")
            gdna_path = os.path.join(temp_dir, "history.gdna.fa")
            bed_path = os.path.join(temp_dir, "history.bed")
            bam_path = os.path.join(temp_dir, "history.filtered.bam")
            manifest_path = os.path.join(temp_dir, "history.graph.json")
            self._write_file(cdna_path, ">Ref.g1.1\nAAAA\n")
            self._write_file(gdna_path, ">Ref.g1\nAAAA\n")
            self._write_file(bed_path, "chr1\t0\t4\tRef.g1\t0\t+\n")
            with pysam.AlignmentFile(
                bam_path,
                "wb",
                header={"HD": {"VN": "1.6"}, "SQ": [{"SN": "Ref.g1", "LN": 4}]},
            ):
                pass
            write_graph_package(
                manifest_path,
                [("Ref.g1", "Ref.g1")],
                {"Ref.g1": 4},
                bed_path,
                bam_path,
                [cdna_path],
                [gdna_path],
                ["Ref"],
                "Ref",
                None,
            )
            self._write_file(cdna_path, ">Ref.g1.1\nAAAAAA\n")
            with self.assertRaisesRegex(ValueError, "source file changed"):
                load_graph_package(manifest_path)

    def test_rejects_edges_outside_bed_nodes(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            cdna_path = os.path.join(temp_dir, "history.cdna.fa")
            gdna_path = os.path.join(temp_dir, "history.gdna.fa")
            bed_path = os.path.join(temp_dir, "history.bed")
            bam_path = os.path.join(temp_dir, "history.filtered.bam")
            manifest_path = os.path.join(temp_dir, "history.graph.json")
            self._write_file(cdna_path, ">Ref.g1.1\nAAAA\n")
            self._write_file(gdna_path, ">Ref.g1\nAAAA\n")
            self._write_file(bed_path, "chr1\t0\t4\tRef.g1\t0\t+\n")
            with pysam.AlignmentFile(
                bam_path,
                "wb",
                header={"HD": {"VN": "1.6"}, "SQ": [{"SN": "Ref.g1", "LN": 4}]},
            ):
                pass

            with self.assertRaisesRegex(ValueError, "outside graph BED nodes"):
                write_graph_package(
                    manifest_path,
                    [("Missing.g1", "Ref.g1")],
                    {"Ref.g1": 4},
                    bed_path,
                    bam_path,
                    [cdna_path],
                    [gdna_path],
                    ["Ref"],
                    "Ref",
                    None,
                )

    def test_validates_and_preserves_each_edge_generation(self):
        current = {
            "minimap2_version": "2.28-r1209",
            "minimap2_options": "-ax splice:hq -uf -N 100",
            "filter_thresholds": {"coverage_min": 0.8},
            "filter_logic_id": "pantrans-filter-v1",
            "filter_thresholds_assumed": False,
        }
        historical = {
            "edge_generations": [
                dict(current, scope="historical_all_to_all"),
                dict(current, scope="first_append_cross_edges"),
            ]
        }
        generations = validate_edge_provenance(historical, current)
        self.assertEqual(
            [generation["scope"] for generation in generations],
            ["historical_all_to_all", "first_append_cross_edges"],
        )
        self.assertEqual(edge_provenance_generations(historical), generations)

        incompatible = {
            "edge_generations": [
                dict(current, minimap2_version="2.27-r1193")
            ]
        }
        with self.assertRaisesRegex(ValueError, "incompatible edge provenance"):
            validate_edge_provenance(incompatible, current)


class BedMergeTest(unittest.TestCase):
    def test_merges_full_history_with_representative_plus_query_bed(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            history_path = os.path.join(temp_dir, "history.bed")
            append_path = os.path.join(temp_dir, "append.bed")
            output_path = os.path.join(temp_dir, "merged.bed")
            with open(history_path, "w", encoding="utf-8") as handle:
                handle.write(
                    "chr1\t0\t4\tRef.g1\t0\t+\n"
                    "chr1\t10\t14\tRef.g2\t0\t+\n"
                )
            with open(append_path, "w", encoding="utf-8") as handle:
                handle.write(
                    "chr1\t0\t4\tRef.g1\t0\t+\n"
                    "chr1\t20\t24\tJM22.g1\t0\t-\n"
                )

            merge_history_and_query_bed(
                history_path,
                append_path,
                history_gene_ids={"Ref.g1", "Ref.g2"},
                query_gene_ids={"JM22.g1"},
                output_path=output_path,
            )

            with open(output_path, encoding="utf-8") as handle:
                rows = [line.rstrip("\n") for line in handle]
            self.assertEqual(
                rows,
                [
                    "chr1\t0\t4\tRef.g1\t0\t+",
                    "chr1\t10\t14\tRef.g2\t0\t+",
                    "chr1\t20\t24\tJM22.g1\t0\t-",
                ],
            )

    def test_rejects_append_bed_gene_outside_history_or_query(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            history_path = os.path.join(temp_dir, "history.bed")
            append_path = os.path.join(temp_dir, "append.bed")
            output_path = os.path.join(temp_dir, "merged.bed")
            with open(history_path, "w", encoding="utf-8") as handle:
                handle.write("chr1\t0\t4\tRef.g1\t0\t+\n")
            with open(append_path, "w", encoding="utf-8") as handle:
                handle.write("chr1\t20\t24\tUnexpected.g1\t0\t+\n")

            with self.assertRaisesRegex(ValueError, "not present in history or query"):
                merge_history_and_query_bed(
                    history_path,
                    append_path,
                    history_gene_ids={"Ref.g1"},
                    query_gene_ids={"JM22.g1"},
                    output_path=output_path,
                )


class BamMergeTest(unittest.TestCase):
    def _header(self, names):
        return {
            "HD": {"VN": "1.6"},
            "SQ": [{"SN": name, "LN": 100} for name in names],
        }

    def _write_alignment(self, path, header, reference_name):
        with pysam.AlignmentFile(path, "wb", header=header) as output:
            self._write_alignment_record(
                output,
                "Query.g1.1",
                reference_id=next(
                    index
                    for index, item in enumerate(header["SQ"])
                    if item["SN"] == reference_name
                ),
            )

    def _write_alignment_record(
        self, output, query_name, reference_id=None, unmapped=False
    ):
        read = pysam.AlignedSegment()
        read.query_name = query_name
        read.query_sequence = "A" * 20
        read.flag = 4 if unmapped else 0
        read.reference_id = -1 if unmapped else reference_id
        read.reference_start = -1 if unmapped else 0
        read.mapping_quality = 0 if unmapped else 60
        read.cigar = [] if unmapped else [(0, 20)]
        read.query_qualities = pysam.qualitystring_to_array("I" * 20)
        output.write(read)

    def test_remaps_subset_reference_ids_by_name(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            source_path = os.path.join(temp_dir, "subset.bam")
            header_path = os.path.join(temp_dir, "full.bam")
            output_path = os.path.join(temp_dir, "merged.bam")
            subset_header = self._header(["New.g1"])
            full_header = self._header(["Old.g1", "New.g1"])
            self._write_alignment(source_path, subset_header, "New.g1")
            with pysam.AlignmentFile(header_path, "wb", header=full_header):
                pass

            merge_bams_with_full_header([source_path], header_path, output_path)

            with pysam.AlignmentFile(output_path, "rb") as merged:
                self.assertEqual(list(merged.references), ["Old.g1", "New.g1"])
                record = next(merged.fetch(until_eof=True))
                self.assertEqual(merged.get_reference_name(record.reference_id), "New.g1")
                self.assertEqual(record.reference_id, 1)

    def test_rejects_reference_length_mismatch_during_merge(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            source_path = os.path.join(temp_dir, "subset.bam")
            header_path = os.path.join(temp_dir, "full.bam")
            output_path = os.path.join(temp_dir, "merged.bam")
            self._write_alignment(
                source_path,
                {"HD": {"VN": "1.6"}, "SQ": [{"SN": "New.g1", "LN": 101}]},
                "New.g1",
            )
            with pysam.AlignmentFile(
                header_path,
                "wb",
                header={"HD": {"VN": "1.6"}, "SQ": [{"SN": "New.g1", "LN": 100}]},
            ):
                pass

            with self.assertRaisesRegex(ValueError, "length_mismatch"):
                merge_bams_with_full_header([source_path], header_path, output_path)

    def test_streams_gene_edges_from_filtered_bam(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            bam_path = os.path.join(temp_dir, "filtered.bam")
            header = self._header(["Ref.g1", "Ref.g2"])
            with pysam.AlignmentFile(bam_path, "wb", header=header) as output:
                for query_name, reference_id in (
                    ("Query.g1.1", 0),
                    ("Query.g1.2", 1),
                ):
                    read = pysam.AlignedSegment()
                    read.query_name = query_name
                    read.query_sequence = "A" * 20
                    read.flag = 0
                    read.reference_id = reference_id
                    read.reference_start = 0
                    read.mapping_quality = 60
                    read.cigar = [(0, 20)]
                    read.query_qualities = pysam.qualitystring_to_array("I" * 20)
                    output.write(read)

            self.assertEqual(
                list(iter_filtered_bam_gene_edges(bam_path)),
                [("Query.g1", "Ref.g1"), ("Query.g1", "Ref.g2")],
            )

    def test_filter_bam_can_skip_edge_list_collection(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            input_path = os.path.join(temp_dir, "input.bam")
            output_path = os.path.join(temp_dir, "filtered.bam")
            header = self._header(["Ref.g1"])
            self._write_alignment(input_path, header, "Ref.g1")

            edges, lengths = filter_bam(
                input_path,
                output_path,
                {"Ref.g1": ["chr1", 0, 100]},
                collect_edges=False,
            )

            self.assertEqual(edges, [])
            self.assertEqual(lengths, {"Ref.g1": 100})
            with pysam.AlignmentFile(output_path, "rb") as filtered:
                self.assertEqual(sum(1 for _ in filtered.fetch(until_eof=True)), 1)
            with self.assertRaisesRegex(ValueError, "already PanTrans-filtered"):
                validate_resume_bam(
                    output_path,
                    {"Query.g1.1"},
                    {"Ref.g1": 100},
                    "query-to-all BAM",
                )

    def test_resume_bam_rejects_missing_expected_primary_query(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            bam_path = os.path.join(temp_dir, "partial.bam")
            self._write_alignment(bam_path, self._header(["Ref.g1"]), "Ref.g1")

            with self.assertRaisesRegex(ValueError, "missing expected query"):
                validate_resume_bam(
                    bam_path,
                    {"Query.g1.1", "Query.g2.1"},
                    {"Ref.g1": 100},
                    "query-to-all BAM",
                )

    def test_resume_bam_counts_primary_unmapped_query_as_present(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            bam_path = os.path.join(temp_dir, "with-unmapped.bam")
            header = self._header(["Ref.g1"])
            with pysam.AlignmentFile(bam_path, "wb", header=header) as output:
                self._write_alignment_record(output, "Query.g1.1", reference_id=0)
                self._write_alignment_record(output, "Query.g2.1", unmapped=True)

            validate_resume_bam(
                bam_path,
                {"Query.g1.1", "Query.g2.1"},
                {"Ref.g1": 100},
                "query-to-all BAM",
            )

    def test_resume_bam_rejects_same_id_with_different_query_sequence(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            bam_path = os.path.join(temp_dir, "stale-sequence.bam")
            self._write_alignment(bam_path, self._header(["Ref.g1"]), "Ref.g1")

            with self.assertRaisesRegex(ValueError, "query sequence does not match"):
                validate_resume_bam(
                    bam_path,
                    {
                        "Query.g1.1": {
                            "length": 20,
                            "sha256": hashlib.sha256(b"C" * 20).hexdigest(),
                        }
                    },
                    {"Ref.g1": 100},
                    "query-to-all BAM",
                )

    def test_resume_bam_accepts_matching_reverse_and_unmapped_sequences(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            bam_path = os.path.join(temp_dir, "reverse-and-unmapped.bam")
            header = self._header(["Ref.g1"])
            with pysam.AlignmentFile(bam_path, "wb", header=header) as output:
                reverse_read = pysam.AlignedSegment()
                reverse_read.query_name = "Query.g1.1"
                reverse_read.query_sequence = "GACT"
                reverse_read.flag = 16
                reverse_read.reference_id = 0
                reverse_read.reference_start = 0
                reverse_read.mapping_quality = 60
                reverse_read.cigar = [(0, 4)]
                reverse_read.query_qualities = pysam.qualitystring_to_array("IIII")
                output.write(reverse_read)

                unmapped_read = pysam.AlignedSegment()
                unmapped_read.query_name = "Query.g2.1"
                unmapped_read.query_sequence = "AACG"
                unmapped_read.flag = 4
                unmapped_read.reference_id = -1
                unmapped_read.reference_start = -1
                unmapped_read.mapping_quality = 0
                unmapped_read.cigar = []
                unmapped_read.query_qualities = pysam.qualitystring_to_array("IIII")
                output.write(unmapped_read)

            validate_resume_bam(
                bam_path,
                {
                    "Query.g1.1": {
                        "length": 4,
                        "sha256": hashlib.sha256(b"AGTC").hexdigest(),
                    },
                    "Query.g2.1": {
                        "length": 4,
                        "sha256": hashlib.sha256(b"AACG").hexdigest(),
                    },
                },
                {"Ref.g1": 100},
                "query-to-all BAM",
            )

    def test_reads_minimap2_provenance_from_bam_header(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            bam_path = os.path.join(temp_dir, "mapped.bam")
            header = self._header(["Ref.g1"])
            header["PG"] = [
                {
                    "ID": "minimap2",
                    "PN": "minimap2",
                    "VN": "2.28-r1209",
                    "CL": "minimap2 -ax splice:hq -uf -I 16G --secondary=yes -N 100 -t 32 gdna.fa cdna.fa",
                }
            ]
            with pysam.AlignmentFile(bam_path, "wb", header=header):
                pass

            provenance = bam_alignment_provenance(bam_path)
            self.assertEqual(provenance["minimap2_version"], "2.28-r1209")
            self.assertEqual(
                provenance["minimap2_options"],
                "-ax splice:hq -uf -I 16G --secondary=yes -N 100",
            )
            self.assertIsNone(provenance["pantrans_version"])


class ConstructGraphPackageTest(unittest.TestCase):
    def test_construct_writes_reusable_graph_package(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            cdna_path = os.path.join(temp_dir, "all.cdna.fa")
            gdna_path = os.path.join(temp_dir, "all.gdna.fa")
            bed_path = os.path.join(temp_dir, "all.bed")
            bam_path = os.path.join(temp_dir, "input.bam")
            for path, content in (
                (cdna_path, ">Ref.g1.1\nAAAA\n"),
                (gdna_path, ">Ref.g1\nAAAA\n"),
                (bed_path, "chr1\t0\t4\tRef.g1\t0\t+\n"),
            ):
                with open(path, "w", encoding="utf-8") as handle:
                    handle.write(content)
            with open(bam_path, "wb"):
                pass

            graph = mock.Mock()
            graph.number_of_nodes.return_value = 1
            graph.edges.return_value = [("Ref.g1", "Ref.g1")]
            with mock.patch.object(
                pipeline,
                "filter_bam",
                return_value=([("Ref.g1", "Ref.g1")], {"Ref.g1": 4}),
            ), mock.patch.object(
                pipeline,
                "_derive_graph_and_clusters",
                return_value=(graph, [["Ref.g1"]], [["Ref.g1"]]),
            ), mock.patch.object(
                pipeline,
                "_write_cluster_outputs",
                return_value=("out.gtf", "out.cdna", "out.gdna", "out.bed"),
            ), mock.patch.object(
                pipeline,
                "bam_alignment_provenance",
                return_value={
                    "graph_source": "provided_bam",
                    "minimap2_version": "test",
                    "minimap2_options": "test-options",
                    "filter_thresholds": {},
                },
            ), mock.patch.object(
                pipeline,
                "write_graph_package",
            ) as write_package:
                pipeline.unit_construct(
                    cdna_path,
                    gdna_path,
                    bed_path,
                    bam_path,
                    None,
                    ["Ref"],
                    1,
                    temp_dir,
                    "Pan",
                    "Ref",
                )

            write_package.assert_called_once()
            call = write_package.call_args.kwargs
            self.assertEqual(
                call["manifest_path"], os.path.join(temp_dir, "Pan.graph.json")
            )
            self.assertEqual(call["cdna_paths"], [cdna_path])
            self.assertEqual(call["gdna_paths"], [gdna_path])
            self.assertEqual(call["variety_names"], ["Ref"])
            self.assertEqual(call["reference_name"], "Ref")
            self.assertEqual(
                call["provenance"]["edge_generations"][0]["graph_source"],
                "provided_bam",
            )


class IncrementalAppendTest(unittest.TestCase):
    def test_validates_historical_gtf_and_cdna_identity(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            history_gtf = os.path.join(temp_dir, "history.gtf")
            history_cdna = os.path.join(temp_dir, "history.cdna.fa")
            with open(history_gtf, "w", encoding="utf-8") as handle:
                handle.write(
                    'Ref.g1\tPan\ttranscript\t1\t4\t.\t+\t.\tgene_id "Ref.g1"; transcript_id "Ref.g1.1";\n'
                    'Ref.g1\tPan\texon\t1\t4\t.\t+\t.\tgene_id "Ref.g1"; transcript_id "Ref.g1.1";\n'
                )
            with open(history_cdna, "w", encoding="utf-8") as handle:
                handle.write(">Ref.g1.1\nAAAA\n")

            pipeline._validate_history_transcriptome(
                history_gtf, history_cdna, {"Ref.g1", "Ref.g2"}
            )

            with open(history_cdna, "w", encoding="utf-8") as handle:
                handle.write(">Ref.g1.2\nAAAA\n")
            with self.assertRaisesRegex(ValueError, "transcript IDs do not match"):
                pipeline._validate_history_transcriptome(
                    history_gtf, history_cdna, {"Ref.g1", "Ref.g2"}
                )

    def test_rejects_history_gtf_genes_outside_graph_or_already_renamed(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            history_cdna = os.path.join(temp_dir, "history.cdna.fa")
            for gene_id, graph_genes, message in (
                ("Other.g1", {"Ref.g1"}, "absent from historical graph"),
                ("Pan1A000001", {"Pan1A000001"}, "renamed Pan"),
            ):
                with self.subTest(gene_id=gene_id):
                    history_gtf = os.path.join(temp_dir, f"{gene_id}.gtf")
                    transcript_id = f"{gene_id}.1"
                    with open(history_gtf, "w", encoding="utf-8") as handle:
                        handle.write(
                            f'{gene_id}\tPan\ttranscript\t1\t4\t.\t+\t.\tgene_id "{gene_id}"; transcript_id "{transcript_id}";\n'
                            f'{gene_id}\tPan\texon\t1\t4\t.\t+\t.\tgene_id "{gene_id}"; transcript_id "{transcript_id}";\n'
                        )
                    with open(history_cdna, "w", encoding="utf-8") as handle:
                        handle.write(f">{transcript_id}\nAAAA\n")
                    with self.assertRaisesRegex(ValueError, message):
                        pipeline._validate_history_transcriptome(
                            history_gtf, history_cdna, graph_genes
                        )

    def test_validates_gene_prefixes_before_append_alignment(self):
        pipeline._validate_variety_gene_ids({"JM22A1.g1"}, ["JM22"], "query")
        with self.assertRaisesRegex(ValueError, "does not match exactly one"):
            pipeline._validate_variety_gene_ids({"Other.g1"}, ["JM22"], "query")
        with self.assertRaisesRegex(ValueError, "duplicate variety"):
            pipeline._validate_variety_gene_ids({"JM22A1.g1"}, ["JM22", "JM22"], "query")
        with self.assertRaisesRegex(ValueError, "ambiguous variety prefixes"):
            pipeline._validate_variety_gene_ids({"JM22A1.g1"}, ["JM", "JM22"], "query")

    def test_append_uses_history_graph_and_two_cross_alignments(self):
        expected_parameters = [
            "history_cdna_path",
            "history_gtf_path",
            "query_cdna_path",
            "query_gdna_path",
            "all_bed_path",
            "history_graph_path",
            "variety_name",
            "threads",
            "out_dir",
            "prefix",
            "query_to_all_bam",
            "history_to_query_bam",
        ]
        self.assertEqual(
            list(inspect.signature(pipeline.unit_append).parameters),
            expected_parameters,
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            query_cdna = os.path.join(temp_dir, "query.cdna.fa")
            query_gdna = os.path.join(temp_dir, "query.gdna.fa")
            history_cdna = os.path.join(temp_dir, "history.final.cdna.fa")
            history_gtf = os.path.join(temp_dir, "history.final.gtf")
            raw_history_cdna = os.path.join(temp_dir, "history.raw.cdna.fa")
            history_gdna = os.path.join(temp_dir, "history.gdna.fa")
            history_bed = os.path.join(temp_dir, "history.bed")
            append_bed = os.path.join(temp_dir, "append.bed")
            history_bam = os.path.join(temp_dir, "history.filtered.bam")
            history_graph = os.path.join(temp_dir, "history.graph.json")
            out_dir = os.path.join(temp_dir, "out")
            for path in (
                query_cdna,
                query_gdna,
                history_cdna,
                raw_history_cdna,
                history_gdna,
                history_bed,
                append_bed,
                history_bam,
                history_graph,
            ):
                with open(path, "w", encoding="utf-8"):
                    pass
            with open(history_gtf, "w", encoding="utf-8") as handle:
                handle.write(
                    'Ref.g1\tPan\ttranscript\t1\t100\t.\t+\t.\tgene_id "Ref.g1"; transcript_id "Ref.g1.1";\n'
                    'Ref.g1\tPan\texon\t1\t100\t.\t+\t.\tgene_id "Ref.g1"; transcript_id "Ref.g1.1";\n'
                )

            package = {
                "cdna_paths": [raw_history_cdna],
                "gdna_paths": [history_gdna],
                "bed_path": history_bed,
                "filtered_bam_path": history_bam,
                "edges_path": os.path.join(temp_dir, "history.graph.edges.tsv"),
                "gene_len_dic": {"Ref.g1": 100},
                "history_gene_ids": {"Ref.g1"},
                "variety_names": ["Ref"],
                "reference_name": "Ref",
                "main_chroms": ["chr1"],
                "provenance": {
                    "edge_generations": [
                        {
                            "scope": "historical_all_to_all",
                            "minimap2_version": "test",
                            "minimap2_options": "test-options",
                            "filter_thresholds": {
                                "coverage_min": 0.8,
                                "identity_min": 0.9,
                                "soft_clip_max": 0.1,
                            },
                            "filter_logic_id": "pantrans-filter-v1",
                            "filter_thresholds_assumed": False,
                        }
                    ]
                },
            }
            current_provenance = {
                "minimap2_version": "test",
                "minimap2_options": "test-options",
                "filter_thresholds": {
                    "coverage_min": 0.8,
                    "identity_min": 0.9,
                    "soft_clip_max": 0.1,
                },
                "filter_logic_id": "pantrans-filter-v1",
                "filter_thresholds_assumed": False,
            }

            def concat(paths, output_path):
                with open(output_path, "w", encoding="utf-8") as output:
                    for path in paths:
                        with open(path, encoding="utf-8") as source:
                            output.write(source.read())

            def merge_bed(*args, **kwargs):
                output_path = kwargs.get("output_path", args[-1])
                with open(output_path, "w", encoding="utf-8") as output:
                    output.write("chr1\t0\t100\tRef.g1\t0\t+\n")
                    output.write("chr1\t100\t200\tJM22.g1\t0\t+\n")
                return output_path

            graph = mock.Mock()
            graph.number_of_nodes.return_value = 2
            with mock.patch.object(pipeline, "load_graph_package", return_value=package), \
                 mock.patch.object(pipeline, "iter_graph_edges", return_value=iter([("Ref.g1", "Ref.g1")])), \
                 mock.patch.object(
                     pipeline,
                     "iter_filtered_bam_gene_edges",
                     side_effect=[
                         iter([("JM22.g1", "Ref.g1")]),
                         iter([("Ref.g1", "JM22.g1")]),
                     ],
                 ), \
                 mock.patch.object(
                     pipeline, "concat_fasta_files", side_effect=concat
                 ) as concat_mock, \
                 mock.patch.object(pipeline, "merge_history_and_query_bed", side_effect=merge_bed), \
                 mock.patch.object(
                     pipeline,
                     "get_fasta_len",
                     side_effect=lambda path: {"JM22.g1": 100}
                     if path.endswith("query.gdna.fa")
                     else {"JM22.g1.1": 100}
                     if path.endswith("query.cdna.fa")
                     else {"Ref.g1.1": 100}
                     if path.endswith("history.final.cdna.fa")
                     else (
                         {"Ref.g1": 100, "JM22.g1": 100}
                         if path.endswith("merged.gdna.fasta")
                         else {"Ref.g1.1": 100, "JM22.g1.1": 100}
                     ),
                 ), \
                 mock.patch.object(
                     pipeline,
                     "get_bed",
                     return_value=({"Ref.g1": ["chr1", 0, 100], "JM22.g1": ["chr1", 100, 200]}, {"Ref.g1": "+", "JM22.g1": "+"}),
                 ), \
                 mock.patch.object(
                     pipeline,
                     "minimap2_map",
                 ) as minimap2_map, \
                 mock.patch.object(
                     pipeline,
                     "filter_bam",
                     side_effect=[
                         ([], {"Ref.g1": 100, "JM22.g1": 100}),
                         ([], {"JM22.g1": 100}),
                     ],
                 ) as filter_bam, \
                 mock.patch.object(pipeline, "merge_bams_with_full_header") as merge_bams, \
                 mock.patch.object(
                     pipeline,
                     "alignment_provenance",
                     return_value=current_provenance,
                 ), \
                 mock.patch.object(
                     pipeline,
                     "_derive_graph_and_clusters",
                     return_value=(graph, [["Ref.g1", "JM22.g1"]], [["Ref.g1", "JM22.g1"]]),
                 ) as derive_graph, \
                 mock.patch.object(
                     pipeline,
                     "_write_cluster_outputs",
                     return_value=("gtf", "cdna", "gdna", "bed"),
                 ) as write_outputs, \
                 mock.patch.object(pipeline, "write_graph_package") as write_package:
                pipeline.unit_append(
                    history_cdna_path=history_cdna,
                    history_gtf_path=history_gtf,
                    query_cdna_path=query_cdna,
                    query_gdna_path=query_gdna,
                    all_bed_path=append_bed,
                    history_graph_path=history_graph,
                    variety_name=["JM22"],
                    threads=1,
                    out_dir=out_dir,
                    prefix="Append",
                )

            self.assertEqual(minimap2_map.call_count, 2)
            self.assertEqual(
                concat_mock.call_args_list,
                [
                    mock.call(
                        [history_cdna, query_cdna],
                        os.path.join(out_dir, "Append_merged.cdna.fasta"),
                    ),
                    mock.call(
                        [history_gdna, query_gdna],
                        os.path.join(out_dir, "Append_merged.gdna.fasta"),
                    ),
                ],
            )
            self.assertEqual(
                minimap2_map.call_args_list,
                [
                    mock.call(
                        query_cdna,
                        os.path.join(out_dir, "Append_merged.gdna.fasta"),
                        1,
                        os.path.join(out_dir, "Append_query_to_all.bam"),
                    ),
                    mock.call(
                        history_cdna,
                        query_gdna,
                        1,
                        os.path.join(out_dir, "Append_history_to_query.bam"),
                    ),
                ],
            )
            self.assertEqual(filter_bam.call_count, 2)
            self.assertTrue(filter_bam.call_args_list[0].kwargs["collect_edges"] is False)
            self.assertTrue(filter_bam.call_args_list[1].kwargs["collect_edges"] is False)
            merge_bams.assert_called_once()
            merge_args = merge_bams.call_args.args
            self.assertEqual(
                merge_args[0],
                [
                    history_bam,
                    os.path.join(out_dir, "Append_query_to_all.filtered.bam"),
                    os.path.join(out_dir, "Append_history_to_query.filtered.bam"),
                ],
            )
            self.assertEqual(
                merge_args[1], os.path.join(out_dir, "Append_query_to_all.bam")
            )
            self.assertEqual(
                merge_args[2],
                os.path.join(out_dir, "Append_merged_cdna_align_gdna.filtered.bam"),
            )
            derive_graph.assert_called_once()
            self.assertEqual(
                list(derive_graph.call_args.kwargs["aligned_gene_li"]),
                [
                    ("Ref.g1", "Ref.g1"),
                    ("JM22.g1", "Ref.g1"),
                    ("Ref.g1", "JM22.g1"),
                ],
            )
            self.assertEqual(
                derive_graph.call_args.kwargs["variety_li"], ["Ref", "JM22"]
            )
            self.assertEqual(derive_graph.call_args.kwargs["refer_name"], "Ref")
            self.assertEqual(
                derive_graph.call_args.kwargs["eligible_gene_set"],
                {"Ref.g1", "JM22.g1"},
            )
            self.assertEqual(
                write_outputs.call_args_list[0].kwargs.get("seed_gtf_path"),
                None,
            )
            self.assertEqual(
                write_outputs.call_args_list[1].kwargs["seed_gtf_path"],
                history_gtf,
            )
            write_package.assert_called_once()
            self.assertEqual(
                write_package.call_args.kwargs["cdna_paths"],
                [history_cdna, query_cdna],
            )
            written_provenance = write_package.call_args.kwargs["provenance"]
            self.assertEqual(
                [generation["scope"] for generation in written_provenance["edge_generations"]],
                ["historical_all_to_all", "append_cross_alignments"],
            )


if __name__ == "__main__":
    unittest.main()
