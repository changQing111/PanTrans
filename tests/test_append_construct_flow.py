import os
import inspect
import sys
import tempfile
import types
import unittest
from unittest import mock


ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")


def dependency_stubs():
    pysam_stub = types.ModuleType("pysam")
    networkx_stub = types.ModuleType("networkx")
    networkx_stub.DiGraph = lambda: None
    networkx_stub.strongly_connected_components = lambda graph: []

    bio_stub = types.ModuleType("Bio")
    bio_stub.__path__ = []
    bio_seqio_stub = types.ModuleType("Bio.SeqIO")
    bio_seq_stub = types.ModuleType("Bio.Seq")
    bio_seq_stub.Seq = str
    bio_seqrecord_stub = types.ModuleType("Bio.SeqRecord")
    bio_seqrecord_stub.SeqRecord = object
    return {
        "pysam": pysam_stub,
        "networkx": networkx_stub,
        "Bio": bio_stub,
        "Bio.SeqIO": bio_seqio_stub,
        "Bio.Seq": bio_seq_stub,
        "Bio.SeqRecord": bio_seqrecord_stub,
    }


sys.path.insert(0, SRC)


def import_pipeline_with_dependency_stubs():
    preexisting_pantrans_modules = {
        module_name: module
        for module_name, module in sys.modules.items()
        if module_name == "pantrans" or module_name.startswith("pantrans.")
    }
    try:
        with mock.patch.dict(sys.modules, dependency_stubs()):
            from pantrans import pipeline as imported_pipeline
    finally:
        # Retain the direct reference below, but remove only PanTrans modules
        # introduced with fake globals and restore entries that predated this test.
        for module_name in tuple(sys.modules):
            is_pantrans_module = (
                module_name == "pantrans" or module_name.startswith("pantrans.")
            )
            if is_pantrans_module and module_name not in preexisting_pantrans_modules:
                del sys.modules[module_name]
        sys.modules.update(preexisting_pantrans_modules)
    return imported_pipeline


pipeline = import_pipeline_with_dependency_stubs()


def import_main_with_pipeline_stub():
    preexisting_pantrans_modules = {
        module_name: module
        for module_name, module in sys.modules.items()
        if module_name == "pantrans" or module_name.startswith("pantrans.")
    }
    try:
        with mock.patch.dict(sys.modules, {"pantrans.pipeline": pipeline}):
            from pantrans import main as imported_main
    finally:
        for module_name in tuple(sys.modules):
            is_pantrans_module = (
                module_name == "pantrans" or module_name.startswith("pantrans.")
            )
            if is_pantrans_module and module_name not in preexisting_pantrans_modules:
                del sys.modules[module_name]
        sys.modules.update(preexisting_pantrans_modules)
    return imported_main


pantrans_main = import_main_with_pipeline_stub()


class AppendCliContractTest(unittest.TestCase):
    def test_append_accepts_combined_bed_and_optional_bam(self):
        args = pantrans_main.read_parameters(
            [
                "append",
                "--name", "JM22",
                "--cdna", "query.cdna.fasta",
                "--gdna", "query.gdna.fasta",
                "--bed", "combined.bed",
                "--refer_cdna", "reference.cdna.fasta",
                "--refer_gdna", "reference.gdna.fasta",
                "--bam", "existing.bam",
                "--output", "append_out",
            ]
        )

        self.assertEqual(args.bed, "combined.bed")
        self.assertEqual(args.bam, "existing.bam")
        self.assertFalse(hasattr(args, "refer_bed"))
        self.assertFalse(hasattr(args, "refer_cluster"))

    def test_append_rejects_stale_reference_bed_option(self):
        with self.assertRaises(SystemExit):
            pantrans_main.read_parameters(
                [
                    "append",
                    "--name", "JM22",
                    "--cdna", "query.cdna.fasta",
                    "--gdna", "query.gdna.fasta",
                    "--bed", "combined.bed",
                    "--refer_cdna", "reference.cdna.fasta",
                    "--refer_gdna", "reference.gdna.fasta",
                    "--refer_bed", "obsolete.bed",
                    "--output", "append_out",
                ]
            )

    def test_construct_cli_contract_is_unchanged(self):
        args = pantrans_main.read_parameters(
            [
                "construct",
                "--name", "RefA",
                "--cdna", "all.cdna.fasta",
                "--gdna", "all.gdna.fasta",
                "--bed", "all.bed",
                "--reference", "RefA",
                "--prefix", "Pan",
                "--output", "construct_out",
            ]
        )

        self.assertEqual(args.subcommand, "construct")
        self.assertIsNone(args.bam)


class AppendFlowUnitTest(unittest.TestCase):
    def test_shared_clustering_recovers_pre_and_last_independently(self):
        graph = mock.Mock()
        graph.number_of_nodes.return_value = 2
        with mock.patch.object(pipeline, "di_graph_from_pair", return_value=graph), \
             mock.patch.object(
                 pipeline,
                 "get_conn_comp",
                 return_value=[("Ref.g1", "New.g1")],
             ), \
             mock.patch.object(
                 pipeline,
                 "assign_sccs",
                 return_value=([['Ref.g1', 'New.g1']], [['Ref.g1']]),
             ) as assign_mock:
            pre, last = pipeline._derive_clusters_from_alignment(
                aligned_gene_li=[("Ref.g1", "New.g1")],
                gene_len_dic={"Ref.g1": 100, "New.g1": 90, "New.g2": 80},
                bed_dic={
                    "Ref.g1": [],
                    "New.g1": [],
                    "New.g2": [],
                    "Ref.nonrep": [],
                },
                variety_li=["New"],
                refer_name="Ref",
                eligible_gene_set={"Ref.g1", "New.g1", "New.g2"},
            )

        assign_mock.assert_called_once()
        self.assertEqual(pre, [["Ref.g1", "New.g1"], ["New.g2"]])
        self.assertEqual(last, [["Ref.g1"], ["New.g1"], ["New.g2"]])
        self.assertNotIn("Ref.nonrep", [gene for cluster in pre + last for gene in cluster])

    def test_construct_recovers_pre_and_last_clusters_independently(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            cdna_path = os.path.join(temp_dir, "all.cdna.fasta")
            gdna_path = os.path.join(temp_dir, "all.gdna.fasta")
            bed_path = os.path.join(temp_dir, "all.bed")
            bam_path = os.path.join(temp_dir, "input.bam")
            with open(cdna_path, "w", encoding="utf-8") as handle:
                handle.write(
                    ">Ref.hub.1\nAAAA\n"
                    ">Ref.member.1\nAAAA\n"
                    ">Ref.absent.1\nAAAA\n"
                )
            with open(gdna_path, "w", encoding="utf-8") as handle:
                handle.write(
                    ">Ref.hub\nAAAA\n"
                    ">Ref.member\nAAAA\n"
                    ">Ref.absent\nAAAA\n"
                )
            with open(bed_path, "w", encoding="utf-8") as handle:
                handle.write(
                    "chr1\t0\t4\tRef.hub\t0\t+\n"
                    "chr1\t10\t14\tRef.member\t0\t+\n"
                    "chr1\t20\t24\tRef.absent\t0\t+\n"
                )

            graph = mock.Mock()
            graph.number_of_nodes.return_value = 2
            with mock.patch.object(
                pipeline,
                "filter_bam",
                return_value=([("Ref.hub", "Ref.member")], {"Ref.hub": 4}),
            ), mock.patch.object(
                pipeline,
                "di_graph_from_pair",
                return_value=graph,
            ), mock.patch.object(
                pipeline,
                "get_conn_comp",
                return_value=[("Ref.hub", "Ref.member")],
            ), mock.patch.object(
                pipeline,
                "assign_sccs",
                return_value=([['Ref.hub', 'Ref.member']], [['Ref.hub']]),
            ), mock.patch.object(
                pipeline,
                "_write_cluster_outputs",
                return_value=("output.gtf", "output.cdna", "output.gdna", "output.bed"),
            ) as write_outputs:
                pipeline.unit_construct(
                    cdna_path,
                    gdna_path,
                    bed_path,
                    bam_path,
                    None,
                    ["Ref"],
                    1,
                    temp_dir,
                    "test",
                    "Ref",
                )

            pre_clusters = list(write_outputs.call_args_list[0].args[0].values())
            last_clusters = list(write_outputs.call_args_list[1].args[0].values())
            self.assertIn(["Ref.hub", "Ref.member"], pre_clusters)
            self.assertNotIn(["Ref.member"], pre_clusters)
            self.assertIn(["Ref.absent"], pre_clusters)
            self.assertIn(["Ref.member"], last_clusters)
            self.assertIn(["Ref.absent"], last_clusters)

    def test_reference_prefix_inference_excludes_query_prefixes(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            bed_path = os.path.join(temp_dir, "combined.bed")
            with open(bed_path, "w", encoding="utf-8") as handle:
                handle.write(
                    "A1\t1\t100\tRefA.g1\t0\t+\n"
                    "A1\t101\t200\tJM22A1.g1\t0\t+\n"
                    "A1\t201\t300\tJM22Ctg52.g2\t0\t+\n"
                )

            reference_varieties = pipeline._infer_reference_variety_names_from_bed(
                bed_path, ["JM22"]
            )

        self.assertEqual(reference_varieties, ["RefA"])

    def test_reference_prefix_inference_rejects_bed_without_reference_prefixes(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            bed_path = os.path.join(temp_dir, "query-only.bed")
            with open(bed_path, "w", encoding="utf-8") as handle:
                handle.write(
                    "A1\t1\t100\tJM22A1.g1\t0\t+\n"
                    "A1\t101\t200\tJM22Ctg52.g2\t0\t+\n"
                )

            with self.assertRaises(ValueError) as error:
                pipeline._infer_reference_variety_names_from_bed(bed_path, ["JM22"])

        self.assertIn("reference", str(error.exception).lower())

    def test_append_gene_set_excludes_bed_only_genes(self):
        eligible_gene_set = pipeline._validate_append_gene_set(
            {"RefA.g1", "JM22.g1"},
            {"RefA.g1": [], "JM22.g1": [], "RefA.nonrep": []},
        )

        self.assertEqual(eligible_gene_set, {"RefA.g1", "JM22.g1"})
        with self.assertRaisesRegex(ValueError, "missing from append BED"):
            pipeline._validate_append_gene_set(
                {"RefA.g1", "JM22.g1"}, {"RefA.g1": []}
            )

    def test_append_uses_combined_bed_and_shared_clusters(self):
        expected_parameters = [
            "query_cdna_path",
            "query_gdna_path",
            "all_bed_path",
            "refer_cdna_path",
            "refer_gdna_path",
            "bam_path",
            "variety_name",
            "threads",
            "out_dir",
            "prefix",
        ]
        self.assertEqual(
            list(inspect.signature(pipeline.unit_append).parameters), expected_parameters
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            query_cdna_path = os.path.join(temp_dir, "query.cdna.fasta")
            query_gdna_path = os.path.join(temp_dir, "query.gdna.fasta")
            refer_cdna_path = os.path.join(temp_dir, "reference.cdna.fasta")
            refer_gdna_path = os.path.join(temp_dir, "reference.gdna.fasta")
            all_bed_path = os.path.join(temp_dir, "combined.bed")
            bam_path = os.path.join(temp_dir, "existing.bam")
            out_dir = os.path.join(temp_dir, "out")
            merged_cdna_path = os.path.join(out_dir, "Append_merged.cdna.fasta")
            merged_gdna_path = os.path.join(out_dir, "Append_merged.gdna.fasta")
            for fasta_path in (
                query_cdna_path,
                query_gdna_path,
                refer_cdna_path,
                refer_gdna_path,
            ):
                with open(fasta_path, "w", encoding="utf-8"):
                    pass
            with open(all_bed_path, "w", encoding="utf-8") as handle:
                handle.write(
                    "A1\t1\t100\tRefA.g1\t0\t+\n"
                    "A1\t101\t200\tJM22.g1\t0\t+\n"
                    "A1\t201\t300\tRefA.nonrep\t0\t+\n"
                )
            with open(bam_path, "w", encoding="utf-8"):
                pass

            def get_fasta_len_for_merged_inputs(fasta_path):
                if fasta_path == merged_cdna_path:
                    return {"RefA.g1.1": 100, "JM22.g1.1": 90}
                if fasta_path == merged_gdna_path:
                    return {"RefA.g1": 100, "JM22.g1": 90}
                self.fail(f"unexpected FASTA length lookup: {fasta_path}")

            pre_clusters = [["PRE_MARKER"]]
            last_clusters = [["LAST_MARKER"]]
            with mock.patch.object(pipeline, "concat_fasta_files"), \
                 mock.patch.object(
                     pipeline,
                     "get_fasta_len",
                     side_effect=get_fasta_len_for_merged_inputs,
                 ), \
                 mock.patch.object(
                     pipeline,
                     "filter_bam",
                     return_value=([], {"RefA.g1": 100, "JM22.g1": 90}),
                 ) as filter_bam, \
                 mock.patch.object(
                     pipeline,
                     "_derive_clusters_from_alignment",
                     return_value=(pre_clusters, last_clusters),
                 ) as derive_clusters, \
                 mock.patch.object(
                     pipeline,
                     "_write_cluster_outputs",
                     return_value=("output.gtf", "output.cdna", "output.gdna", "output.bed"),
                 ) as write_outputs, \
                 mock.patch.object(
                     pipeline, "derive_last_clusters_from_pre", create=True
                 ) as derive_last_clusters_from_pre, \
                 mock.patch.object(pipeline, "minimap2_map") as minimap2_map:
                pipeline.unit_append(
                    query_cdna_path=query_cdna_path,
                    query_gdna_path=query_gdna_path,
                    all_bed_path=all_bed_path,
                    refer_cdna_path=refer_cdna_path,
                    refer_gdna_path=refer_gdna_path,
                    bam_path=bam_path,
                    variety_name=["JM22"],
                    threads=1,
                    out_dir=out_dir,
                    prefix="Append",
                )

            derive_clusters.assert_called_once()
            self.assertEqual(derive_clusters.call_args.kwargs["refer_prefixes"], ["RefA"])
            self.assertEqual(derive_clusters.call_args.kwargs["variety_li"], ["JM22"])
            self.assertEqual(
                derive_clusters.call_args.kwargs["eligible_gene_set"],
                {"RefA.g1", "JM22.g1"},
            )
            filter_bam.assert_called_once()
            self.assertEqual(filter_bam.call_args.args[0], bam_path)
            minimap2_map.assert_not_called()
            derive_last_clusters_from_pre.assert_not_called()
            self.assertEqual(write_outputs.call_count, 2)
            self.assertEqual(
                write_outputs.call_args_list[0].args[0], {"PRE_MARKER": ["PRE_MARKER"]}
            )
            self.assertEqual(
                write_outputs.call_args_list[1].args[0], {"LAST_MARKER": ["LAST_MARKER"]}
            )
            with open(
                os.path.join(out_dir, "Append_merged.bed"), encoding="utf-8"
            ) as merged_handle, open(all_bed_path, encoding="utf-8") as input_handle:
                self.assertEqual(merged_handle.read(), input_handle.read())

    def test_append_rejects_missing_bam_before_filtering(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            query_cdna_path = os.path.join(temp_dir, "query.cdna.fasta")
            query_gdna_path = os.path.join(temp_dir, "query.gdna.fasta")
            refer_cdna_path = os.path.join(temp_dir, "reference.cdna.fasta")
            refer_gdna_path = os.path.join(temp_dir, "reference.gdna.fasta")
            all_bed_path = os.path.join(temp_dir, "combined.bed")
            bam_path = os.path.join(temp_dir, "missing.bam")
            for fasta_path in (
                query_cdna_path,
                query_gdna_path,
                refer_cdna_path,
                refer_gdna_path,
            ):
                with open(fasta_path, "w", encoding="utf-8"):
                    pass
            with open(all_bed_path, "w", encoding="utf-8") as handle:
                handle.write("A1\t1\t100\tRefA.g1\t0\t+\n")

            with mock.patch.object(pipeline, "filter_bam") as filter_bam:
                with self.assertRaises(FileNotFoundError) as error:
                    pipeline.unit_append(
                        query_cdna_path=query_cdna_path,
                        query_gdna_path=query_gdna_path,
                        all_bed_path=all_bed_path,
                        refer_cdna_path=refer_cdna_path,
                        refer_gdna_path=refer_gdna_path,
                        bam_path=bam_path,
                        variety_name=["JM22"],
                        threads=1,
                        out_dir=os.path.join(temp_dir, "out"),
                    )

            self.assertIn(bam_path, str(error.exception))
            filter_bam.assert_not_called()
