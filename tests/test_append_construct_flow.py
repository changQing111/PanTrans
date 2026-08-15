import os
import inspect
import subprocess
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
        with mock.patch.dict(
            sys.modules,
            {"pantrans.pipeline": pipeline, **dependency_stubs()},
        ):
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


class MergeBamFilesTest(unittest.TestCase):
    def test_empty_input_returns_false_without_creating_destination(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            destination = os.path.join(temp_dir, "merged.bam")

            with mock.patch.object(pipeline.subprocess, "run") as run_mock:
                result = pipeline._merge_bam_files([], destination)

            self.assertFalse(result)
            self.assertFalse(os.path.exists(destination))
            run_mock.assert_not_called()

    def test_single_input_is_copied_verbatim(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            source = os.path.join(temp_dir, "input.bam")
            destination = os.path.join(temp_dir, "merged.bam")
            expected_contents = b"BAM\x01\x00test contents"
            with open(source, "wb") as handle:
                handle.write(expected_contents)

            with mock.patch.object(pipeline.subprocess, "run") as run_mock:
                result = pipeline._merge_bam_files([source], destination)

            self.assertTrue(result)
            with open(destination, "rb") as handle:
                self.assertEqual(handle.read(), expected_contents)
            run_mock.assert_not_called()

    def test_129_inputs_are_merged_in_two_bounded_rounds(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            bam_paths = []
            for index in range(129):
                bam_path = os.path.join(temp_dir, f"input-{index:03d}.bam")
                with open(bam_path, "wb") as handle:
                    handle.write(str(index).encode("ascii"))
                bam_paths.append(bam_path)

            destination = os.path.join(temp_dir, "merged.bam")

            def create_merge_output(command, check):
                self.assertTrue(check)
                with open(command[3], "wb") as handle:
                    handle.write(b"merged")

            with mock.patch.object(
                pipeline.subprocess, "run", side_effect=create_merge_output
            ) as run_mock, mock.patch.object(
                pipeline.shutil,
                "copyfile",
                wraps=pipeline.shutil.copyfile,
            ) as copy_mock:
                result = pipeline._merge_bam_files(bam_paths, destination)

            commands = [call.args[0] for call in run_mock.call_args_list]
            copy_mock.assert_called_once()
            copied_source, copied_singleton = copy_mock.call_args.args
            self.assertTrue(result)
            self.assertTrue(os.path.exists(destination))
            self.assertEqual(len(commands), 2)
            self.assertTrue(all(len(command) - 4 <= 128 for command in commands))
            self.assertEqual(len(commands[0]) - 4, 128)
            self.assertEqual(len(commands[1]) - 4, 2)
            self.assertEqual(copied_source, bam_paths[-1])
            self.assertCountEqual(
                commands[1][4:],
                [commands[0][3], copied_singleton],
            )
            self.assertNotEqual(commands[1][3], destination)
            self.assertFalse(
                any(name.startswith("merged.merge-") for name in os.listdir(temp_dir))
            )

    def test_called_process_error_propagates_and_cleans_temporary_directory(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            bam_paths = []
            for index in range(2):
                bam_path = os.path.join(temp_dir, f"input-{index}.bam")
                with open(bam_path, "wb") as handle:
                    handle.write(str(index).encode("ascii"))
                bam_paths.append(bam_path)
            destination = os.path.join(temp_dir, "merged.bam")
            requested_outputs = []

            def fail_merge(command, check):
                self.assertTrue(check)
                requested_outputs.append(command[3])
                raise subprocess.CalledProcessError(1, command)

            with mock.patch.object(
                pipeline.subprocess, "run", side_effect=fail_merge
            ), self.assertRaises(subprocess.CalledProcessError):
                pipeline._merge_bam_files(bam_paths, destination)

            self.assertEqual(len(requested_outputs), 1)
            merge_directory = os.path.dirname(requested_outputs[0])
            self.assertNotEqual(merge_directory, temp_dir)
            self.assertTrue(
                os.path.basename(merge_directory).startswith("merged.merge-")
            )
            self.assertFalse(os.path.exists(merge_directory))
            self.assertFalse(os.path.exists(destination))


class AppendCliContractTest(unittest.TestCase):
    def test_append_accepts_history_graph_and_combined_bed(self):
        args = pantrans_main.read_parameters(
            [
                "append",
                "--name", "JM22",
                "--cdna", "history.unrenamed.cdna.fasta",
                "--history-gtf", "history.unrenamed.gtf",
                "--query-cdna", "query.cdna.fasta",
                "--gdna", "query.gdna.fasta",
                "--bed", "combined.bed",
                "--history-graph", "history.graph.json",
                "--output", "append_out",
            ]
        )

        self.assertEqual(args.bed, "combined.bed")
        self.assertEqual(args.cdna, "history.unrenamed.cdna.fasta")
        self.assertEqual(args.history_gtf, "history.unrenamed.gtf")
        self.assertEqual(args.query_cdna, "query.cdna.fasta")
        self.assertEqual(args.history_graph, "history.graph.json")
        self.assertIsNone(args.query_to_all_bam)
        self.assertIsNone(args.history_to_query_bam)
        self.assertFalse(hasattr(args, "refer_bed"))
        self.assertFalse(hasattr(args, "refer_cluster"))
        self.assertFalse(hasattr(args, "refer_cdna"))
        self.assertFalse(hasattr(args, "refer_gdna"))

    def test_append_requires_history_gtf_and_query_cdna(self):
        base_args = [
            "append",
            "--name", "JM22",
            "--cdna", "history.unrenamed.cdna.fasta",
            "--history-gtf", "history.unrenamed.gtf",
            "--query-cdna", "query.cdna.fasta",
            "--gdna", "query.gdna.fasta",
            "--bed", "combined.bed",
            "--history-graph", "history.graph.json",
            "--output", "append_out",
        ]
        for option in ("--history-gtf", "--query-cdna"):
            with self.subTest(option=option), self.assertRaises(SystemExit):
                option_index = base_args.index(option)
                pantrans_main.read_parameters(
                    base_args[:option_index] + base_args[option_index + 2 :]
                )

    def test_append_rejects_stale_representative_sequence_options(self):
        with self.assertRaises(SystemExit):
            pantrans_main.read_parameters(
                [
                    "append",
                    "--name", "JM22",
                    "--cdna", "history.unrenamed.cdna.fasta",
                    "--history-gtf", "history.unrenamed.gtf",
                    "--query-cdna", "query.cdna.fasta",
                    "--gdna", "query.gdna.fasta",
                    "--bed", "combined.bed",
                    "--history-graph", "history.graph.json",
                    "--refer_cdna", "obsolete.cdna.fasta",
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
    def test_write_cluster_outputs_creates_unrenamed_and_official_transcriptomes(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            rename_map = {"Ref.g1": "Pan1A000001"}
            history_gtf = os.path.join(temp_dir, "history.gtf")

            with mock.patch.object(pipeline, "transcript_dedup") as dedup_mock, \
                 mock.patch.object(
                     pipeline, "_rescue_missing_cluster_genes"
                 ) as rescue_mock, \
                 mock.patch.object(pipeline, "rename_gtf_ids") as rename_mock, \
                 mock.patch.object(pipeline, "sort_gtf_by_gene_id") as sort_mock, \
                 mock.patch.object(
                     pipeline, "extract_fasta_subset_by_names"
                 ) as gdna_mock, \
                 mock.patch.object(pipeline, "get_cdna_from_gtf") as cdna_mock, \
                 mock.patch.object(pipeline, "get_subset_bed") as bed_mock:
                result = pipeline._write_cluster_outputs(
                    cluster_dic={"Ref.g1": ["Ref.g1", "JM22.g1"]},
                    all_gdna_path="merged.gdna.fasta",
                    all_bed_path="merged.bed",
                    filter_bam_path="merged.filtered.bam",
                    trans_len_dic={"Ref.g1.1": 20},
                    gene_len_dic={"Ref.g1": 100},
                    gene_strand_dic={"Ref.g1": "+"},
                    rename_map=rename_map,
                    out_dir=temp_dir,
                    prefix="Pan",
                    label="",
                    all_cdna_path="merged.cdna.fasta",
                    threads=4,
                    enable_rescue=True,
                    pre_gtf_path="pre.gtf",
                    seed_gtf_path=history_gtf,
                )

            official_gtf = os.path.join(temp_dir, "Pan.gtf")
            unrenamed_gtf = os.path.join(temp_dir, "Pan_unrenamed.gtf")
            official_cdna = os.path.join(temp_dir, "Pan_cdna.refer.fasta")
            unrenamed_cdna = os.path.join(
                temp_dir, "Pan_unrenamed_cdna.refer.fasta"
            )

            dedup_mock.assert_called_once_with(
                "merged.filtered.bam",
                cluster_dic={"Ref.g1": ["Ref.g1", "JM22.g1"]},
                trans_len_dic={"Ref.g1.1": 20},
                gene_len_dic={"Ref.g1": 100},
                gene_strand_dic={"Ref.g1": "+"},
                rename_map=None,
                gtf_path=unrenamed_gtf,
                seed_gtf_path=history_gtf,
            )
            self.assertEqual(rescue_mock.call_args.kwargs["gtf_path"], unrenamed_gtf)
            self.assertIsNone(rescue_mock.call_args.kwargs["rename_map"])
            rename_mock.assert_called_once_with(
                unrenamed_gtf, official_gtf, rename_map
            )
            sort_mock.assert_called_once_with(official_gtf)
            self.assertEqual(
                cdna_mock.call_args_list,
                [
                    mock.call("merged.gdna.fasta", unrenamed_gtf, unrenamed_cdna),
                    mock.call("merged.gdna.fasta", official_gtf, official_cdna),
                ],
            )
            gdna_mock.assert_called_once()
            bed_mock.assert_called_once()
            self.assertEqual(
                result,
                (
                    official_gtf,
                    official_cdna,
                    os.path.join(temp_dir, "Pan_gdna.refer.fasta"),
                    os.path.join(temp_dir, "Pan.refer.bed"),
                ),
            )

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
            ) as write_outputs, mock.patch.object(
                pipeline,
                "bam_alignment_provenance",
                return_value={
                    "minimap2_version": "test",
                    "minimap2_options": "test-options",
                    "filter_thresholds": {},
                },
            ), mock.patch.object(
                pipeline,
                "write_graph_package",
            ):
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
