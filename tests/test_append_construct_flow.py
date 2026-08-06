import os
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
