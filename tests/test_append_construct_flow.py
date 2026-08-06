import os
import sys
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
    try:
        with mock.patch.dict(sys.modules, dependency_stubs()):
            from pantrans import pipeline as imported_pipeline
    finally:
        # Retain the direct reference below, but make later test modules import
        # PanTrans normally instead of reusing modules loaded with fake globals.
        for module_name in tuple(sys.modules):
            if module_name == "pantrans" or module_name.startswith("pantrans."):
                del sys.modules[module_name]
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
