import os
import sys
import tempfile
import types
import unittest
from unittest import mock


ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")


def install_dependency_stubs():
    sys.modules["pysam"] = types.ModuleType("pysam")

    networkx_stub = types.ModuleType("networkx")
    networkx_stub.DiGraph = lambda: None
    networkx_stub.strongly_connected_components = lambda graph: []
    sys.modules["networkx"] = networkx_stub

    bio_stub = types.ModuleType("Bio")
    bio_stub.__path__ = []
    sys.modules["Bio"] = bio_stub
    for name in ("Bio.SeqIO", "Bio.Seq", "Bio.SeqRecord"):
        sys.modules[name] = types.ModuleType(name)
    sys.modules["Bio.Seq"].Seq = str
    sys.modules["Bio.SeqRecord"].SeqRecord = object


install_dependency_stubs()
sys.path.insert(0, SRC)

from pantrans import main as pantrans_main
from pantrans import pipeline


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

