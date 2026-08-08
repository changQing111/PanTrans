import unittest

from pantrans import construct_di_graph


class FakeGraph:
    def __init__(self, edges):
        self.edges = set(edges)
        self.nodes = {
            node
            for edge in self.edges
            for node in edge
        }

    def __contains__(self, node):
        return node in self.nodes

    def subgraph(self, nodes):
        node_set = set(nodes)
        graph = FakeGraph(
            (source, target)
            for source, target in self.edges
            if source in node_set and target in node_set
        )
        graph.nodes = node_set
        return graph

    def in_degree(self, node):
        return sum(1 for _, target in self.edges if target == node)

    def has_edge(self, source, target):
        return (source, target) in self.edges


class DeriveLastClustersTest(unittest.TestCase):
    def derive_last_clusters(self, *args, **kwargs):
        self.assertTrue(
            hasattr(construct_di_graph, "derive_last_clusters_from_pre"),
            "derive_last_clusters_from_pre is required",
        )
        return construct_di_graph.derive_last_clusters_from_pre(*args, **kwargs)

    def test_preserves_pre_member_without_bed_annotation(self):
        graph = FakeGraph(
            [("Q.g1", "RefA.g1"), ("RefA.g1", "Q.g1")]
        )
        result = self.derive_last_clusters(
            [["RefA.g1", "RefA.old", "Q.g1"]],
            graph,
            {"RefA.g1": 100, "Q.g1": 100},
            {
                "RefA.g1": ["A1", 1, 100],
                "Q.g1": ["A1", 10, 110],
            },
            ["Q"],
            "RefA",
            refer_prefixes=["RefA"],
        )

        self.assertEqual(result, [["RefA.g1", "Q.g1", "RefA.old"]])

    def test_preserves_members_omitted_by_coordinate_assignment(self):
        graph = FakeGraph([])
        cluster = ["Q.long", "Q.same_start", "RefA.hub"]
        result = self.derive_last_clusters(
            [cluster],
            graph,
            {gene: 100 for gene in cluster},
            {
                "Q.long": ["A1", 1, 100],
                "Q.same_start": ["A1", 1, 100],
                "RefA.hub": ["A1", 50, 150],
            },
            ["Q"],
            "RefA",
            refer_prefixes=["RefA"],
        )

        result_genes = {gene for derived in result for gene in derived}
        self.assertEqual(result_genes, set(cluster))


if __name__ == "__main__":
    unittest.main()
