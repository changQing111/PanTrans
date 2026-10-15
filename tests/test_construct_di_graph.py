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

    def test_ctg_cannot_be_scc_representative_when_main_gene_exists(self):
        graph = FakeGraph([
            ("Ctg.g1", "Main.g1"), ("Main.g1", "Ctg.g1"),
            ("Ctg.g1", "Main.g2"), ("Main.g2", "Ctg.g1"),
        ])
        result = construct_di_graph.unit_recursion(
            graph,
            ["Ctg.g1", "Main.g1", "Main.g2"],
            {"Ctg.g1": 1000, "Main.g1": 10, "Main.g2": 10},
            "Main",
            bed_dic={
                "Ctg.g1": ["contig", 1, 100],
                "Main.g1": ["Chr1A", 1, 100],
                "Main.g2": ["Chr1A", 200, 300],
            },
            main_chroms=["Chr1A"],
        )
        self.assertEqual(result[0][0], "Main.g1")
        self.assertIn("Ctg.g1", result[0])

    def test_ctg_is_attached_only_with_reciprocal_edge(self):
        from pantrans.assign_cluster_gene import assign_gene_by_chrom

        bed = {
            "A.main": ["Chr1A", 1, 100],
            "B.main": ["Chr1A", 1, 100],
            "A.ctg": ["contig", 1, 100],
        }
        reciprocal = FakeGraph([("A.ctg", "A.main"), ("A.main", "A.ctg")])
        assigned = assign_gene_by_chrom(
            ["A", "B"], ["A.main", "B.main", "A.ctg"], bed,
            main_chroms=["Chr1A"],
            gene_variety_dic={"A.main": "A", "B.main": "B", "A.ctg": "A"},
            graph=reciprocal,
        )
        self.assertTrue(any({"A.main", "A.ctg"}.issubset(group) for group in assigned))

        one_way = FakeGraph([("A.ctg", "A.main")])
        assigned = assign_gene_by_chrom(
            ["A", "B"], ["A.main", "B.main", "A.ctg"], bed,
            main_chroms=["Chr1A"],
            gene_variety_dic={"A.main": "A", "B.main": "B", "A.ctg": "A"},
            graph=one_way,
        )
        self.assertTrue(any(group == ["A.ctg"] for group in assigned))

    def test_last_clusters_absorb_ctg_with_direct_reciprocal_homolog(self):
        graph = FakeGraph([("A.ctg", "A.main"), ("A.main", "A.ctg")])
        result = self.derive_last_clusters(
            [["A.main", "A.ctg"]], graph,
            {"A.main": 100, "A.ctg": 1000},
            {"A.main": ["Chr1A", 1, 100], "A.ctg": ["contig", 1, 100]},
            ["A"], "A", main_chroms=["Chr1A"],
            gene_variety_dic={"A.main": "A", "A.ctg": "A"},
        )
        self.assertEqual(result, [["A.main", "A.ctg"]])

    def test_last_clusters_keep_unrelated_ctg_as_singleton(self):
        graph = FakeGraph([("A.ctg", "A.main")])
        result = self.derive_last_clusters(
            [["A.main", "A.ctg"]], graph,
            {"A.main": 100, "A.ctg": 1000},
            {"A.main": ["Chr1A", 1, 100], "A.ctg": ["contig", 1, 100]},
            ["A"], "A", main_chroms=["Chr1A"],
            gene_variety_dic={"A.main": "A", "A.ctg": "A"},
        )
        self.assertCountEqual(result, [["A.main"], ["A.ctg"]])


if __name__ == "__main__":
    unittest.main()
