import random
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

    def successors(self, node):
        return (target for source, target in self.edges if source == node)


def full_scan_attachment(clusters, graph, bed, varieties):
    """Reference behavior before the neighbor-index optimization."""
    for source_index, cluster in enumerate(clusters):
        if len(cluster) != 1 or bed.get(cluster[0], (None,))[0] == "A1":
            continue
        gene = cluster[0]
        candidates = []
        for index, target in enumerate(clusters):
            if index == source_index or not target:
                continue
            linked = [member for member in target
                      if graph.has_edge(gene, member) and graph.has_edge(member, gene)]
            if linked:
                same = any(varieties.get(member) == varieties.get(gene)
                           for member in linked) if varieties else False
                candidates.append((not same, len(target), index))
        if candidates:
            clusters[min(candidates)[2]].append(gene)
            clusters[source_index] = []
    return [cluster for cluster in clusters if cluster]


class IncrementalUnitRecursionTest(unittest.TestCase):
    def test_exact_legacy_equivalence(self):
        from legacy_unit_recursion import legacy_unit_recursion
        rng = random.Random(18543099)
        for case in range(400):
            nodes = [f"{prefix}.g{i}" for i, prefix in enumerate(rng.choices(['A', 'B', 'Ref'], k=25))]
            graph = FakeGraph((a, b) for a in nodes for b in nodes if rng.random() < 0.16)
            lengths = {g: rng.choice([10, 20, 30]) for g in nodes if rng.random() < 0.9}
            bed = {g: [rng.choice(['A1', 'Chr1A', 'ctg'])] for g in nodes if rng.random() < 0.9}
            labels = {g: rng.choice(['A', 'B', 'Ref']) for g in nodes if rng.random() < 0.7}
            rng.shuffle(nodes)
            kwargs = dict(
                bed_dic=bed if case % 4 else None,
                main_chroms=['Chr1A'] if case % 3 else None,
                gene_variety_dic=labels if case % 2 else None,
                refer_prefixes=['Ref'] if case % 5 else None,
            )
            with self.subTest(case=case):
                self.assertEqual(
                    construct_di_graph.unit_recursion(graph, nodes + ['missing'], lengths, 'Ref', **kwargs),
                    legacy_unit_recursion(graph, nodes + ['missing'], lengths, 'Ref', **kwargs),
                )

    def test_empty_singleton_and_contig_chain(self):
        from legacy_unit_recursion import legacy_unit_recursion
        graph = FakeGraph([('a', 'b'), ('b', 'a'), ('b', 'c'), ('c', 'b'), ('c', 'c')])
        for nodes in ([], ['a'], ['missing'], ['c', 'b', 'a']):
            self.assertEqual(
                construct_di_graph.unit_recursion(graph, nodes, {}, '', bed_dic={}),
                legacy_unit_recursion(graph, nodes, {}, '', bed_dic={}),
            )


class IndexedAttachmentTest(unittest.TestCase):
    def test_matches_full_scan_with_random_partitions_and_partial_labels(self):
        rng = random.Random(18535481)
        for case in range(200):
            genes = [f"g{i}" for i in range(30)]
            rng.shuffle(genes)
            clusters = []
            while genes:
                size = rng.choice([1, 1, 1, 2, 4])
                clusters.append(genes[:size])
                del genes[:size]
            nodes = [gene for cluster in clusters for gene in cluster]
            # Include self-loops, one-way edges and neighbors outside clusters.
            graph = FakeGraph((a, b) for a in nodes + ["outside"]
                              for b in nodes + ["outside"] if rng.random() < 0.18)
            bed = {gene: [rng.choice(["A1", "contig"])]
                   for gene in nodes if rng.random() < 0.9}
            varieties = {gene: rng.choice(["A", "B"])
                         for gene in nodes if rng.random() < 0.7} if case % 2 else {}
            with self.subTest(case=case):
                expected = full_scan_attachment([c[:] for c in clusters], graph, bed, varieties)
                actual = construct_di_graph._attach_reciprocal_contig_singletons(
                    [c[:] for c in clusters], graph, bed, None, varieties)
                self.assertEqual(actual, expected)

    def test_tracks_moved_neighbors_and_skips_grown_singletons(self):
        graph = FakeGraph([("a", "b"), ("b", "a"), ("c", "a"), ("a", "c")])
        result = construct_di_graph._attach_reciprocal_contig_singletons(
            [["a"], ["b"], ["c"], ["missing"]], graph, {}, None, {})
        self.assertEqual(result, [["b", "a", "c"], ["missing"]])

    def test_sparse_graph_does_not_probe_unrelated_genes(self):
        class CountingGraph(FakeGraph):
            probes = 0

            def has_edge(self, source, target):
                self.probes += 1
                return super().has_edge(source, target)

        graph = CountingGraph([("ctg", "main"), ("main", "ctg")])
        clusters = [[f"unrelated{i}"] for i in range(1000)] + [["main"], ["ctg"]]
        bed = {cluster[0]: ["A1"] for cluster in clusters[:-1]}
        result = construct_di_graph._attach_reciprocal_contig_singletons(
            clusters, graph, bed, None, {})
        self.assertEqual(result[-1], ["main", "ctg"])
        self.assertEqual(graph.probes, 1)


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
