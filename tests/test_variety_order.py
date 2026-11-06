import itertools
import random

import networkx as nx

from pantrans.assign_cluster_gene import assign_gene_by_chrom
from pantrans.construct_di_graph import derive_last_clusters_from_pre


def test_variety_permutations_preserve_assignment_and_last_clusters():
    rng = random.Random(20261102)
    varieties = ['Abo', 'AMN', 'CS']
    permutations = list(itertools.permutations(varieties))
    for case in range(40):
        genes = [f'{v}.g{i}' for v in varieties for i in range(4)]
        labels = {g: g.split('.')[0] for g in genes}
        bed = {g: [rng.choice(['Chr1A', 'Chr1B', 'contig']), 100 * i, 100 * i + 90]
               for i, g in enumerate(genes)}
        lengths = {g: 100 for g in genes}
        graph = nx.DiGraph()
        graph.add_nodes_from(genes)
        graph.add_edges_from((a, b) for a in genes for b in genes if rng.random() < 0.5)
        expected_assignment = expected_last = None
        for order in permutations:
            order = list(order)
            kwargs = dict(main_chroms=['Chr1A', 'Chr1B'], gene_variety_dic=labels)
            assigned = assign_gene_by_chrom(order, genes, bed, graph=graph, **kwargs)
            last = derive_last_clusters_from_pre([genes], graph, lengths, bed, order, 'CS', **kwargs)
            assert sorted(g for c in last for g in c) == sorted(genes)
            if expected_assignment is None:
                expected_assignment, expected_last = assigned, last
            else:
                assert assigned == expected_assignment, (case, order)
                assert last == expected_last, (case, order)


def test_reference_alias_priority_and_input_list_are_preserved():
    genes = ['A.g1', 'B.g1', 'A.g2', 'B.g2']
    bed = {g: ['Chr1A', i * 100, i * 100 + 90] for i, g in enumerate(genes)}
    graph = nx.complete_graph(genes, create_using=nx.DiGraph)
    kwargs = dict(main_chroms=['Chr1A'], refer_prefixes=['A'], graph=graph,
                  gene_variety_dic={g: g.split('.')[0] for g in genes})
    a = ['B', 'A']; original = a[:]
    assert assign_gene_by_chrom(a, genes, bed, **kwargs) == assign_gene_by_chrom(['A', 'B'], genes, bed, **kwargs)
    assert a == original
