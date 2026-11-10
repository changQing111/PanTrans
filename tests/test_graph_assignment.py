"""Regression tests for graph-supported last assignment, independent of pre generation."""
import itertools
import random
import unittest
from collections import Counter
import networkx as nx
from pantrans.graph_assignment import maximum_assignment
from pantrans.assign_cluster_gene import assign_gene_by_chrom
from pantrans.construct_di_graph import derive_last_clusters_from_pre

class GraphMatchingTests(unittest.TestCase):
    def test_matching_matches_exhaustive_lexicographic_optimum(self):
        rng=random.Random(4881)
        for _ in range(100):
            n=rng.randint(1,4);m=rng.randint(n,5)
            scores=[[tuple(rng.randint(-4,4) for _ in range(4)) for _ in range(m)] for _ in range(n)]
            def value(cols):return tuple(sum(scores[i][j][k] for i,j in enumerate(cols)) for k in range(4))
            chosen=maximum_assignment(scores)
            self.assertEqual(len(set(chosen)),n)
            self.assertEqual(value(chosen),max(value(c) for c in itertools.permutations(range(m),n)))

    def test_crossed_copy_assignment_avoids_extra_singletons(self):
        nodes=['A1','A2','B1','B2'];graph=nx.DiGraph();graph.add_nodes_from(nodes)
        for a,b in [('B1','A1'),('B1','A2'),('B1','B2'),('A1','B2')]:graph.add_edges_from([(a,b),(b,a)])
        bed={x:['Chr1A',i*100,i*100+90] for i,x in enumerate(nodes)}
        labels={x:x[0] for x in nodes}
        last=derive_last_clusters_from_pre([['B1','A1','A2','B2']],graph,{x:100 for x in nodes},bed,['A','B'],'CS',main_chroms=['Chr1A'],gene_variety_dic=labels)
        self.assertEqual({frozenset(row) for row in last},{frozenset(['A1','B2']),frozenset(['A2','B1'])})

    def test_later_variety_evidence_repairs_earlier_tied_assignment(self):
        genes=['A1','A2','B1','B2','C1','C2'];graph=nx.DiGraph();graph.add_nodes_from(genes)
        pairs=[('A1','B1'),('A1','B2'),('A2','B1'),('A2','B2'),
               ('C1','A1'),('C1','B2'),('C2','A2'),('C2','B1')]
        for x,y in pairs:graph.add_edges_from([(x,y),(y,x)])
        bed={g:['Chr1A',i*100,i*100+90] for i,g in enumerate(genes)}
        groups=assign_gene_by_chrom(['A','B','C'],genes,bed,main_chroms=['Chr1A'],gene_variety_dic={g:g[0] for g in genes},graph=graph)
        self.assertEqual({frozenset(x) for x in groups},{frozenset(['A1','B2','C1']),frozenset(['A2','B1','C2'])})

    def test_pre_hub_is_preserved_despite_longer_incoming_only_candidate(self):
        graph=nx.DiGraph([('H','A'),('A','H'),('H','B'),('B','H'),('B','A')])
        bed={x:['Chr1A',i*100,i*100+90] for i,x in enumerate(['H','A','B'])}
        last=derive_last_clusters_from_pre([['H','A','B']],graph,{'H':100,'A':1000,'B':100},bed,['H','A','B'],'CS',main_chroms=['Chr1A'],gene_variety_dic={x:x for x in bed})
        self.assertEqual(last,[['H','A','B']])

    def test_contig_prefers_its_own_variety_not_last_loop_variety(self):
        bed={'A1':['Chr1A',0,10],'A2':['Chr1A',20,30],'B1':['Chr1A',0,10],'B2':['Chr1A',20,30],'Actg':['contig',0,10]}
        labels={g:g[0] for g in bed}
        graph=nx.DiGraph([('Actg','A1'),('A1','Actg'),('Actg','B2'),('B2','Actg')])
        groups=assign_gene_by_chrom(['A','B'],list(bed),bed,main_chroms=['Chr1A'],gene_variety_dic=labels,graph=graph,graph_assignment=False)
        self.assertTrue(any('Actg' in group and 'A1' in group for group in groups))
        self.assertFalse(any('Actg' in group and 'B2' in group for group in groups))

    def test_input_member_order_does_not_change_matching(self):
        graph=nx.complete_graph(['A1','A2','B1','B2','C1'],create_using=nx.DiGraph)
        bed={g:['Chr1A',i*100,i*100+90] for i,g in enumerate(graph)};labels={g:g[0] for g in graph}
        expected=None;rng=random.Random(487)
        for _ in range(20):
            members=list(graph);rng.shuffle(members)
            groups=assign_gene_by_chrom(['C','B','A'],members,bed,main_chroms=['Chr1A'],gene_variety_dic=labels,graph=graph)
            sig={frozenset(x) for x in groups}
            if expected is None:expected=sig
            self.assertEqual(sig,expected)
            self.assertCountEqual([g for row in groups for g in row],list(graph))
            for group in groups:self.assertEqual(max(Counter(labels[g] for g in group).values()),1)

    def test_missing_original_hub_support_does_not_force_merge(self):
        graph=nx.DiGraph([('H','A'),('A','H'),('A','B')]);graph.add_node('B')
        bed={x:['Chr1A',i*100,i*100+90] for i,x in enumerate(['H','A','B'])}
        last=derive_last_clusters_from_pre([['H','A','B']],graph,{g:100 for g in bed},bed,['H','A','B'],'CS',main_chroms=['Chr1A'],gene_variety_dic={x:x for x in bed})
        self.assertGreater(len(last),1)
        for row in last:
            for g in row[1:]:self.assertTrue(graph.has_edge(row[0],g) and graph.has_edge(g,row[0]))

    def test_same_start_recovery_preserves_genes_without_duplicate_variety(self):
        genes=['A1','Adup','B1'];graph=nx.complete_graph(genes,create_using=nx.DiGraph)
        bed={'A1':['Chr1A',0,100],'Adup':['Chr1A',0,90],'B1':['Chr1A',200,300]}
        last=derive_last_clusters_from_pre([genes],graph,{g:100 for g in genes},bed,['A','B'],'CS',main_chroms=['Chr1A'],gene_variety_dic={'A1':'A','Adup':'A','B1':'B'})
        self.assertCountEqual([g for row in last for g in row],genes)
        for row in last:
            self.assertFalse('A1' in row and 'Adup' in row)
            for g in row[1:]:self.assertTrue(graph.has_edge(row[0],g) and graph.has_edge(g,row[0]))
