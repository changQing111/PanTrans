import logging
from time import perf_counter

import networkx as nx
from .assign_cluster_gene import MAIN_CHROMS, assign_gene_by_chrom

logger = logging.getLogger(__name__)

def calculate_in_degree_in_scc(graph, scc):
    """计算SCC中节点的内部入度（仅来自SCC内部的边）"""
    scc_subgraph = graph.subgraph(scc)
    in_degrees = {node: scc_subgraph.in_degree(node) for node in sorted(scc)}
    
    return in_degrees

def find_max_in_degree_nodes(
    in_degrees, gene_len_dic, refer_name, refer_prefixes=None, gene_variety_dic=None
):
    """找到入度最大的节点，若有多个则选长度最长的，若仍有多个则优先参考集基因。"""
    max_degree = max(in_degrees.values())
    max_nodes = sorted(node for node, degree in in_degrees.items() if degree == max_degree)
    len_li = [gene_len_dic.get(i, 0) for i in max_nodes]
    max_len = max(len_li)
    max_nodes = [node for node, length in zip(max_nodes, len_li) if length == max_len]
    refer_prefixes = tuple(refer_prefixes or ())
    refer_nodes = []
    if refer_prefixes:
        refer_nodes = [node for node in max_nodes if node.startswith(refer_prefixes)]
    elif refer_name:
        if gene_variety_dic is not None:
            refer_nodes = [
                node for node in max_nodes if gene_variety_dic.get(node) == refer_name
            ]
            has_labeled_candidates = any(node in gene_variety_dic for node in max_nodes)
        else:
            has_labeled_candidates = False
        if not refer_nodes and not has_labeled_candidates:
            refer_nodes = [node for node in max_nodes if refer_name in node]
    return sorted(refer_nodes or max_nodes)[0]

def is_bidirectionally_connected(graph, nodeA, nodeB):
    """判断两个节点是否有双向连接"""
    return graph.has_edge(nodeA, nodeB) and graph.has_edge(nodeB, nodeA)


def unit_recursion(
    G, node_li, gene_len_dic, refer_name, refer_prefixes=None, gene_variety_dic=None,
    bed_dic=None, main_chroms=None,
):
    """Split a component, maintaining induced in-degrees as genes are removed.

    Representative tie-breaking and output ordering match the full-scan version.
    Input nodes are unique gene IDs, as supplied by SCC and chromosome grouping.
    """
    results = [[node] for node in sorted(node for node in node_li if node not in G)]
    pending = {node for node in node_li if node in G}
    if not pending:
        return results
    if len(pending) == 1:
        results.append([next(iter(pending))])
        return results

    if bed_dic is None:
        eligible = pending.copy()
    else:
        from .assign_cluster_gene import build_main_chrom_order
        main_chrom_set = set(build_main_chrom_order(main_chroms))
        eligible = {node for node in pending
                    if bed_dic.get(node, (None,))[0] in main_chrom_set}

    # Count every edge within the initial component once. Removing a source
    # changes only the in-degrees of its still-pending successors.
    degrees = {node: 0 for node in pending}
    if eligible:
        for node in pending:
            for target in G.successors(node):
                if target in pending:
                    degrees[target] += 1

    while pending:
        if len(pending) == 1:
            results.append([next(iter(pending))])
            break
        if not eligible:
            # Preserve the first-compatible-cluster rule, including attachments
            # to contigs processed earlier in this loop.
            gene_to_cluster = {gene: i for i, cluster in enumerate(results) for gene in cluster}
            for node in sorted(pending):
                targets = [gene_to_cluster[member] for member in G.successors(node)
                           if member in gene_to_cluster and G.has_edge(member, node)]
                if targets:
                    index = min(targets)
                    results[index].append(node)
                else:
                    index = len(results)
                    results.append([node])
                gene_to_cluster[node] = index
            break

        max_node = find_max_in_degree_nodes(
            {node: degrees[node] for node in eligible}, gene_len_dic, refer_name,
            refer_prefixes=refer_prefixes, gene_variety_dic=gene_variety_dic,
        )
        members = sorted(node for node in G.successors(max_node)
                         if node != max_node and node in pending and G.has_edge(node, max_node))
        cluster = [max_node] + members
        results.append(cluster)
        pending.difference_update(cluster)
        eligible.difference_update(cluster)
        for node in cluster:
            for target in G.successors(node):
                if target in pending:
                    degrees[target] -= 1
    return results


def di_graph_from_pair(gene_pair_li):
    G = nx.DiGraph()
    G.add_edges_from(gene_pair_li)
    return G

def get_conn_comp(G):
    return sorted(
        (tuple(sorted(component)) for component in nx.strongly_connected_components(G)),
        key=lambda component: component
    )


def _attach_missing_cluster_members(derived_clusters, genes_to_attach, preferred_gene):
    if not genes_to_attach:
        return derived_clusters
    if not derived_clusters:
        return [genes_to_attach[:]]

    target_index = 0
    for index, cluster in enumerate(derived_clusters):
        if preferred_gene in cluster:
            target_index = index
            break

    seen = set(derived_clusters[target_index])
    for gene in genes_to_attach:
        if gene not in seen:
            derived_clusters[target_index].append(gene)
            seen.add(gene)
    return derived_clusters


def _cluster_gene_set(clusters):
    return {gene for cluster in clusters for gene in cluster}


def _attach_reciprocal_contig_singletons(clusters, graph, bed_dic, main_chroms, gene_variety_dic):
    """Attach singleton unplaced genes to any cluster with direct reciprocal support."""
    started = perf_counter()
    logger.info("Start reciprocal contig singleton attachment: %d clusters", len(clusters))
    main_chrom_set = set(main_chroms or MAIN_CHROMS)
    if main_chroms:
        from .assign_cluster_gene import build_main_chrom_order
        main_chrom_set = set(build_main_chrom_order(main_chroms))
    # Clusters partition the genes. Keep stable cluster indexes while moving
    # genes so later singletons see the same state as the sequential full scan.
    gene_to_cluster = {
        gene: index for index, cluster in enumerate(clusters) for gene in cluster
    }
    examined = attached = 0
    for source_index, source_cluster in enumerate(clusters):
        if len(source_cluster) != 1:
            continue
        gene = source_cluster[0]
        if bed_dic.get(gene, (None,))[0] in main_chrom_set:
            continue
        examined += 1
        if gene not in graph:
            continue
        same_variety_by_cluster = {}
        for member in graph.successors(gene):
            if not graph.has_edge(member, gene):
                continue
            target_index = gene_to_cluster.get(member)
            if target_index is None or target_index == source_index:
                continue
            same_variety = (
                gene_variety_dic.get(member) == gene_variety_dic.get(gene)
            ) if gene_variety_dic else False
            same_variety_by_cluster[target_index] = (
                same_variety_by_cluster.get(target_index, False) or same_variety
            )
        if same_variety_by_cluster:
            target_index = min(
                (not same_variety, len(clusters[index]), index)
                for index, same_variety in same_variety_by_cluster.items()
            )[2]
            clusters[target_index].append(gene)
            clusters[source_index] = []
            gene_to_cluster[gene] = target_index
            attached += 1
    result = [cluster for cluster in clusters if cluster]
    logger.info(
        "Finish reciprocal contig singleton attachment in %.3f s: "
        "examined %d, attached %d, remaining clusters %d",
        perf_counter() - started, examined, attached, len(result),
    )
    return result


def derive_last_clusters_from_pre(
    pre_clusters, G, gene_len_dic, bed_dic, variety_li, refer_name,
    main_chroms=None, refer_prefixes=None, gene_variety_dic=None,
):
    started = perf_counter()
    logger.info("Start chromosome assignment and last-cluster splitting")
    last_clusters = []
    for cluster in pre_clusters:
        cluster = list(cluster)
        annotated_cluster = [gene for gene in cluster if gene in bed_dic]
        if not annotated_cluster:
            last_clusters.append(cluster)
            continue
        if len(annotated_cluster) == 1:
            derived_clusters = [annotated_cluster]
        else:
            assigned_cluster = assign_gene_by_chrom(
                variety_li,
                annotated_cluster,
                bed_dic,
                main_chroms=main_chroms,
                refer_prefixes=refer_prefixes,
                gene_variety_dic=gene_variety_dic,
                graph=G,
            )
            derived_clusters = []
            for ac in assigned_cluster:
                if len(ac) == 1:
                    derived_clusters.append(ac)
                else:
                    derived_clusters.extend(
                        unit_recursion(
                            G,
                            ac,
                            gene_len_dic,
                            refer_name,
                            refer_prefixes=refer_prefixes,
                            gene_variety_dic=gene_variety_dic,
                            bed_dic=bed_dic,
                            main_chroms=main_chroms,
                        )
                    )

        derived_gene_set = _cluster_gene_set(derived_clusters)
        omitted_genes = [gene for gene in cluster if gene not in derived_gene_set]
        last_clusters.extend(
            _attach_missing_cluster_members(derived_clusters, omitted_genes, cluster[0])
        )
    logger.info(
        "Finish chromosome assignment and last-cluster splitting in %.3f s: %d clusters",
        perf_counter() - started, len(last_clusters),
    )
    return _attach_reciprocal_contig_singletons(
        last_clusters, G, bed_dic, main_chroms, gene_variety_dic or {}
    )

def assign_sccs(sccs, G, gene_len_dic, bed_dic, variety_li, refer_name, main_chroms=None, refer_prefixes=None, gene_variety_dic=None):
    """
    Split strongly connected components into final clusters by recursively
    resolving graph structure and assigning genes by chromosome.
    """
    started = perf_counter()
    logger.info("Start SCC splitting into pre-clusters")
    pre_cluster = []
    for scc in sccs:
        li = list(scc)
        if len(li) == 1:
            pre_cluster.append(li)
        else:
            clusters = unit_recursion(
                G,
                li,
                gene_len_dic,
                refer_name,
                refer_prefixes=refer_prefixes,
                gene_variety_dic=gene_variety_dic,
                bed_dic=bed_dic,
                main_chroms=main_chroms,
            )
            pre_cluster.extend(clusters)
    logger.info(
        "Finish SCC splitting into pre-clusters in %.3f s: %d clusters",
        perf_counter() - started, len(pre_cluster),
    )
    last_clusters = derive_last_clusters_from_pre(
        pre_cluster, G, gene_len_dic, bed_dic, variety_li, refer_name,
        main_chroms=main_chroms,
        refer_prefixes=refer_prefixes,
        gene_variety_dic=gene_variety_dic,
    )
    return pre_cluster, last_clusters

