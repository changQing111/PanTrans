import networkx as nx
from .assign_cluster_gene import assign_gene_by_chrom

def calculate_in_degree_in_scc(graph, scc):
    """计算SCC中节点的内部入度（仅来自SCC内部的边）"""
    scc_subgraph = graph.subgraph(scc)
    in_degrees = {node: scc_subgraph.in_degree(node) for node in sorted(scc)}
    
    return in_degrees

def find_max_in_degree_nodes(in_degrees, gene_len_dic, refer_name, refer_prefixes=None):
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
        refer_nodes = [node for node in max_nodes if refer_name in node]
    return sorted(refer_nodes or max_nodes)[0]

def is_bidirectionally_connected(graph, nodeA, nodeB):
    """判断两个节点是否有双向连接"""
    return graph.has_edge(nodeA, nodeB) and graph.has_edge(nodeB, nodeA)

def unit_recursion(G, node_li, gene_len_dic, refer_name, refer_prefixes=None):
    """分组：挑选代表节点并输出与其双向连接的簇（迭代实现）

    Args:
        G: 有向图
        node_li: 节点列表
        gene_len_dic: 基因长度字典
        refer_name: 参考基因组名称
    """
    if not node_li:
        return []

    results = []

    # 确保节点都在图中；无效节点直接作为单节点簇输出
    valid_nodes = sorted(node for node in node_li if node in G)
    invalid_nodes = sorted(node for node in node_li if node not in G)
    results.extend([[node] for node in invalid_nodes])

    pending = list(valid_nodes)
    if not pending:
        return results

    while pending:
        if len(pending) == 1:
            results.append([pending[0]])
            break

        node_degree_dic = calculate_in_degree_in_scc(G, set(pending))
        if not node_degree_dic:
            # 无法计算入度时，将剩余节点全部作为单节点簇输出
            results.extend([[node] for node in sorted(pending)])
            break

        max_node = find_max_in_degree_nodes(node_degree_dic, gene_len_dic, refer_name, refer_prefixes=refer_prefixes)
        if max_node not in pending:
            # 理论上不应发生，作为兜底：将剩余节点全部作为单节点簇输出
            results.extend([[node] for node in sorted(pending)])
            break

        di_node_li = [max_node]
        non_node_li = []
        for node in sorted(pending):
            if node == max_node:
                continue
            if is_bidirectionally_connected(G, max_node, node):
                di_node_li.append(node)
            else:
                non_node_li.append(node)

        results.append(di_node_li)
        pending = sorted(non_node_li)

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


def derive_last_clusters_from_pre(
    pre_clusters, G, gene_len_dic, bed_dic, variety_li, refer_name,
    main_chroms=None, refer_prefixes=None,
):
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
                variety_li, annotated_cluster, bed_dic, main_chroms=main_chroms, refer_prefixes=refer_prefixes
            )
            derived_clusters = []
            for ac in assigned_cluster:
                if len(ac) == 1:
                    derived_clusters.append(ac)
                else:
                    derived_clusters.extend(
                        unit_recursion(G, ac, gene_len_dic, refer_name, refer_prefixes=refer_prefixes)
                    )

        derived_gene_set = _cluster_gene_set(derived_clusters)
        omitted_genes = [gene for gene in cluster if gene not in derived_gene_set]
        last_clusters.extend(
            _attach_missing_cluster_members(derived_clusters, omitted_genes, cluster[0])
        )
    return last_clusters

def assign_sccs(sccs, G, gene_len_dic, bed_dic, variety_li, refer_name, main_chroms=None, refer_prefixes=None):
    """
    Split strongly connected components into final clusters by recursively
    resolving graph structure and assigning genes by chromosome.
    """
    pre_cluster = []
    for scc in sccs:
        li = list(scc)
        if len(li) == 1:
            pre_cluster.append(li)
        else:
            clusters = unit_recursion(G, li, gene_len_dic, refer_name, refer_prefixes=refer_prefixes)
            pre_cluster.extend(clusters)
    last_clusters = derive_last_clusters_from_pre(
        pre_cluster, G, gene_len_dic, bed_dic, variety_li, refer_name,
        main_chroms=main_chroms, refer_prefixes=refer_prefixes
    )
    return pre_cluster, last_clusters

