from typing import Any
import networkx as nx
from .assign_cluster_gene import assign_gene_by_chrom

def calculate_in_degree_in_scc(graph, scc):
    """计算SCC中节点的内部入度（仅来自SCC内部的边）"""
    scc_subgraph = graph.subgraph(scc)
    in_degrees = {node: scc_subgraph.in_degree(node) for node in scc}
    
    return in_degrees

def find_max_in_degree_nodes(in_degrees, gene_len_dic, refer_name, refer_prefixes=None):
    """找到入度最大的节点，若有多个则选长度最长的，若仍有多个则优先参考集基因。"""
    max_degree = max(in_degrees.values())
    max_nodes = [node for node, degree in in_degrees.items() if degree == max_degree]
    len_li = [gene_len_dic.get(i, 0) for i in max_nodes]
    max_len = max(len_li)
    max_index_li = [index for index, l in enumerate(len_li) if l == max_len]
    refer_prefixes = tuple(refer_prefixes or ())
    refer_index_li = []
    if refer_prefixes:
        refer_index_li = [i for i in max_index_li if max_nodes[i].startswith(refer_prefixes)]
    elif refer_name:
        refer_index_li = [i for i in max_index_li if refer_name in max_nodes[i]]
    if refer_index_li:
        max_node = max_nodes[refer_index_li[0]]
    else:
        max_node = max_nodes[max_index_li[0]]
    return max_node

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
    valid_nodes = [node for node in node_li if node in G]
    invalid_nodes = [node for node in node_li if node not in G]
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
            results.extend([[node] for node in pending])
            break

        max_node = find_max_in_degree_nodes(node_degree_dic, gene_len_dic, refer_name, refer_prefixes=refer_prefixes)
        if max_node not in pending:
            # 理论上不应发生，作为兜底：将剩余节点全部作为单节点簇输出
            results.extend([[node] for node in pending])
            break

        di_node_li = [max_node]
        non_node_li = []
        for node in pending:
            if node == max_node:
                continue
            if is_bidirectionally_connected(G, max_node, node):
                di_node_li.append(node)
            else:
                non_node_li.append(node)

        results.append(di_node_li)
        pending = non_node_li

    return results

def di_graph_from_pair(gene_pair_li):
    G = nx.DiGraph()
    G.add_edges_from(gene_pair_li)
    return G

def construct_di_graph_from_cluster(cluster_li):
    """根据最终簇列表构建双向图.

    参数:
        cluster_li: List[List[str]]，每个子列表是一个簇中的基因 ID。
    """
    G = nx.DiGraph()
    for cluster in cluster_li:
        if not cluster:
            continue
        hub = cluster[0]
        others = cluster[1:]
        for other in others:
            G.add_edge(hub, other)
            G.add_edge(other, hub)
    return G

def write_di_graph(G, out_path):
    nx.write_edgelist(G, out_path, data=False)

def read_di_graph(in_path):
    G = nx.read_edgelist(in_path, create_using=nx.DiGraph())
    return G

def get_conn_comp(G):
    return list(nx.strongly_connected_components(G))

def assign_sccs(sccs, G, gene_len_dic, bed_dic, variety_li, refer_name, main_chroms=None, refer_prefixes=None):
    """
    Split strongly connected components into final clusters by recursively
    resolving graph structure and assigning genes by chromosome.
    """
    pre_cluster = []
    last_clusters = []
    for scc in sccs:
        li = list(scc)
        if len(li) == 1:
            pre_cluster.append(li)
            last_clusters.append(li)
        else:
            clusters = unit_recursion(G, li, gene_len_dic, refer_name, refer_prefixes=refer_prefixes)
            pre_cluster.extend(clusters)
            for cluster in clusters:
                if len(cluster) == 1:
                    last_clusters.append(cluster)
                else:
                    assigned_cluster = assign_gene_by_chrom(
                        variety_li, cluster, bed_dic, main_chroms=main_chroms, refer_prefixes=refer_prefixes
                    )
                    for ac in assigned_cluster:
                        if len(ac) == 1:
                            last_clusters.append(ac)
                        else:
                            last_clusters.extend(
                                unit_recursion(G, ac, gene_len_dic, refer_name, refer_prefixes=refer_prefixes)
                            )
    return pre_cluster, last_clusters

