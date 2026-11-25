"""Frozen reference for exact clustering compatibility tests."""
from pantrans.construct_di_graph import (calculate_in_degree_in_scc, find_max_in_degree_nodes, is_bidirectionally_connected)
from pantrans.assign_cluster_gene import MAIN_CHROMS

def legacy_unit_recursion(
    G, node_li, gene_len_dic, refer_name, refer_prefixes=None, gene_variety_dic=None,
    bed_dic=None, main_chroms=None,
):
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

        eligible_representatives = pending
        if len(pending) > 1 and bed_dic is not None:
            main_chrom_set = set(main_chroms or MAIN_CHROMS)
            if main_chroms:
                from pantrans.assign_cluster_gene import build_main_chrom_order
                main_chrom_set = set(build_main_chrom_order(main_chroms))
            eligible_representatives = [
                node for node in pending
                if bed_dic.get(node, (None,))[0] in main_chrom_set
            ]
            if not eligible_representatives:
                for node in pending:
                    compatible = [
                        index
                        for index, cluster in enumerate(results)
                        if cluster and any(is_bidirectionally_connected(G, node, member) for member in cluster)
                    ]
                    if compatible:
                        results[compatible[0]].append(node)
                    else:
                        results.append([node])
                break
        node_degree_dic = calculate_in_degree_in_scc(G, set(pending))
        if len(pending) > 1 and bed_dic is not None:
            node_degree_dic = {
                node: degree for node, degree in node_degree_dic.items()
                if node in eligible_representatives
            }
        if not node_degree_dic:
            # 无法计算入度时，将剩余节点全部作为单节点簇输出
            results.extend([[node] for node in sorted(pending)])
            break

        max_node = find_max_in_degree_nodes(
            node_degree_dic,
            gene_len_dic,
            refer_name,
            refer_prefixes=refer_prefixes,
            gene_variety_dic=gene_variety_dic,
        )
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
