import re

MAIN_CHROMS = tuple(
    [f"A{i}" for i in range(1, 8)] +
    [f"B{i}" for i in range(1, 8)] +
    [f"D{i}" for i in range(1, 8)]
)

def build_main_chrom_order(main_chroms=None):
    chroms = tuple(main_chroms) if main_chroms else MAIN_CHROMS
    order = {}
    for index, chrom in enumerate(chroms):
        order[chrom] = index
        match = re.match(r"^(?:Chr)?([1-7])([ABD])$", chrom)
        if match:
            number, subgenome = match.groups()
            for alias in (f"{subgenome}{number}", f"{number}{subgenome}", f"Chr{number}{subgenome}"):
                order.setdefault(alias, index)
    return order

def is_main_chrom(chrom, main_chrom_order):
    return chrom in main_chrom_order

def chrom_sort_key(chrom, main_chrom_order):
    if chrom in main_chrom_order:
        return (0, main_chrom_order[chrom], chrom)
    return (1, chrom.lower(), chrom)

def merge_same_start_gene(gene_li, bed_dic):
    """同一染色体起点的多个基因，取长度最长的"""
    gene_len_dic = {i: (int(bed_dic[i][2]) - int(bed_dic[i][1]) + 1) for i in gene_li}
    gene_coord_dic = {}
    for i in gene_li:
        coord = (bed_dic[i][0], bed_dic[i][1])
        gene_coord_dic.setdefault(coord, []).append(i)
    filter_gene_li = []
    for coord, genes in gene_coord_dic.items():
        lengths = [gene_len_dic[g] for g in genes]
        max_index = lengths.index(max(lengths))
        filter_gene_li.append(genes[max_index])
    return filter_gene_li

def sort_gene(gene_li, bed_dic, main_chrom_order):
    """按染色体和起始坐标排序基因"""
    filter_gene_li = merge_same_start_gene(gene_li, bed_dic)
    if len(filter_gene_li) == 1:
        return filter_gene_li
    else:
        sorted_li = sorted(
            filter_gene_li,
            key=lambda g: (
                chrom_sort_key(bed_dic[g][0], main_chrom_order),
                int(bed_dic[g][1]),
                int(bed_dic[g][2]),
                g,
            ),
        )
        return sorted_li

def get_variety_n(query, var_li, refer_prefixes=None, gene_variety_dic=None):
    """严格匹配品种名（前缀匹配），可将参考集前缀折叠为 Refer。"""
    if gene_variety_dic is not None:
        if query not in gene_variety_dic:
            return get_variety_n(query, var_li, refer_prefixes=refer_prefixes)
        mapped = gene_variety_dic[query]
        if refer_prefixes and mapped in refer_prefixes:
            return "Refer"
        return mapped if mapped in var_li else None
    refer_prefixes = tuple(refer_prefixes or [])
    if refer_prefixes and query.startswith(refer_prefixes):
        return "Refer"
    for v in var_li:
        if v == "Refer":
            continue
        if query.startswith(v):
            return v
    return None

def select_seed_variety(variety_li, variety_gene_dic, bed_dic, main_chrom_order):
    present_varieties = [v for v in variety_li if variety_gene_dic[v]]
    return max(
        present_varieties,
        key=lambda v: (
            len(variety_gene_dic[v]),
            sum(1 for g in variety_gene_dic[v] if is_main_chrom(bed_dic[g][0], main_chrom_order)),
            -variety_li.index(v),
        ),
    )

def choose_group_for_gene(groups, candidate_indexes, gene, bed_dic, main_chrom_order):
    gene_chrom = bed_dic[gene][0]
    gene_is_main = is_main_chrom(gene_chrom, main_chrom_order)

    def score(index):
        group = groups[index]
        same_main = 1 if gene_is_main and group["main_chrom"] == gene_chrom else 0
        can_anchor = 1 if gene_is_main and group["main_chrom"] is None else 0
        contig_flexible = 1 if (not gene_is_main) and group["main_chrom"] is None else 0
        return (
            same_main,
            can_anchor,
            contig_flexible,
            -len(group["genes"]),
            -index,
        )

    return max(candidate_indexes, key=score)

def assign_variety_genes(groups, genes, bed_dic, main_chrom_order):
    available_indexes = list(range(len(groups)))
    main_genes = [g for g in genes if is_main_chrom(bed_dic[g][0], main_chrom_order)]
    contig_genes = [g for g in genes if not is_main_chrom(bed_dic[g][0], main_chrom_order)]

    for gene in main_genes + contig_genes:
        if not available_indexes:
            break
        chosen_index = choose_group_for_gene(groups, available_indexes, gene, bed_dic, main_chrom_order)
        groups[chosen_index]["genes"].append(gene)
        if groups[chosen_index]["main_chrom"] is None and is_main_chrom(bed_dic[gene][0], main_chrom_order):
            groups[chosen_index]["main_chrom"] = bed_dic[gene][0]
        available_indexes.remove(chosen_index)

def assign_gene_by_chrom(variety_li, cluster_li, bed_dic, main_chroms=None, refer_prefixes=None, gene_variety_dic=None, graph=None):
    """Assign main-chromosome genes first, then attach reciprocal Ctg homologs.
    
    参数:
        variety_li: List[str], 品种名称列表
        cluster_li: List[str], 一个基因簇中的基因ID列表（不是嵌套列表）
        bed_dic: Dict[str, Tuple], 基因位置信息字典，key为基因ID

    返回值:
        List[List[str]]: 每个子列表是一组同簇基因。
    """
    main_chrom_order = build_main_chrom_order(main_chroms)
    normalized_variety_li = ["Refer"] + [v for v in variety_li if v != "Refer"] if refer_prefixes else variety_li[:]
    variety_gene_dic = {v: [] for v in normalized_variety_li}
    unknown_genes = []
    for gene in cluster_li:
        variety = get_variety_n(
            gene,
            normalized_variety_li,
            refer_prefixes=refer_prefixes,
            gene_variety_dic=gene_variety_dic,
        )
        if variety is None:
            unknown_genes.append(gene)
            continue
        variety_gene_dic[variety].append(gene)

    contig_genes_by_variety = {variety: [] for variety in normalized_variety_li}
    for variety, genes in variety_gene_dic.items():
        if genes:
            main_genes = [gene for gene in genes if is_main_chrom(bed_dic[gene][0], main_chrom_order)]
            contig_genes_by_variety[variety] = sort_gene(
                [gene for gene in genes if not is_main_chrom(bed_dic[gene][0], main_chrom_order)],
                bed_dic,
                main_chrom_order,
            ) if any(not is_main_chrom(bed_dic[gene][0], main_chrom_order) for gene in genes) else []
            variety_gene_dic[variety] = sort_gene(main_genes, bed_dic, main_chrom_order) if main_genes else []

    max_group_num = max((len(genes) for genes in variety_gene_dic.values()), default=0)
    contig_genes = [gene for variety in normalized_variety_li for gene in contig_genes_by_variety[variety]]
    if max_group_num <= 1 and not contig_genes and not unknown_genes:
        return [cluster_li]

    groups = [{"genes": [], "main_chrom": None} for _ in range(max_group_num)]

    seed_variety = None
    if max_group_num:
        seed_variety = select_seed_variety(normalized_variety_li, variety_gene_dic, bed_dic, main_chrom_order)
        for index, gene in enumerate(variety_gene_dic[seed_variety]):
            groups[index]["genes"].append(gene)
            groups[index]["main_chrom"] = bed_dic[gene][0]

    for variety in normalized_variety_li:
        if variety == seed_variety or not variety_gene_dic[variety] or not groups:
            continue
        assign_variety_genes(groups, variety_gene_dic[variety], bed_dic, main_chrom_order)

    for gene in contig_genes:
        matching_groups = []
        if graph is not None:
            for index, group in enumerate(groups):
                linked_members = [
                    member for member in group["genes"]
                    if graph.has_edge(gene, member) and graph.has_edge(member, gene)
                ]
                if linked_members:
                    for member in linked_members:
                        same_variety = (
                            (gene_variety_dic or {}).get(member) == variety
                            if gene_variety_dic is not None
                            else get_variety_n(
                                member, normalized_variety_li,
                                refer_prefixes=refer_prefixes,
                            ) == variety
                        )
                        matching_groups.append((not same_variety, len(group["genes"]), index))
        if matching_groups:
            chosen_index = min(matching_groups)[2]
            groups[chosen_index]["genes"].append(gene)
        else:
            groups.append({"genes": [gene], "main_chrom": None})

    for gene in sort_gene(unknown_genes, bed_dic, main_chrom_order) if unknown_genes else []:
        if not groups:
            groups.append({"genes": [], "main_chrom": None})
        chosen_index = min(range(len(groups)), key=lambda i: (len(groups[i]["genes"]), i))
        groups[chosen_index]["genes"].append(gene)

    return [group["genes"] for group in groups if group["genes"]]

