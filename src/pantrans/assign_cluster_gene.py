MAIN_CHROMS = tuple(
    [f"A{i}" for i in range(1, 8)] +
    [f"B{i}" for i in range(1, 8)] +
    [f"D{i}" for i in range(1, 8)]
)

def build_main_chrom_order(main_chroms=None):
    chroms = tuple(main_chroms) if main_chroms else MAIN_CHROMS
    return {chrom: index for index, chrom in enumerate(chroms)}

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

def get_variety_n(query, var_li, refer_prefixes=None):
    """严格匹配品种名（前缀匹配），可将参考集前缀折叠为 Refer。"""
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

def assign_gene_by_chrom(variety_li, cluster_li, bed_dic, main_chroms=None, refer_prefixes=None):
    """按“品种覆盖优先、主染色体优先、contig 灵活补位”策略输出簇列表.
    
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
        variety = get_variety_n(gene, normalized_variety_li, refer_prefixes=refer_prefixes)
        if variety is None:
            unknown_genes.append(gene)
            continue
        variety_gene_dic[variety].append(gene)

    for variety, genes in variety_gene_dic.items():
        if genes:
            variety_gene_dic[variety] = sort_gene(genes, bed_dic, main_chrom_order)

    max_group_num = max((len(genes) for genes in variety_gene_dic.values()), default=0)
    if max_group_num <= 1 and not unknown_genes:
        return [cluster_li]

    if max_group_num == 0:
        return [[gene] for gene in unknown_genes]

    seed_variety = select_seed_variety(normalized_variety_li, variety_gene_dic, bed_dic, main_chrom_order)
    groups = [{"genes": [], "main_chrom": None} for _ in range(max_group_num)]

    for index, gene in enumerate(variety_gene_dic[seed_variety]):
        groups[index]["genes"].append(gene)
        chrom = bed_dic[gene][0]
        if is_main_chrom(chrom, main_chrom_order):
            groups[index]["main_chrom"] = chrom

    for variety in normalized_variety_li:
        if variety == seed_variety or not variety_gene_dic[variety]:
            continue
        assign_variety_genes(groups, variety_gene_dic[variety], bed_dic, main_chrom_order)

    for gene in sort_gene(unknown_genes, bed_dic, main_chrom_order) if unknown_genes else []:
        chosen_index = min(range(len(groups)), key=lambda i: (len(groups[i]["genes"]), i))
        groups[chosen_index]["genes"].append(gene)

    return [group["genes"] for group in groups if group["genes"]]

