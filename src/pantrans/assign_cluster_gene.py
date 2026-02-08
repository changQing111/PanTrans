fill_char = "NA"

def traverse_pseu_mat(li_li):
    """逐列输出矩阵，去掉填充字符"""
    for i in range(len(li_li[0])):
        out_li = [li_li[j][i] for j in range(len(li_li))]
        new_li = [a for a in out_li if a != fill_char]
        if new_li:  # 只有非空时才输出
            return new_li

def fill_NA(li_li):
    """对齐不同长度的列表，用 fill_char 填充"""
    nrow = get_max_len(li_li)
    new_li_li = []
    for li in li_li:
        new_li = li[:]  # 拷贝，避免原地修改
        num = nrow - len(li)
        if num > 0:
            new_li.extend([fill_char] * num)
        new_li_li.append(new_li)
    return new_li_li

def get_max_len(li_li):
    return max(len(li) for li in li_li)

def split_gene_by_chrom(gene_li, bed_dic, chrom_s):
    """按染色体分组基因"""
    chrom_gene_dic = {}
    for i in gene_li:
        chrom = bed_dic[i][0]
        chrom_s.add(chrom)
        chrom_gene_dic.setdefault(chrom, []).append(i)
    return chrom_gene_dic

def merge_same_start_gene(gene_li, bed_dic):
    """同一染色体起点的多个基因，取长度最长的"""
    gene_len_dic = {i: (int(bed_dic[i][2]) - int(bed_dic[i][1]) + 1) for i in gene_li}
    gene_coord_dic = {}
    for i in gene_li:
        coord = bed_dic[i][0] + "_" + bed_dic[i][1]
        gene_coord_dic.setdefault(coord, []).append(i)
    filter_gene_li = []
    for coord, genes in gene_coord_dic.items():
        lengths = [gene_len_dic[g] for g in genes]
        max_index = lengths.index(max(lengths))
        filter_gene_li.append(genes[max_index])
    return filter_gene_li

def sort_gene(gene_li, bed_dic):
    """按染色体和起始坐标排序基因"""
    filter_gene_li = merge_same_start_gene(gene_li, bed_dic)
    if len(filter_gene_li) == 1:
        return filter_gene_li
    else:
        sorted_li = sorted(filter_gene_li, key=lambda g: (bed_dic[g][0], int(bed_dic[g][1])))
        return sorted_li

def assign_gene(sorted_chrom_gene_dic_li, all_chrom_s):
    """对齐不同品种的染色体基因列表并输出成簇列表.

    返回值:
        List[List[str]]: 每个子列表是一组在不同品种上按染色体位置对齐的基因 ID。
    """
    results = []
    for chrom in all_chrom_s:
        li_li = []
        for dic in sorted_chrom_gene_dic_li:
            if chrom in dic:
                li_li.append(dic[chrom])
        if li_li:
            fill_li_li = fill_NA(li_li)
            new_li = traverse_pseu_mat(fill_li_li)
            if new_li:
                results.append(new_li)
    return results

def find_duplicate_gene(variety_set, variety_list, gene_list):
    """找到同一品种的重复基因"""
    duplicate_gene_dic = {i: [] for i in variety_set}
    for i in variety_set:
        for index, v in enumerate(variety_list):
            if i == v:
                duplicate_gene_dic[i].append(gene_list[index])
    return duplicate_gene_dic

def get_variety_n(query, var_li):
    """严格匹配品种名（前缀匹配）"""
    for v in var_li:
        if query.startswith(v):
            return v
    return None

def assign_gene_by_chrom(variety_li, cluster_li, bed_dic):
    """按染色体对齐不同品种的基因并输出簇列表.
    
    参数:
        variety_li: List[str], 品种名称列表
        cluster_li: List[str], 一个基因簇中的基因ID列表（不是嵌套列表）
        bed_dic: Dict[str, Tuple], 基因位置信息字典，key为基因ID

    返回值:
        List[List[str]]: 每个子列表是一组同簇基因。
    """
    results = []
    my_variety_li = [get_variety_n(g, variety_li) for g in cluster_li]
    variety_s = set(my_variety_li)

    # 如果同一个品种在当前 cluster 中出现多次，则需要按染色体进一步拆分
    if len(variety_s) != len(my_variety_li):
        duplicate_gene_dic = find_duplicate_gene(variety_s, my_variety_li, cluster_li)
        if len(duplicate_gene_dic.keys()) == 1:
            for genes in duplicate_gene_dic.values():
                for g in genes:
                    results.append([g])
        else:
            chrom_s = set()
            chrom_gene_dic_li = []
            for genes in duplicate_gene_dic.values():
                sort_li = sort_gene(genes, bed_dic)
                chrom_gene_dic_li.append(split_gene_by_chrom(sort_li, bed_dic, chrom_s))
            results.extend(assign_gene(chrom_gene_dic_li, chrom_s))
    else:
        # 每个品种在该簇中只出现一次，直接整体作为一个簇返回
        results.append(cluster_li)

    return results

