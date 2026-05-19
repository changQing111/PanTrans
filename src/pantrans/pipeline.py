import os
import logging
from typing import Iterable
from .common_func import (
    get_fasta_len,
    get_bed,
    get_bed_rows,
    cluster2dic,
    extract_fasta_subset_by_names,
    concat_fasta_files,
    concat_text_files,
)
from .align_filter import minimap2_map, filter_bam
from .construct_di_graph import assign_sccs, di_graph_from_pair, get_conn_comp
from .transcript_processor import transcript_dedup, get_cdna_from_gtf, get_subset_bed

logger = logging.getLogger(__name__)


def _normalize_variety_names(variety_input):
    """Normalize variety input into a non-empty list of names."""
    if isinstance(variety_input, str):
        if os.path.isfile(variety_input):
            with open(variety_input, "rt", encoding="utf-8") as f:
                variety_li = [line.strip() for line in f if line.strip()]
        else:
            raw_items = variety_input.replace(";", ",")
            variety_li = [item.strip() for item in raw_items.split(",") if item.strip()]
            if len(variety_li) == 1 and " " in variety_li[0]:
                variety_li = [item for item in variety_li[0].split() if item]
    elif isinstance(variety_input, Iterable):
        variety_li = [str(item).strip() for item in variety_input if str(item).strip()]
    else:
        raise TypeError("variety_li must be a string, path, or iterable of names.")

    if not variety_li:
        raise ValueError("No valid variety names were provided.")
    return variety_li


def _collect_input_genes(bed_dic, variety_li):
    """Collect genes that belong to the requested varieties."""
    variety_prefixes = tuple(variety_li)
    return {
        gene_id
        for gene_id in bed_dic
        if gene_id.startswith(variety_prefixes)
    }

def _load_main_chroms(main_chrom_path):
    if not main_chrom_path:
        return None
    with open(main_chrom_path, "rt", encoding="utf-8") as handle:
        main_chroms = [line.strip() for line in handle if line.strip()]
    if not main_chroms:
        raise ValueError(f"No valid main chromosome names were found in {main_chrom_path}.")
    return main_chroms

def _infer_variety_names_from_bed(bed_path):
    variety_li = []
    with open(bed_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            fields = line.strip().split("\t")
            if len(fields) < 4:
                continue
            gene_id = fields[3]
            variety = gene_id.split(".", 1)[0]
            if variety and variety not in variety_li:
                variety_li.append(variety)
    if not variety_li:
        raise ValueError(f"No valid variety names could be inferred from {bed_path}.")
    return variety_li

def _normalize_chrom_label(chrom):
    chrom = chrom.strip()
    if len(chrom) == 2 and chrom[0] in "ABD" and chrom[1].isdigit():
        return f"{chrom[1]}{chrom[0]}"
    if len(chrom) == 2 and chrom[1] in "ABD" and chrom[0].isdigit():
        return chrom
    return "Ctg"

def _build_gene_rename_map(cluster_dic, bed_path, variety_li, prefix, refer_prefixes=None):
    bed_rows = get_bed_rows(bed_path)
    bed_chrom_by_gene = {}
    chrom_variety_rows = {}
    normalized_variety_li = ["Refer"] + [v for v in variety_li if v != "Refer"] if refer_prefixes else variety_li[:]
    for chrom, start, end, gene_id, strand in bed_rows:
        bed_chrom_by_gene[gene_id] = chrom
        chrom_label = _normalize_chrom_label(chrom)
        chrom_variety_rows.setdefault(chrom_label, {variety: [] for variety in normalized_variety_li})
        if refer_prefixes and gene_id.startswith(tuple(refer_prefixes)):
            chrom_variety_rows[chrom_label]["Refer"].append(gene_id)
            continue
        for variety in normalized_variety_li:
            if variety == "Refer":
                continue
            if gene_id.startswith(variety):
                chrom_variety_rows[chrom_label][variety].append(gene_id)
                break

    row_gene_index = {}
    for chrom_label, variety_gene_rows in chrom_variety_rows.items():
        max_len = max((len(rows) for rows in variety_gene_rows.values()), default=0)
        for row_idx in range(max_len):
            for variety in normalized_variety_li:
                rows = variety_gene_rows[variety]
                if row_idx < len(rows):
                    row_gene_index[rows[row_idx]] = row_idx + 1

    rename_map = {}
    for refer_gene in cluster_dic.keys():
        chrom = bed_chrom_by_gene.get(refer_gene, "Ctg")
        chrom_label = _normalize_chrom_label(chrom)
        if refer_gene in row_gene_index:
            row_num = row_gene_index[refer_gene]
        else:
            used_row_nums = []
            for renamed_gene in rename_map.values():
                if renamed_gene.startswith(f"{prefix}{chrom_label}"):
                    used_row_nums.append(int(renamed_gene[len(prefix + chrom_label):]))
            row_num = max(used_row_nums, default=0) + 1
        rename_map[refer_gene] = f"{prefix}{chrom_label}{row_num:06d}"
    return rename_map

def _write_cluster_outputs(cluster_dic, all_gdna_path, all_bed_path, filter_bam_path,
                           trans_len_dic, gene_len_dic, gene_strand_dic, rename_map, out_dir, prefix, label):
    label_suffix = f"_{label}" if label else ""
    gtf_path = os.path.join(out_dir, f"{prefix}{label_suffix}.gtf")
    gdna_path = os.path.join(out_dir, f"{prefix}{label_suffix}_gdna.refer.fasta")
    cdna_path = os.path.join(out_dir, f"{prefix}{label_suffix}_cdna.refer.fasta")
    bed_path = os.path.join(out_dir, f"{prefix}{label_suffix}.refer.bed")
    refer_gene_li = list(cluster_dic.keys())

    transcript_dedup(
        filter_bam_path,
        cluster_dic=cluster_dic,
        trans_len_dic=trans_len_dic,
        gene_len_dic=gene_len_dic,
        gene_strand_dic=gene_strand_dic,
        rename_map=rename_map,
        gtf_path=gtf_path,
    )
    extract_fasta_subset_by_names(all_gdna_path, refer_gene_li, gdna_path)
    get_cdna_from_gtf(all_gdna_path, gtf_path, cdna_path)
    get_subset_bed(refer_gene_li, all_bed_path, bed_path)
    return gtf_path, cdna_path, gdna_path, bed_path


def unit_construct(all_cdna_path, all_gdna_path, all_bed_path, bam_path, main_chrom_path, variety_li, threads, out_dir, prefix, refer_name):
    variety_li = _normalize_variety_names(variety_li)
    main_chroms = _load_main_chroms(main_chrom_path)
    if not refer_name:
        refer_name = variety_li[0]
    elif refer_name not in variety_li:
        logger.warning("Reference name '%s' is not present in variety list: %s", refer_name, variety_li)

    os.makedirs(out_dir, exist_ok=True)
    bed_dic, gene_strand_dic = get_bed(all_bed_path)
    trans_len_dic = get_fasta_len(all_cdna_path)
    align_bam_path = bam_path or os.path.join(out_dir, f"{prefix}_cdna_align_gdna.bam")
    filter_bam_path = os.path.join(out_dir, f"{prefix}_cdna_align_gdna.filtered.bam")

    if bam_path:
        logger.info("Skip sequence alignment; using existing BAM: %s", bam_path)
    else:
        logger.info("Start sequence alignment: all cDNA vs all gDNA")
        minimap2_map(all_cdna_path, all_gdna_path, threads, align_bam_path)
        logger.info("Finish sequence alignment")
    aligned_gene_li, gene_len_dic = filter_bam(align_bam_path, filter_bam_path, bed_dic)
    logger.info("Finish BAM filtering, generated %s", filter_bam_path)
    # Step 2: Build graph and derive final gene clusters
    logger.info("Start graph building and gene assignment")
    G = di_graph_from_pair(aligned_gene_li)
    sccs = get_conn_comp(G)

    pre_clusters, last_clusters = assign_sccs(
        sccs, G, gene_len_dic, bed_dic, variety_li, refer_name, main_chroms=main_chroms
    )
    input_gene_set = _collect_input_genes(bed_dic, variety_li)
    clustered_gene_set = {gene_id for cluster in last_clusters for gene_id in cluster}
    missing_gene_clusters = [[gene_id] for gene_id in sorted(input_gene_set - clustered_gene_set)]
    if missing_gene_clusters:
        logger.warning(
            "Recovered %d genes that were absent from the filtered alignment graph.",
            len(missing_gene_clusters),
        )
        pre_clusters.extend(missing_gene_clusters)
        last_clusters.extend(missing_gene_clusters)
    last_cluster_dic = cluster2dic(last_clusters)
    pre_cluster_dic = cluster2dic(pre_clusters)
    last_rename_map = _build_gene_rename_map(last_cluster_dic, all_bed_path, variety_li, "Pan")
    # write pre last clusters into file
    pre_cluster_path = os.path.join(out_dir, f"{prefix}_pre.tmp.cluster")
    with open(pre_cluster_path, "w") as f:
        for li in pre_cluster_dic.values():
            f.write("\t".join(li) + "\n")

    last_cluster_path = os.path.join(out_dir, f"{prefix}_last.tmp.cluster")
    with open(last_cluster_path, "w") as f:
        for li in last_cluster_dic.values():
            f.write("\t".join(li) + "\n")

    logger.info("Produce pre-cluster transcript and sequence outputs")
    _write_cluster_outputs(
        pre_cluster_dic,
        all_gdna_path,
        all_bed_path,
        filter_bam_path,
        trans_len_dic,
        gene_len_dic,
        gene_strand_dic,
        None,
        out_dir,
        prefix,
        "pre",
    )

    logger.info("Produce last-cluster transcript and sequence outputs")
    gtf_path, new_cdna_path, new_gdna_path, new_bed_path = _write_cluster_outputs(
        last_cluster_dic,
        all_gdna_path,
        all_bed_path,
        filter_bam_path,
        trans_len_dic,
        gene_len_dic,
        gene_strand_dic,
        last_rename_map,
        out_dir,
        prefix,
        "",
    )
    logger.info(f"{variety_li[-1]} construct end")
    return new_cdna_path, new_gdna_path, new_bed_path

def unit_append(query_cdna_path, query_gdna_path, query_bed_path, refer_cdna_path, refer_gdna_path, refer_bed_path,
                variety_name, threads, out_dir, prefix="Append"):
    os.makedirs(out_dir, exist_ok=True)
    new_variety_li = _normalize_variety_names(variety_name)

    merged_cdna_path = os.path.join(out_dir, f"{prefix}_merged.cdna.fasta")
    merged_gdna_path = os.path.join(out_dir, f"{prefix}_merged.gdna.fasta")
    merged_bed_path = os.path.join(out_dir, f"{prefix}_merged.bed")
    merged_bam_path = os.path.join(out_dir, f"{prefix}_merged_cdna_align_gdna.bam")
    filtered_bam_path = os.path.join(out_dir, f"{prefix}_merged_cdna_align_gdna.filtered.bam")

    logger.info("Merge query and reference cDNA into %s", merged_cdna_path)
    concat_fasta_files([refer_cdna_path, query_cdna_path], merged_cdna_path)
    logger.info("Merge query and reference gDNA into %s", merged_gdna_path)
    concat_fasta_files([refer_gdna_path, query_gdna_path], merged_gdna_path)
    logger.info("Merge query and reference BED into %s", merged_bed_path)
    concat_text_files([refer_bed_path, query_bed_path], merged_bed_path)

    logger.info("Start append alignment: merged cDNA vs merged gDNA")
    minimap2_map(merged_cdna_path, merged_gdna_path, threads, merged_bam_path)
    logger.info("Finish append alignment")

    bed_dic, gene_strand_dic = get_bed(merged_bed_path)
    trans_len_dic = get_fasta_len(merged_cdna_path)
    logger.info("Start append BAM filtering")
    aligned_gene_li, gene_len_dic = filter_bam(merged_bam_path, filtered_bam_path, bed_dic)
    logger.info("Finish append BAM filtering, generated %s", filtered_bam_path)

    refer_variety_li = _infer_variety_names_from_bed(refer_bed_path)
    variety_li = new_variety_li[:]
    refer_name = refer_variety_li[0]

    logger.info("Start append graph building and gene assignment")
    G = di_graph_from_pair(aligned_gene_li)
    sccs = get_conn_comp(G)
    pre_clusters, last_clusters = assign_sccs(
        sccs, G, gene_len_dic, bed_dic, variety_li, refer_name, main_chroms=None, refer_prefixes=refer_variety_li
    )

    input_gene_set = _collect_input_genes(bed_dic, variety_li)
    clustered_gene_set = {gene_id for cluster in last_clusters for gene_id in cluster}
    missing_gene_clusters = [[gene_id] for gene_id in sorted(input_gene_set - clustered_gene_set)]
    if missing_gene_clusters:
        logger.warning(
            "Recovered %d genes that were absent from the filtered alignment graph.",
            len(missing_gene_clusters),
        )
        pre_clusters.extend(missing_gene_clusters)
        last_clusters.extend(missing_gene_clusters)

    pre_cluster_dic = cluster2dic(pre_clusters)
    last_cluster_dic = cluster2dic(last_clusters)
    last_rename_map = _build_gene_rename_map(
        last_cluster_dic, merged_bed_path, variety_li, "Pan", refer_prefixes=refer_variety_li
    )

    pre_cluster_path = os.path.join(out_dir, f"{prefix}_pre.tmp.cluster")
    with open(pre_cluster_path, "w") as handle:
        for li in pre_cluster_dic.values():
            handle.write("\t".join(li) + "\n")

    last_cluster_path = os.path.join(out_dir, f"{prefix}_last.tmp.cluster")
    with open(last_cluster_path, "w") as handle:
        for li in last_cluster_dic.values():
            handle.write("\t".join(li) + "\n")

    logger.info("Produce append pre-cluster transcript and sequence outputs")
    _write_cluster_outputs(
        pre_cluster_dic,
        merged_gdna_path,
        merged_bed_path,
        filtered_bam_path,
        trans_len_dic,
        gene_len_dic,
        gene_strand_dic,
        None,
        out_dir,
        prefix,
        "pre",
    )

    logger.info("Produce append last-cluster transcript and sequence outputs")
    _write_cluster_outputs(
        last_cluster_dic,
        merged_gdna_path,
        merged_bed_path,
        filtered_bam_path,
        trans_len_dic,
        gene_len_dic,
        gene_strand_dic,
        last_rename_map,
        out_dir,
        prefix,
        "",
    )

    return merged_cdna_path, merged_gdna_path, merged_bed_path, merged_bam_path, filtered_bam_path
