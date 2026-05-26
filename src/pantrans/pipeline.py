import os
import shutil
import subprocess
import logging
from typing import Iterable
from .common_func import (
    get_fasta_len,
    get_bed,
    get_bed_rows,
    cluster2dic,
    extract_fasta_subset_by_names,
    extract_transcripts_by_gene_names,
    extract_fasta_records_by_exact_names,
    build_transcript_index,
    write_transcripts_from_index,
    extract_fasta_records_by_exact_names_indexed,
    concat_fasta_files,
    concat_text_files,
)
from .align_filter import minimap2_map, minimap2_map_rescue, filter_bam
from .construct_di_graph import assign_sccs, di_graph_from_pair, get_conn_comp
from .transcript_processor import transcript_dedup, get_cdna_from_gtf, get_subset_bed, sort_gtf_by_gene_id

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
    bed_order_by_gene = {}
    chrom_variety_rows = {}
    normalized_variety_li = ["Refer"] + [v for v in variety_li if v != "Refer"] if refer_prefixes else variety_li[:]
    for index, (chrom, start, end, gene_id, strand) in enumerate(bed_rows):
        bed_chrom_by_gene[gene_id] = chrom
        bed_order_by_gene[gene_id] = index
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

    chrom_hub_dic = {}
    for refer_gene in cluster_dic.keys():
        chrom = bed_chrom_by_gene.get(refer_gene, "Ctg")
        chrom_label = _normalize_chrom_label(chrom)
        chrom_hub_dic.setdefault(chrom_label, []).append(refer_gene)

    rename_map = {}
    for chrom_label in sorted(chrom_hub_dic.keys()):
        sorted_hubs = sorted(
            chrom_hub_dic[chrom_label],
            key=lambda gene_id: (
                row_gene_index.get(gene_id, float("inf")),
                bed_order_by_gene.get(gene_id, float("inf")),
                gene_id,
            ),
        )
        for index, refer_gene in enumerate(sorted_hubs, start=1):
            rename_map[refer_gene] = f"{prefix}{chrom_label}{index:06d}"
    return rename_map

def _write_cluster_outputs(cluster_dic, all_gdna_path, all_bed_path, filter_bam_path,
                           trans_len_dic, gene_len_dic, gene_strand_dic, rename_map, out_dir, prefix, label,
                           all_cdna_path=None, threads=8, enable_rescue=False, pre_gtf_path=None):
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
    if enable_rescue:
        _rescue_missing_cluster_genes(
            cluster_dic=cluster_dic,
            gtf_path=gtf_path,
            all_cdna_path=all_cdna_path,
            all_gdna_path=all_gdna_path,
            all_bed_path=all_bed_path,
            rename_map=rename_map,
            out_dir=out_dir,
            prefix=prefix,
            label=label,
            threads=threads,
            pre_gtf_path=pre_gtf_path,
        )
    if rename_map:
        sort_gtf_by_gene_id(gtf_path)
    extract_fasta_subset_by_names(all_gdna_path, refer_gene_li, gdna_path)
    get_cdna_from_gtf(all_gdna_path, gtf_path, cdna_path)
    get_subset_bed(refer_gene_li, all_bed_path, bed_path)
    return gtf_path, cdna_path, gdna_path, bed_path

def _get_gtf_hub_genes(gtf_path):
    genes = set()
    with open(gtf_path, "rt") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) == 9 and fields[2] == "transcript":
                genes.add(fields[0])
    return genes

def _append_file_contents(src_path, dst_path):
    with open(dst_path, "a") as dst, open(src_path, "rt") as src:
        for line in src:
            dst.write(line)

def _merge_bam_files(bam_paths, merged_bam_path):
    if not bam_paths:
        return False
    if len(bam_paths) == 1:
        shutil.copyfile(bam_paths[0], merged_bam_path)
        return True
    subprocess.run(["samtools", "merge", "-f", merged_bam_path] + bam_paths, check=True)
    return True

def _parse_gtf_attrs(attrs):
    parsed = {}
    for item in attrs.strip().split(";"):
        item = item.strip()
        if not item or " " not in item:
            continue
        key, value = item.replace('"', "").split(" ", 1)
        parsed[key] = value
    return parsed

def _format_gtf_attrs(transcript_id, gene_id):
    return f'transcript_id "{transcript_id}"; gene_id "{gene_id}";'

def _extract_gtf_records_by_gene(gtf_path, gene_ids):
    records_by_gene = {gene_id: [] for gene_id in gene_ids}
    with open(gtf_path, "rt") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] not in {"transcript", "exon"}:
                continue
            gene_id = fields[0]
            if gene_id in records_by_gene:
                records_by_gene[gene_id].append(fields)
    return records_by_gene

def _append_rescue_records_from_gtf(dst_gtf_path, rescue_records_by_gene, rename_map):
    rescued_genes = set()
    with open(dst_gtf_path, "a") as dst:
        for gene_id in sorted(rescue_records_by_gene):
            records = rescue_records_by_gene[gene_id]
            if not records:
                continue
            output_gene_id = rename_map.get(gene_id, gene_id) if rename_map else gene_id
            transcript_counter = 0
            transcript_id_map = {}
            for fields in records:
                new_fields = fields[:]
                attrs = _parse_gtf_attrs(new_fields[8])
                old_transcript_id = attrs.get("transcript_id", "")
                if new_fields[2] == "transcript":
                    transcript_counter += 1
                    transcript_id_map[old_transcript_id] = f"{output_gene_id}.{transcript_counter}"
                new_transcript_id = transcript_id_map.get(old_transcript_id)
                if not new_transcript_id:
                    transcript_counter = max(transcript_counter, 1)
                    new_transcript_id = f"{output_gene_id}.{transcript_counter}"
                    transcript_id_map[old_transcript_id] = new_transcript_id
                new_fields[8] = _format_gtf_attrs(new_transcript_id, output_gene_id)
                dst.write("\t".join(new_fields) + "\n")
            rescued_genes.add(gene_id)
    return rescued_genes

def _rescue_missing_cluster_genes(cluster_dic, gtf_path, all_cdna_path, all_gdna_path, all_bed_path,
                                  rename_map, out_dir, prefix, label, threads, pre_gtf_path=None):
    gtf_gene_set = _get_gtf_hub_genes(gtf_path)
    missing_cluster_dic = {
        hub: members
        for hub, members in cluster_dic.items()
        if hub not in gtf_gene_set
    }
    if not missing_cluster_dic:
        return 0

    rescued_from_pre = set()
    if pre_gtf_path and os.path.exists(pre_gtf_path):
        pre_records_by_gene = _extract_gtf_records_by_gene(pre_gtf_path, missing_cluster_dic.keys())
        rescued_from_pre = _append_rescue_records_from_gtf(gtf_path, pre_records_by_gene, rename_map)
        if rescued_from_pre:
            logger.info(
                "Recovered %d missing last-cluster genes directly from pre-cluster GTF.",
                len(rescued_from_pre),
            )
    missing_cluster_dic = {
        hub: members
        for hub, members in missing_cluster_dic.items()
        if hub not in rescued_from_pre
    }
    if not missing_cluster_dic:
        return len(rescued_from_pre)

    rescue_gene_names = sorted({gene for members in missing_cluster_dic.values() for gene in members})

    label_suffix = f"_{label}" if label else ""
    rescue_cdna_path = os.path.join(out_dir, f"{prefix}{label_suffix}_rescue.cdna.fasta")
    rescue_gdna_path = os.path.join(out_dir, f"{prefix}{label_suffix}_rescue.gdna.fasta")
    rescue_bed_path = os.path.join(out_dir, f"{prefix}{label_suffix}_rescue.bed")
    rescue_bam_path = os.path.join(out_dir, f"{prefix}{label_suffix}_rescue_cdna_align_gdna.bam")
    rescue_filtered_bam_path = os.path.join(out_dir, f"{prefix}{label_suffix}_rescue_cdna_align_gdna.filtered.bam")
    rescue_gtf_path = os.path.join(out_dir, f"{prefix}{label_suffix}_rescue.gtf")
    rescue_work_dir = os.path.join(out_dir, f"{prefix}{label_suffix}_rescue_work")
    transcript_index = build_transcript_index(all_cdna_path)

    kept_transcripts = write_transcripts_from_index(transcript_index, rescue_gene_names, rescue_cdna_path)
    if kept_transcripts == 0:
        logger.warning(
            "No rescue transcripts were found for %d missing clusters (%d genes).",
            len(missing_cluster_dic),
            len(rescue_gene_names),
        )
        return 0
    kept_genes = extract_fasta_records_by_exact_names_indexed(all_gdna_path, rescue_gene_names, rescue_gdna_path)
    if kept_genes == 0:
        logger.warning(
            "No rescue genomic sequences were found for %d missing clusters (%d genes).",
            len(missing_cluster_dic),
            len(rescue_gene_names),
        )
        return 0
    get_subset_bed(rescue_gene_names, all_bed_path, rescue_bed_path)

    logger.info(
        "Run rescue alignment for %d missing clusters covering %d genes",
        len(missing_cluster_dic),
        len(rescue_gene_names),
    )
    if os.path.isdir(rescue_work_dir):
        shutil.rmtree(rescue_work_dir)
    os.makedirs(rescue_work_dir)

    single_cluster_dic = {hub: members for hub, members in missing_cluster_dic.items() if len(members) == 1}
    multi_cluster_dic = {hub: members for hub, members in missing_cluster_dic.items() if len(members) > 1}
    batch_raw_bams = []
    batch_filtered_bams = []

    for index, (hub, members) in enumerate(sorted(single_cluster_dic.items()), start=1):
        batch_prefix = f"single_{index:06d}"
        batch_cdna_path = os.path.join(rescue_work_dir, f"{batch_prefix}.cdna.fasta")
        batch_gdna_path = os.path.join(rescue_work_dir, f"{batch_prefix}.gdna.fasta")
        batch_bed_path = os.path.join(rescue_work_dir, f"{batch_prefix}.bed")
        batch_bam_path = os.path.join(rescue_work_dir, f"{batch_prefix}.bam")
        batch_filtered_bam_path = os.path.join(rescue_work_dir, f"{batch_prefix}.filtered.bam")

        write_transcripts_from_index(transcript_index, members, batch_cdna_path)
        extract_fasta_records_by_exact_names_indexed(all_gdna_path, members, batch_gdna_path)
        get_subset_bed(members, all_bed_path, batch_bed_path)
        minimap2_map_rescue(batch_cdna_path, batch_gdna_path, threads, batch_bam_path)
        batch_bed_dic, _ = get_bed(batch_bed_path)
        filter_bam(batch_bam_path, batch_filtered_bam_path, batch_bed_dic)
        batch_raw_bams.append(batch_bam_path)
        batch_filtered_bams.append(batch_filtered_bam_path)

    if multi_cluster_dic:
        multi_gene_names = sorted({gene for members in multi_cluster_dic.values() for gene in members})
        batch_cdna_path = os.path.join(rescue_work_dir, "multi.cdna.fasta")
        batch_gdna_path = os.path.join(rescue_work_dir, "multi.gdna.fasta")
        batch_bed_path = os.path.join(rescue_work_dir, "multi.bed")
        batch_bam_path = os.path.join(rescue_work_dir, "multi.bam")
        batch_filtered_bam_path = os.path.join(rescue_work_dir, "multi.filtered.bam")

        write_transcripts_from_index(transcript_index, multi_gene_names, batch_cdna_path)
        extract_fasta_records_by_exact_names_indexed(all_gdna_path, multi_gene_names, batch_gdna_path)
        get_subset_bed(multi_gene_names, all_bed_path, batch_bed_path)
        minimap2_map_rescue(batch_cdna_path, batch_gdna_path, threads, batch_bam_path)
        batch_bed_dic, _ = get_bed(batch_bed_path)
        filter_bam(batch_bam_path, batch_filtered_bam_path, batch_bed_dic)
        batch_raw_bams.append(batch_bam_path)
        batch_filtered_bams.append(batch_filtered_bam_path)

    if not _merge_bam_files(batch_raw_bams, rescue_bam_path):
        shutil.rmtree(rescue_work_dir)
        return len(rescued_from_pre)
    if not _merge_bam_files(batch_filtered_bams, rescue_filtered_bam_path):
        shutil.rmtree(rescue_work_dir)
        return len(rescued_from_pre)

    for bam_path in batch_raw_bams + batch_filtered_bams:
        if os.path.exists(bam_path):
            os.remove(bam_path)

    rescue_gene_strand_dic = get_bed(rescue_bed_path)[1]
    rescue_gene_len_dic = get_fasta_len(rescue_gdna_path)
    rescue_trans_len_dic = get_fasta_len(rescue_cdna_path)
    rescue_cluster_dic = missing_cluster_dic
    rescue_rename_map = {
        gene: rename_map[gene]
        for gene in missing_cluster_dic
        if rename_map and gene in rename_map
    }

    transcript_dedup(
        rescue_filtered_bam_path,
        cluster_dic=rescue_cluster_dic,
        trans_len_dic=rescue_trans_len_dic,
        gene_len_dic=rescue_gene_len_dic,
        gene_strand_dic=rescue_gene_strand_dic,
        rename_map=rescue_rename_map if rename_map else None,
        gtf_path=rescue_gtf_path,
    )

    rescued_gene_set = _get_gtf_hub_genes(rescue_gtf_path) if os.path.exists(rescue_gtf_path) else set()
    if rescued_gene_set:
        _append_file_contents(rescue_gtf_path, gtf_path)
        logger.info("Rescued %d previously missing genes into %s", len(rescued_gene_set), gtf_path)
    else:
        logger.warning("Rescue alignment produced no additional GTF entries.")
    shutil.rmtree(rescue_work_dir)
    return len(rescued_gene_set) + len(rescued_from_pre)


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
    pre_gtf_path, _, _, _ = _write_cluster_outputs(
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
        all_cdna_path=all_cdna_path,
        threads=threads,
        enable_rescue=True,
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
        all_cdna_path=all_cdna_path,
        threads=threads,
        enable_rescue=True,
        pre_gtf_path=pre_gtf_path,
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
    pre_gtf_path, _, _, _ = _write_cluster_outputs(
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
        all_cdna_path=merged_cdna_path,
        threads=threads,
        enable_rescue=True,
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
        all_cdna_path=merged_cdna_path,
        threads=threads,
        enable_rescue=True,
        pre_gtf_path=pre_gtf_path,
    )

    return merged_cdna_path, merged_gdna_path, merged_bed_path, merged_bam_path, filtered_bam_path
