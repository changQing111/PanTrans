import os
import shutil
import subprocess
import logging
from itertools import chain
from typing import Iterable
from .common_func import (
    get_fasta_len,
    get_fasta_sequence_identities,
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
)
from .align_filter import (
    minimap2_map,
    minimap2_map_rescue,
    filter_bam,
    bam_alignment_provenance,
    merge_bams_with_full_header,
    alignment_provenance,
    transcript_to_gene_id,
    iter_filtered_bam_gene_edges,
    FILTER_LOGIC_ID,
    validate_resume_bam,
)
from .construct_di_graph import assign_sccs, di_graph_from_pair, get_conn_comp
from .graph_package import (
    iter_graph_edges,
    load_graph_package,
    merge_history_and_query_bed,
    validate_edge_provenance,
    write_graph_package,
)
from .transcript_processor import (
    get_cdna_from_gtf,
    get_subset_bed,
    load_gtf_transcript_models,
    rename_gtf_ids,
    sort_gtf_by_gene_id,
    transcript_dedup,
)

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


def _validate_variety_gene_ids(gene_ids, variety_li, label):
    """Require every gene ID to match one unambiguous variety prefix."""
    variety_li = _normalize_variety_names(variety_li)
    if len(set(variety_li)) != len(variety_li):
        raise ValueError(f"{label} contains duplicate variety names")
    ambiguous_pairs = sorted(
        {
            tuple(sorted((left, right)))
            for index, left in enumerate(variety_li)
            for right in variety_li[index + 1 :]
            if left.startswith(right) or right.startswith(left)
        }
    )
    if ambiguous_pairs:
        formatted = ", ".join(f"{left}/{right}" for left, right in ambiguous_pairs)
        raise ValueError(f"{label} has ambiguous variety prefixes: {formatted}")

    matched_varieties = set()
    invalid_genes = []
    for gene_id in gene_ids:
        matches = [variety for variety in variety_li if gene_id.startswith(variety)]
        if len(matches) != 1:
            invalid_genes.append(gene_id)
        else:
            matched_varieties.add(matches[0])
    if invalid_genes:
        preview = ", ".join(sorted(invalid_genes)[:10])
        suffix = "..." if len(invalid_genes) > 10 else ""
        raise ValueError(
            f"{label} gene ID does not match exactly one variety prefix: "
            f"{preview}{suffix}"
        )
    missing_varieties = sorted(set(variety_li) - matched_varieties)
    if missing_varieties:
        raise ValueError(
            f"{label} has no genes for varieties: {', '.join(missing_varieties)}"
        )


def _derive_graph_and_clusters(
    aligned_gene_li,
    gene_len_dic,
    bed_dic,
    variety_li,
    refer_name,
    eligible_gene_set,
    main_chroms=None,
    refer_prefixes=None,
):
    graph = di_graph_from_pair(aligned_gene_li)
    sccs = get_conn_comp(graph)
    eligible_gene_set = set(eligible_gene_set)
    logger.info(
        "Cluster graph contains %d nodes; eligible gene set contains %d genes.",
        graph.number_of_nodes(),
        len(eligible_gene_set),
    )
    pre_clusters, last_clusters = assign_sccs(
        sccs,
        graph,
        gene_len_dic,
        bed_dic,
        variety_li,
        refer_name,
        main_chroms=main_chroms,
        refer_prefixes=refer_prefixes,
    )
    pre_clustered = {gene for cluster in pre_clusters for gene in cluster}
    last_clustered = {gene for cluster in last_clusters for gene in cluster}
    missing_pre = [[gene] for gene in sorted(eligible_gene_set - pre_clustered)]
    missing_last = [[gene] for gene in sorted(eligible_gene_set - last_clustered)]
    if missing_pre or missing_last:
        logger.warning(
            "Recovered %d pre-cluster genes and %d last-cluster genes absent from the filtered graph.",
            len(missing_pre),
            len(missing_last),
        )
        pre_clusters.extend(missing_pre)
        last_clusters.extend(missing_last)
    return graph, pre_clusters, last_clusters


def _derive_clusters_from_alignment(
    aligned_gene_li,
    gene_len_dic,
    bed_dic,
    variety_li,
    refer_name,
    eligible_gene_set,
    main_chroms=None,
    refer_prefixes=None,
):
    _, pre_clusters, last_clusters = _derive_graph_and_clusters(
        aligned_gene_li=aligned_gene_li,
        gene_len_dic=gene_len_dic,
        bed_dic=bed_dic,
        variety_li=variety_li,
        refer_name=refer_name,
        eligible_gene_set=eligible_gene_set,
        main_chroms=main_chroms,
        refer_prefixes=refer_prefixes,
    )
    return pre_clusters, last_clusters


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


def _infer_reference_variety_names_from_bed(bed_path, query_variety_li):
    """Infer historical reference prefixes from a representative-plus-query BED.

    Append receives a single BED containing historical representative genes and
    the genes of the newly appended varieties.  Query variety names can be
    prefixes (for example ``JM22`` for ``JM22A1.g1``), so they are removed by
    prefix match before collecting the first dot-delimited component of each
    remaining gene identifier.
    """
    query_variety_li = _normalize_variety_names(query_variety_li)
    query_prefixes = tuple(query_variety_li)
    reference_variety_li = []
    with open(bed_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4 or not fields[3]:
                continue
            gene_id = fields[3]
            if gene_id.startswith(query_prefixes):
                continue
            variety = gene_id.split(".", 1)[0]
            if variety and variety not in reference_variety_li:
                reference_variety_li.append(variety)
    if not reference_variety_li:
        raise ValueError(
            "No reference variety prefixes could be inferred from append BED: "
            f"{bed_path}."
        )
    return reference_variety_li


def _validate_append_gene_set(sequence_gene_set, bed_dic):
    """Return sequence-backed append genes after enforcing BED consistency."""
    sequence_gene_set = set(sequence_gene_set)
    bed_gene_set = set(bed_dic)
    missing_bed_genes = sorted(sequence_gene_set - bed_gene_set)
    if missing_bed_genes:
        preview = ", ".join(missing_bed_genes[:10])
        suffix = "..." if len(missing_bed_genes) > 10 else ""
        raise ValueError(
            "Merged gDNA genes are missing from append BED: "
            f"{preview}{suffix}"
        )

    bed_only_genes = bed_gene_set - sequence_gene_set
    if bed_only_genes:
        logger.warning(
            "Excluded %d append BED-only genes that have no merged gDNA sequence.",
            len(bed_only_genes),
        )
    return sequence_gene_set


def _validate_history_transcriptome(
    history_gtf_path, history_cdna_path, history_gene_ids
):
    """Validate an unrenamed final GTF/cDNA pair against the graph package."""
    models = load_gtf_transcript_models(history_gtf_path)
    gtf_transcript_ids = set(models["transcript_gene"])
    cdna_transcript_ids = set(get_fasta_len(history_cdna_path))
    if gtf_transcript_ids != cdna_transcript_ids:
        missing_from_cdna = sorted(gtf_transcript_ids - cdna_transcript_ids)
        missing_from_gtf = sorted(cdna_transcript_ids - gtf_transcript_ids)
        raise ValueError(
            "Historical GTF and cDNA transcript IDs do not match; "
            f"missing from cDNA: {missing_from_cdna[:10]}; "
            f"missing from GTF: {missing_from_gtf[:10]}"
        )

    history_gtf_genes = set(models["transcript_gene"].values())
    renamed_genes = sorted(
        gene_id for gene_id in history_gtf_genes if gene_id.startswith("Pan")
    )
    if renamed_genes:
        raise ValueError(
            "Historical GTF contains renamed Pan gene IDs; use the previous "
            f"unrenamed GTF: {', '.join(renamed_genes[:10])}"
        )

    genes_outside_graph = sorted(history_gtf_genes - set(history_gene_ids))
    if genes_outside_graph:
        raise ValueError(
            "Historical GTF contains genes absent from historical graph: "
            + ", ".join(genes_outside_graph[:10])
        )
    return models


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
                           all_cdna_path=None, threads=8, enable_rescue=False, pre_gtf_path=None,
                           seed_gtf_path=None):
    label_suffix = f"_{label}" if label else ""
    gtf_path = os.path.join(out_dir, f"{prefix}{label_suffix}.gtf")
    gdna_path = os.path.join(out_dir, f"{prefix}{label_suffix}_gdna.refer.fasta")
    cdna_path = os.path.join(out_dir, f"{prefix}{label_suffix}_cdna.refer.fasta")
    bed_path = os.path.join(out_dir, f"{prefix}{label_suffix}.refer.bed")
    refer_gene_li = list(cluster_dic.keys())
    unrenamed_gtf_path = os.path.join(
        out_dir, f"{prefix}{label_suffix}_unrenamed.gtf"
    )
    transcript_gtf_path = unrenamed_gtf_path if rename_map else gtf_path

    transcript_dedup(
        filter_bam_path,
        cluster_dic=cluster_dic,
        trans_len_dic=trans_len_dic,
        gene_len_dic=gene_len_dic,
        gene_strand_dic=gene_strand_dic,
        rename_map=None,
        gtf_path=transcript_gtf_path,
        seed_gtf_path=seed_gtf_path,
    )
    if enable_rescue:
        _rescue_missing_cluster_genes(
            cluster_dic=cluster_dic,
            gtf_path=transcript_gtf_path,
            all_cdna_path=all_cdna_path,
            all_gdna_path=all_gdna_path,
            all_bed_path=all_bed_path,
            rename_map=None,
            out_dir=out_dir,
            prefix=prefix,
            label=label,
            threads=threads,
            pre_gtf_path=pre_gtf_path,
        )
    if rename_map:
        rename_gtf_ids(unrenamed_gtf_path, gtf_path, rename_map)
        sort_gtf_by_gene_id(gtf_path)
    extract_fasta_subset_by_names(all_gdna_path, refer_gene_li, gdna_path)
    if rename_map:
        unrenamed_cdna_path = os.path.join(
            out_dir, f"{prefix}{label_suffix}_unrenamed_cdna.refer.fasta"
        )
        get_cdna_from_gtf(all_gdna_path, unrenamed_gtf_path, unrenamed_cdna_path)
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
    input_gene_set = _collect_input_genes(bed_dic, variety_li)
    graph, pre_clusters, last_clusters = _derive_graph_and_clusters(
        aligned_gene_li=aligned_gene_li,
        gene_len_dic=gene_len_dic,
        bed_dic=bed_dic,
        variety_li=variety_li,
        refer_name=refer_name,
        eligible_gene_set=input_gene_set,
        main_chroms=main_chroms,
    )
    graph_package_path = os.path.join(out_dir, f"{prefix}.graph.json")
    construct_provenance = (
        bam_alignment_provenance(align_bam_path)
        if bam_path
        else alignment_provenance()
    )
    construct_provenance["scope"] = "construct_all_to_all"
    construct_provenance["filter_thresholds_assumed"] = False
    construct_provenance["filter_thresholds_source"] = "current_construct_filter"
    construct_provenance["filter_logic_id"] = FILTER_LOGIC_ID
    write_graph_package(
        manifest_path=graph_package_path,
        edge_iter=graph.edges(),
        gene_len_dic=gene_len_dic,
        bed_path=all_bed_path,
        filtered_bam_path=filter_bam_path,
        cdna_paths=[all_cdna_path],
        gdna_paths=[all_gdna_path],
        variety_names=variety_li,
        reference_name=refer_name,
        main_chroms=main_chroms,
        provenance={"edge_generations": [construct_provenance]},
    )
    logger.info("Wrote reusable graph package: %s", graph_package_path)
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

def unit_append(
    history_cdna_path,
    history_gtf_path,
    query_cdna_path,
    query_gdna_path,
    all_bed_path,
    history_graph_path,
    variety_name,
    threads,
    out_dir,
    prefix="Append",
    query_to_all_bam=None,
    history_to_query_bam=None,
):
    """Append new varieties by reusing the complete historical graph."""
    os.makedirs(out_dir, exist_ok=True)
    new_variety_li = _normalize_variety_names(variety_name)
    history_package = load_graph_package(history_graph_path)
    current_provenance = alignment_provenance()
    history_generations = validate_edge_provenance(
        history_package.get("provenance"), current_provenance
    )
    history_variety_li = list(history_package["variety_names"])
    duplicate_varieties = sorted(set(history_variety_li) & set(new_variety_li))
    if duplicate_varieties:
        raise ValueError(
            "Append variety names already exist in the historical graph: "
            + ", ".join(duplicate_varieties)
        )
    variety_li = history_variety_li + new_variety_li

    merged_cdna_path = os.path.join(out_dir, f"{prefix}_merged.cdna.fasta")
    merged_gdna_path = os.path.join(out_dir, f"{prefix}_merged.gdna.fasta")
    merged_bed_path = os.path.join(out_dir, f"{prefix}_merged.bed")
    query_to_all_bam_path = query_to_all_bam or os.path.join(
        out_dir, f"{prefix}_query_to_all.bam"
    )
    query_to_all_filtered_bam_path = os.path.join(
        out_dir, f"{prefix}_query_to_all.filtered.bam"
    )
    history_to_query_bam_path = history_to_query_bam or os.path.join(
        out_dir, f"{prefix}_history_to_query.bam"
    )
    history_to_query_filtered_bam_path = os.path.join(
        out_dir, f"{prefix}_history_to_query.filtered.bam"
    )
    filtered_bam_path = os.path.join(out_dir, f"{prefix}_merged_cdna_align_gdna.filtered.bam")

    query_gene_len_dic = get_fasta_len(query_gdna_path)
    query_gene_set = set(query_gene_len_dic)
    query_trans_len_dic = get_fasta_len(query_cdna_path)
    query_transcript_ids = set(query_trans_len_dic)
    query_cdna_gene_set = {
        transcript_to_gene_id(transcript_id)
        for transcript_id in query_transcript_ids
    }
    history_gene_set = set(history_package["history_gene_ids"])
    _validate_history_transcriptome(
        history_gtf_path, history_cdna_path, history_gene_set
    )
    _validate_variety_gene_ids(
        history_gene_set, history_variety_li, "historical graph"
    )
    _validate_variety_gene_ids(query_gene_set, new_variety_li, "query gDNA")
    _validate_variety_gene_ids(query_cdna_gene_set, new_variety_li, "query cDNA")
    unexpected_query_cdna = sorted(query_cdna_gene_set - query_gene_set)
    if unexpected_query_cdna:
        preview = ", ".join(unexpected_query_cdna[:10])
        suffix = "..." if len(unexpected_query_cdna) > 10 else ""
        raise ValueError(
            "Query cDNA contains genes absent from query gDNA: "
            f"{preview}{suffix}"
        )
    _validate_variety_gene_ids(
        history_gene_set | query_gene_set, variety_li, "merged append input"
    )
    duplicate_genes = sorted(history_gene_set & query_gene_set)
    if duplicate_genes:
        preview = ", ".join(duplicate_genes[:10])
        suffix = "..." if len(duplicate_genes) > 10 else ""
        raise ValueError(f"Append gene IDs already exist in history: {preview}{suffix}")

    logger.info("Prepare full historical and merged sequence inputs")
    concat_fasta_files([history_cdna_path, query_cdna_path], merged_cdna_path)
    concat_fasta_files(history_package["gdna_paths"] + [query_gdna_path], merged_gdna_path)
    merge_history_and_query_bed(
        history_package["bed_path"],
        all_bed_path,
        history_gene_ids=history_gene_set,
        query_gene_ids=query_gene_set,
        output_path=merged_bed_path,
    )
    bed_dic, gene_strand_dic = get_bed(merged_bed_path)
    trans_len_dic = get_fasta_len(merged_cdna_path)
    eligible_gene_set = history_gene_set | query_gene_set
    _validate_append_gene_set(eligible_gene_set, bed_dic)
    cdna_gene_set = {
        transcript_to_gene_id(transcript_id) for transcript_id in trans_len_dic
    }
    unexpected_cdna_genes = sorted(cdna_gene_set - eligible_gene_set)
    if unexpected_cdna_genes:
        preview = ", ".join(unexpected_cdna_genes[:10])
        suffix = "..." if len(unexpected_cdna_genes) > 10 else ""
        raise ValueError(
            "Merged cDNA contains genes absent from merged gDNA/BED: "
            f"{preview}{suffix}"
        )
    history_transcript_ids = set(get_fasta_len(history_cdna_path))
    query_resume_expectations = (
        get_fasta_sequence_identities(query_cdna_path)
        if query_to_all_bam
        else query_transcript_ids
    )
    history_resume_expectations = (
        get_fasta_sequence_identities(history_cdna_path)
        if history_to_query_bam
        else history_transcript_ids
    )
    expected_merged_gene_lengths = dict(history_package["gene_len_dic"])
    expected_merged_gene_lengths.update(query_gene_len_dic)
    for resume_bam_path, expected_queries, expected_targets, label in (
        (
            query_to_all_bam,
            query_resume_expectations,
            expected_merged_gene_lengths,
            "query-to-all BAM",
        ),
        (
            history_to_query_bam,
            history_resume_expectations,
            query_gene_len_dic,
            "history-to-query BAM",
        ),
    ):
        if resume_bam_path:
            resume_provenance = validate_resume_bam(
                resume_bam_path, expected_queries, expected_targets, label
            )
            # The raw BAM's aligner provenance is checked here. Filtering is
            # performed below by this PanTrans run, so its logic is known.
            resume_provenance["filter_thresholds"] = current_provenance[
                "filter_thresholds"
            ]
            resume_provenance["filter_thresholds_assumed"] = False
            resume_provenance["filter_logic_id"] = FILTER_LOGIC_ID
            validate_edge_provenance(
                {"edge_generations": [resume_provenance]}, current_provenance
            )

    logger.info("Align new cDNA against complete historical-plus-new gDNA")
    if query_to_all_bam:
        logger.info("Reuse query-to-all BAM: %s", query_to_all_bam_path)
    else:
        minimap2_map(query_cdna_path, merged_gdna_path, threads, query_to_all_bam_path)
    _, gene_len_dic = filter_bam(
        query_to_all_bam_path,
        query_to_all_filtered_bam_path,
        bed_dic,
        collect_edges=False,
    )

    logger.info("Align complete historical cDNA against new gDNA")
    if history_to_query_bam:
        logger.info("Reuse history-to-query BAM: %s", history_to_query_bam_path)
    else:
        minimap2_map(history_cdna_path, query_gdna_path, threads, history_to_query_bam_path)
    _, _ = filter_bam(
        history_to_query_bam_path,
        history_to_query_filtered_bam_path,
        bed_dic,
        collect_edges=False,
    )

    merge_bams_with_full_header(
        [
            history_package["filtered_bam_path"],
            query_to_all_filtered_bam_path,
            history_to_query_filtered_bam_path,
        ],
        query_to_all_bam_path,
        filtered_bam_path,
    )

    logger.info("Merge historical and cross-alignment edges, then assign clusters")
    graph, pre_clusters, last_clusters = _derive_graph_and_clusters(
            aligned_gene_li=chain(
                iter_graph_edges(history_package["edges_path"]),
                iter_filtered_bam_gene_edges(query_to_all_filtered_bam_path),
                iter_filtered_bam_gene_edges(history_to_query_filtered_bam_path),
            ),
        gene_len_dic=gene_len_dic,
        bed_dic=bed_dic,
        variety_li=variety_li,
        refer_name=history_package["reference_name"],
        eligible_gene_set=eligible_gene_set,
        main_chroms=history_package.get("main_chroms"),
    )

    pre_cluster_dic = cluster2dic(pre_clusters)
    last_cluster_dic = cluster2dic(last_clusters)
    last_rename_map = _build_gene_rename_map(last_cluster_dic, merged_bed_path, variety_li, "Pan")

    graph_package_path = os.path.join(out_dir, f"{prefix}.graph.json")
    append_provenance = dict(current_provenance)
    append_provenance["scope"] = "append_cross_alignments"
    write_graph_package(
        manifest_path=graph_package_path,
        edge_iter=graph.edges(),
        gene_len_dic=gene_len_dic,
        bed_path=merged_bed_path,
        filtered_bam_path=filtered_bam_path,
        cdna_paths=[history_cdna_path, query_cdna_path],
        gdna_paths=history_package["gdna_paths"] + [query_gdna_path],
        variety_names=variety_li,
        reference_name=history_package["reference_name"],
        main_chroms=history_package.get("main_chroms"),
        provenance={"edge_generations": history_generations + [append_provenance]},
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
        seed_gtf_path=history_gtf_path,
    )

    return (
        merged_cdna_path,
        merged_gdna_path,
        merged_bed_path,
        graph_package_path,
        filtered_bam_path,
    )
