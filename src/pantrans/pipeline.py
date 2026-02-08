import sys
import os
import logging
from .common_func import get_fasta_len, get_bed, cluster2dic, extract_fasta_subset_by_names
from .align_filter import minimap2_map, filter_bam
from .construct_di_graph import *
from .transcript_processor import transcript_dedup, get_cdna_from_gtf, get_subset_bed

logger = logging.getLogger(__name__)

def unit_construct(all_cdna_path, all_gdna_path, all_bed_path, variety_li, threads, out_dir, prefix, refer_name):
    refer_name = variety_li[0]
    bed_dic, gene_strand_dic = get_bed(all_bed_path)
    bam_path = os.path.join(out_dir, f"{prefix}_cdna_align_gdna.bam")
    filter_bam_path = os.path.join(out_dir, f"{prefix}_cdna_align_gdna.filtered.bam")

    # Run minimap2 alignment first
    logger.info("Start sequenct aligment, all cDNA vs all gDNA....................")
    minimap2_map(all_cdna_path, all_gdna_path, threads, bam_path)
    logger.info(f"Finish sequenct aligment, all cDNA vs all gDNA")
    aligned_gene_li, gene_len_dic = filter_bam(bam_path, filter_bam_path)
    logger.info(f"Finish filter bam, {filter_bam_path} file generated")
    # Step 2: Build graph and derive final gene clusters
    logger.info(f"Start build graph and assign gene")
    G = di_graph_from_pair(aligned_gene_li)
    sccs = get_conn_comp(G)

    pre_clusters, last_clusters = assign_sccs(sccs, G, gene_len_dic, bed_dic, variety_li, refer_name)
    last_cluster_dic = cluster2dic(last_clusters)
    pre_cluster_dic = cluster2dic(pre_clusters)
    # write pre last clusters into file
    pre_cluster_path = os.path.join(out_dir, f"{prefix}_pre.tmp.cluster")
    with open(pre_cluster_path, "w") as f:
        for li in pre_cluster_dic.values():
            f.write("\t".join(li) + "\n")

    last_cluster_path = os.path.join(out_dir, f"{prefix}_last.tmp.cluster")
    with open(last_cluster_path, "w") as f:
        for li in last_cluster_dic.values():
            f.write("\t".join(li) + "\n")
    refer_gene_li = list(last_cluster_dic.keys())
    # Step 3: Transcripts deduplication based on splice sites
    gtf_path = os.path.join(out_dir, f"{prefix}.gtf")
    logger.info(f"Produced Pan gtf")
    transcript_dedup(
        filter_bam_path,
        cluster_dic=last_cluster_dic,
        trans_len_dic=trans_len_dic,
        gene_len_dic=gene_len_dic,
        gene_strand_dic=gene_strand_dic,
        gtf_path=gtf_path,
    )
     # write typical gene sequence into file
    new_gdna_path = os.path.join(out_dir, f"{prefix}_gdna.refer.fasta")
    new_cdna_path =  os.path.join(out_dir, f"{prefix}_cdna.refer.fasta")
    new_bed_path = os.path.join(out_dir, f"{prefix}.refer.bed")
    logger.info(f"Produced Pan refer Gene")
    extract_fasta_subset_by_names(all_gdna_path, refer_gene_li, new_gdna_path)
    logger.info(f"Produced Pan refer mRNA")
    get_cdna_from_gtf(all_gdna_path, gtf_path, new_cdna_path)
    logger.info(f"Produced Pan bed")
    get_subset_bed(refer_gene_li, all_bed_path, new_bed_path)
    logger.info(f"{variety_li[-1]} construct end")
    return new_cdna_path, new_gdna_path, new_bed_path
