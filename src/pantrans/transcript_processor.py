import pysam
import Bio.SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import logging

from .align_filter import transcript_to_gene_id
from .gff3 import (load_gff3_transcript_models, load_history_transcript_models,
                   make_row, write_gff3, rename_gff3_ids, sort_gff3_by_gene_id)


logger = logging.getLogger(__name__)


def _get_tag(read, key, default=0):
    try:
        return read.get_tag(key)
    except (KeyError, ValueError):
        return default


def _alignment_rank(read):
    match_len = 0
    for op, length in read.cigartuples or []:
        if op in (0, 7, 8):
            match_len += length
    score = _get_tag(read, "AS", _get_tag(read, "ms", 0))
    nm = _get_tag(read, "NM", 0)
    return (not read.is_secondary, score, match_len, -nm)


def _extract_exon_and_splice(read):
    splice_site_li = []
    exon_coord_li = []
    current_pos = read.reference_start  # BAM reference 坐标是 0-based 半开区间
    exon_start = current_pos

    for op, length in read.cigartuples:
        if op in (0, 2, 7, 8):  # M, D, =, X 消耗 reference 坐标
            current_pos += length
        elif op == 3:  # N 表示 intron / splice gap
            if current_pos > exon_start:
                exon_coord_li.append((exon_start + 1, current_pos))

            splice_start = current_pos + 1
            splice_end = current_pos + length
            splice_site_li.append((splice_start, splice_end))

            current_pos = splice_end
            exon_start = current_pos
        else:
            # I/S/H/P 等不消耗 reference 坐标
            pass

    if current_pos > exon_start:
        exon_coord_li.append((exon_start + 1, current_pos))

    return splice_site_li, exon_coord_li

def extract_exon_coord_splice_site(bam_path, cluster_dic, include_strands=False):
    """
    提取转录本的剪切位点和外显子坐标信息
    Args:
        bam_path (str): BAM文件路径
        cluster_dic (dict): {reference_name: [gene_ids]} 映射关系
    Returns:
        trans_splice_site_dic (dict): {read_id: [(splice_start, splice_end), ...]}
        exon_coord_dic (dict): {read_id: [(exon_start, exon_end), ...]}
    """
    best_alignment_dic = {}

    bam = pysam.AlignmentFile(bam_path, "rb")
    for read in bam:
        if read.is_unmapped or read.is_supplementary:
            continue

        refer = bam.get_reference_name(read.reference_id)
        query = read.query_name
        if refer not in cluster_dic:
            continue

        gene_li = cluster_dic[refer]
        query_parts = query.split(".")
        query_id = ".".join(query_parts[:-1]) if len(query_parts) > 1 else query
        if query_id not in gene_li:
            continue

        splice_site_li, exon_coord_li = _extract_exon_and_splice(read)
        if not exon_coord_li:
            continue
        rank = _alignment_rank(read)
        if query not in best_alignment_dic or rank > best_alignment_dic[query][0]:
            best_alignment_dic[query] = (rank, splice_site_li, exon_coord_li, "-" if read.is_reverse else "+")

    bam.close()
    trans_splice_site_dic = {
        query: splice_site_li
        for query, (_, splice_site_li, _, _) in best_alignment_dic.items()
    }
    exon_coord_dic = {
        query: exon_coord_li
        for query, (_, _, exon_coord_li, _) in best_alignment_dic.items()
    }
    if include_strands:
        return trans_splice_site_dic, exon_coord_dic, {q: value[3] for q, value in best_alignment_dic.items()}
    return trans_splice_site_dic, exon_coord_dic

def group_by_longest_superstring(orig: dict) -> dict:
    """
    将字典按 value 的子串关系分组，使用每组中最长的 value 作为输出字典的 key，
    对应 value 为属于该组的原始 keys 列表。
    """
    # value -> list of keys
    val_to_keys = {}
    for k, v in orig.items():
        val_to_keys.setdefault(v, []).append(k)

    # 所有不同的 value，按长度降序
    values = sorted(val_to_keys.keys(), key=len, reverse=True)

    assigned = set()   # 已经被分配到某组的 value
    result = {}

    for rep in values:
        if rep in assigned:
            continue
        group_keys = []
        # 找出所有尚未分配且是 rep 的子串（包含相等）
        for v in list(val_to_keys.keys()):
            if v in assigned:
                continue
            if v in rep:  # v 是 rep 的子串或相等
                group_keys.extend(val_to_keys[v])
                assigned.add(v)
        if group_keys:
            result[rep] = group_keys

    return result

def dedup_ss(splice_site_dic, trans_li, transcript_strands=None):
    """根据剪切位点对转录本进行去重"""
    trans_intron_coord_dic = {}
    for trans in trans_li:
        # 跳过不在 splice_site_dic 中的转录本（可能在某些情况下缺失）
        if trans not in splice_site_dic:
            continue
        coord_c = splice_site_dic[trans]
    #    print(coord_c)
        if len(coord_c) == 0:
            coord_str = "mono"
        else:
            coord_str = str(coord_c[0][0]) + "-" + str(coord_c[0][1])
        if len(coord_c) > 1:
            for coord in coord_c[1:]:
                tmp = str(coord[0]) + "-" + str(coord[1])
                coord_str = coord_str + ":" + tmp
        trans_intron_coord_dic[trans] = coord_str
    if transcript_strands is None:
        return group_by_longest_superstring(trans_intron_coord_dic)
    groups = {}
    for strand in ("+", "-"):
        subset = {t: pattern for t, pattern in trans_intron_coord_dic.items()
                  if transcript_strands[t] == strand}
        for pattern, members in group_by_longest_superstring(subset).items():
            groups[(strand, pattern)] = members
    return groups

def get_last_trans(splice_site_dic, cluster_dic, gene_trans_dic, transcript_strands=None):
    '''out gene transcript map'''
    # Abo6B402900.1 [(2229, 2339), (2404, 2486), (3760, 3842), (4548, 5809)]
    gene_ss_dic = {}
    for gene, gene_li in cluster_dic.items():
        trans_li = []
        for i in gene_li:
            if i in gene_trans_dic:
                trans_li += gene_trans_dic[i] # [.1, .2]
        gene_ss_dic[gene] = dedup_ss(splice_site_dic, trans_li, transcript_strands)
    return gene_ss_dic

# JM471A001600 {'2811-2915:2957-3221:3289-3367:3410-3511:3543-3922:4093-4222:4400-4493:4691-4780:4989-5086': ['JM471A001600.1', 'ZM22Ctg10006500.1'], '2811-2915:2957-3221:3289-3367:3410-3511:3543-3922:4093-4222': ['Abo1A000300.1']}
def generate_gff3(gene_splice_site_dic, exon_coord_dic, gene_len_dic, trans_len_dic, gene_strand_dic, rename_map, gff3_path, transcript_strands=None):
    """Write gene/mRNA/exon hierarchies for the longest model of each splice pattern."""
    rows = []
    for gene_id, groups in gene_splice_site_dic.items():
        if gene_id not in gene_len_dic:
            continue
        length = gene_len_dic[gene_id]
        output_gene = rename_map.get(gene_id, gene_id) if rename_map else gene_id
        children = []
        for index, candidates in enumerate(groups.values(), 1):
            valid = [t for t in candidates if t in trans_len_dic and exon_coord_dic.get(t)]
            if not valid:
                continue
            representative = max(valid, key=lambda t: trans_len_dic[t])
            if transcript_strands is None or transcript_strands.get(representative) not in {"+", "-"}:
                raise ValueError(f"Missing reference-relative strand for {representative}; BED strand is not a substitute")
            strand = transcript_strands[representative]
            exons = [(max(1, a), min(b, length)) for a, b in sorted(exon_coord_dic[representative])
                     if min(b, length) >= max(1, a)]
            if not exons:
                continue
            tid = f"{output_gene}.{index}"
            children.append(make_row(gene_id, "mRNA", exons[0][0], exons[-1][1], strand, tid, output_gene))
            for exon_index, (a, b) in enumerate(exons, 1):
                children.append(make_row(gene_id, "exon", a, b, strand, f"{tid}.exon{exon_index}", tid))
        if children:
            strands = {row[6] for row in children}
            strand = next(iter(strands)) if len(strands) == 1 else "."
            rows.append(make_row(gene_id, "gene", min(int(r[3]) for r in children),
                                 max(int(r[4]) for r in children), strand, output_gene))
            rows.extend(children)
    write_gff3(gff3_path, rows)


def get_gene_trans_dic(splice_site_dic):
    """从转录本剪切位点字典构建基因到转录本的映射"""
    gene_trans_dic = {}
    for i in splice_site_dic.keys():
        gene_id = ".".join((i.split("."))[:-1])
        if gene_id not in gene_trans_dic:
            gene_trans_dic[gene_id] = [i]
        else:
            gene_trans_dic[gene_id].append(i)
    return gene_trans_dic

def transcript_dedup(
    bam_path,
    cluster_dic,
    trans_len_dic,
    gene_len_dic,
    gene_strand_dic,
    rename_map,
    gff3_path,
    seed_gff3_path=None,
):
    """基于剪切位点进行转录本去重并生成GFF3文件"""
    trans_splice_site_dic, exon_coord_dic, transcript_strands = extract_exon_coord_splice_site(
        bam_path, cluster_dic, include_strands=True)

    if seed_gff3_path is not None:
        seed_models = load_history_transcript_models(seed_gff3_path)
        logger.warning("Historical annotation strands are reused relative to their FASTA. "
                       "Regenerate histories produced before the strand fix; chromosome BED strands cannot correct them.")
        skipped_transcripts = [
            transcript_id
            for transcript_id, gene_id in seed_models["transcript_gene"].items()
            if gene_id not in cluster_dic
        ]
        skipped_genes = sorted(
            {
                seed_models["transcript_gene"][transcript_id]
                for transcript_id in skipped_transcripts
            }
        )
        if skipped_transcripts:
            logger.info(
                "Skipped %d historical GFF3 transcript%s from %d non-current "
                "representative gene%s: %s",
                len(skipped_transcripts),
                "" if len(skipped_transcripts) == 1 else "s",
                len(skipped_genes),
                "" if len(skipped_genes) == 1 else "s",
                ", ".join(skipped_genes),
            )

        for transcript_id, gene_id in seed_models["transcript_gene"].items():
            if gene_id not in cluster_dic:
                continue
            strand = seed_models["transcript_context"][transcript_id][1]
            if strand not in {"+", "-"}:
                raise ValueError(f"Unknown historical transcript strand: {transcript_id}")
            transcript_strands[transcript_id] = strand
            trans_splice_site_dic[transcript_id] = seed_models["splice_sites"][
                transcript_id
            ]
            exon_coord_dic[transcript_id] = seed_models["exon_coords"][
                transcript_id
            ]

    gene_trans_dic = get_gene_trans_dic(trans_splice_site_dic)
    gene_splice_site_dic = get_last_trans(trans_splice_site_dic, cluster_dic, gene_trans_dic, transcript_strands)
    generate_gff3(gene_splice_site_dic, exon_coord_dic, gene_len_dic, trans_len_dic, gene_strand_dic, rename_map, gff3_path, transcript_strands)

def get_subset_bed(gene_li, merge_bed_path, out_bed_path):
    merge_bed_dic = {}
    with open(merge_bed_path, "rt") as merge_bed:
        for line in merge_bed:
            li = line.rstrip().split("\t")
            merge_bed_dic[li[3]] = line
    with open(out_bed_path, "w") as out_bed:
        for i in gene_li:
            if i in merge_bed_dic:
                out_bed.write(merge_bed_dic[i])


def get_cdna_from_gff3(fasta_path, gff3_path, out_path):
    """Extract spliced cDNA using Parent relationships and unchanged sequence IDs."""
    genome = Bio.SeqIO.to_dict(Bio.SeqIO.parse(fasta_path, "fasta"))
    models = load_gff3_transcript_models(gff3_path, require_unrenamed=False)
    records = []
    for tid, exons in models["exon_coords"].items():
        seqid, strand = models["transcript_context"][tid]
        if strand not in {"+", "-"}:
            raise ValueError(f"Cannot orient cDNA for {tid}: GFF3 strand must be + or -")
        if seqid not in genome:
            raise ValueError(f"GFF3 sequence {seqid} is missing from FASTA")
        if any(end > len(genome[seqid]) for _, end in exons):
            raise ValueError(f"GFF3 exon outside FASTA sequence for {tid}")
        cdna = Seq("").join([genome[seqid].seq[a - 1:b] for a, b in exons])
        if strand == "-":
            cdna = cdna.reverse_complement()
        records.append(SeqRecord(cdna, id=tid, description="cDNA"))
    Bio.SeqIO.write(records, out_path, "fasta")
    return len(records)
