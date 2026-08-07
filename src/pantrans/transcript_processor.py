import pysam
import Bio.SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import logging
import os
import re

from .align_filter import transcript_to_gene_id


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

def _parse_gtf_attributes(attrs):
    parsed = {}
    for item in attrs.strip().split(";"):
        item = item.strip()
        if not item or " " not in item:
            continue
        key, value = item.replace('"', "").split(" ", 1)
        parsed[key] = value
    return parsed


def load_gtf_transcript_models(gtf_path):
    """Load validated transcript exon and splice models from an unrenamed GTF."""
    transcript_gene = {}
    transcript_context = {}
    exon_coords = {}
    transcript_features = set()

    with open(gtf_path, "rt") as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.rstrip("\n")
            if not stripped or stripped.startswith("#"):
                continue

            fields = stripped.split("\t")
            if len(fields) != 9:
                raise ValueError(
                    f"{gtf_path}:{line_number}: expected 9 tab-separated GTF fields"
                )

            sequence_name, _, feature, start, end, _, strand, _, attr_text = fields
            if feature not in ("transcript", "exon"):
                continue

            attrs = _parse_gtf_attributes(attr_text)
            missing_attrs = [
                key for key in ("gene_id", "transcript_id") if not attrs.get(key)
            ]
            if missing_attrs:
                raise ValueError(
                    f"{gtf_path}:{line_number}: {feature} is missing "
                    + ", ".join(missing_attrs)
                )

            gene_id = attrs["gene_id"]
            transcript_id = attrs["transcript_id"]
            if sequence_name != gene_id:
                raise ValueError(
                    f"{gtf_path}:{line_number}: sequence name {sequence_name!r} "
                    f"does not match gene_id {gene_id!r}"
                )
            derived_gene_id = transcript_to_gene_id(transcript_id)
            if derived_gene_id != gene_id:
                raise ValueError(
                    f"{gtf_path}:{line_number}: transcript_id {transcript_id!r} "
                    f"maps to gene {derived_gene_id!r}, but gene_id is {gene_id!r}"
                )

            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError as exc:
                raise ValueError(
                    f"{gtf_path}:{line_number}: invalid {feature} coordinates "
                    f"{start!r}-{end!r}"
                ) from exc
            if start_i < 1 or end_i < start_i:
                raise ValueError(
                    f"{gtf_path}:{line_number}: invalid {feature} coordinates "
                    f"{start_i}-{end_i}"
                )

            context = (sequence_name, strand)
            previous_context = transcript_context.setdefault(transcript_id, context)
            previous_gene = transcript_gene.setdefault(transcript_id, gene_id)
            if previous_context != context or previous_gene != gene_id:
                raise ValueError(
                    f"{gtf_path}:{line_number}: transcript {transcript_id!r} has "
                    "mixed sequence names or strands"
                )

            if feature == "exon":
                exon_coords.setdefault(transcript_id, []).append((start_i, end_i))
            else:
                transcript_features.add(transcript_id)

    missing_transcript_features = sorted(set(transcript_gene) - transcript_features)
    if missing_transcript_features:
        raise ValueError(
            f"{gtf_path}: transcript {missing_transcript_features[0]!r} has no "
            "transcript feature"
        )

    splice_sites = {}
    sorted_exon_coords = {}
    for transcript_id in transcript_gene:
        coords = sorted(exon_coords.get(transcript_id, []))
        if not coords:
            raise ValueError(
                f"{gtf_path}: transcript {transcript_id!r} has no exon records"
            )

        transcript_splice_sites = []
        for previous_exon, next_exon in zip(coords, coords[1:]):
            if next_exon[0] <= previous_exon[1] + 1:
                raise ValueError(
                    f"{gtf_path}: transcript {transcript_id!r} has overlapping "
                    "or adjacent exons"
                )
            transcript_splice_sites.append(
                (previous_exon[1] + 1, next_exon[0] - 1)
            )

        sorted_exon_coords[transcript_id] = coords
        splice_sites[transcript_id] = transcript_splice_sites

    return {
        "transcript_gene": transcript_gene,
        "splice_sites": splice_sites,
        "exon_coords": sorted_exon_coords,
    }


def _replace_gtf_attribute_values(attr_text, replacements):
    rewritten = []
    seen = set()
    for raw_item in attr_text.strip().split(";"):
        item = raw_item.strip()
        if not item:
            continue
        if " " not in item:
            raise ValueError(f"Malformed GTF attribute {item!r}")
        key, _ = item.split(None, 1)
        if key in replacements:
            item = f'{key} "{replacements[key]}"'
            seen.add(key)
        rewritten.append(item)
    missing = sorted(set(replacements) - seen)
    if missing:
        raise ValueError("GTF record is missing " + ", ".join(missing))
    return "; ".join(rewritten) + ";"


def rename_gtf_ids(source_path, output_path, rename_map):
    """Rename GTF gene/transcript attributes while keeping source coordinates."""
    output_path = str(output_path)
    temporary_path = output_path + ".tmp"
    try:
        with open(source_path, "rt") as source, open(temporary_path, "wt") as output:
            for line_number, line in enumerate(source, start=1):
                stripped = line.rstrip("\n")
                if not stripped or stripped.startswith("#"):
                    output.write(line)
                    continue

                fields = stripped.split("\t")
                if len(fields) != 9:
                    raise ValueError(
                        f"{source_path}:{line_number}: expected 9 tab-separated GTF fields"
                    )
                attrs = _parse_gtf_attributes(fields[8])
                gene_id = attrs.get("gene_id")
                transcript_id = attrs.get("transcript_id")
                if not gene_id or not transcript_id:
                    raise ValueError(
                        f"{source_path}:{line_number}: record is missing gene_id or transcript_id"
                    )
                if fields[0] != gene_id:
                    raise ValueError(
                        f"{source_path}:{line_number}: sequence name {fields[0]!r} "
                        f"does not match gene_id {gene_id!r}"
                    )
                prefix, separator, suffix = transcript_id.rpartition(".")
                if not separator or prefix != gene_id or not suffix.isdigit():
                    raise ValueError(
                        f"{source_path}:{line_number}: transcript_id {transcript_id!r} "
                        f"must use gene_id {gene_id!r} plus a numeric suffix"
                    )
                if gene_id not in rename_map:
                    raise ValueError(
                        f"{source_path}:{line_number}: no renamed ID for gene {gene_id!r}"
                    )

                renamed_gene = rename_map[gene_id]
                fields[8] = _replace_gtf_attribute_values(
                    fields[8],
                    {
                        "gene_id": renamed_gene,
                        "transcript_id": f"{renamed_gene}.{suffix}",
                    },
                )
                output.write("\t".join(fields) + "\n")
        os.replace(temporary_path, output_path)
    except Exception:
        if os.path.exists(temporary_path):
            os.unlink(temporary_path)
        raise

def _pan_gene_sort_key(gene_id):
    if not gene_id.startswith("Pan"):
        return (999, 999, gene_id, 999999999)
    core = gene_id[3:]
    contig_match = re.match(r"^Ctg(\d+)$", core)
    if contig_match:
        return (999, 999, "Ctg", int(contig_match.group(1)))
    chrom_match = re.match(r"^(\d+)([ABD])(\d+)$", core)
    if chrom_match:
        chrom_num = int(chrom_match.group(1))
        subgenome = {"A": 0, "B": 1, "D": 2}[chrom_match.group(2)]
        order_num = int(chrom_match.group(3))
        return (chrom_num, subgenome, chrom_match.group(2), order_num)
    return (999, 999, core, 999999999)

def sort_gtf_by_gene_id(gtf_path):
    """按 Pan gene_id 顺序重排 GTF。"""
    records = []
    with open(gtf_path, "rt") as handle:
        for line in handle:
            stripped = line.rstrip("\n")
            if not stripped:
                continue
            fields = stripped.split("\t")
            if len(fields) != 9:
                records.append(((999, 999, "misc", 999999999), 999, 999999999, 999999999, stripped))
                continue
            attrs = _parse_gtf_attributes(fields[8])
            gene_id = attrs.get("gene_id", fields[0])
            transcript_id = attrs.get("transcript_id", "")
            transcript_index = int(transcript_id.rsplit(".", 1)[1]) if "." in transcript_id and transcript_id.rsplit(".", 1)[1].isdigit() else 1
            feature_order = 0 if fields[2] == "transcript" else 1
            start_i = int(fields[3]) if fields[3].isdigit() else 0
            end_i = int(fields[4]) if fields[4].isdigit() else 0
            records.append((_pan_gene_sort_key(gene_id), transcript_index, feature_order, start_i, end_i, stripped))
    records.sort(key=lambda item: item[:-1])
    with open(gtf_path, "w") as handle:
        for record in records:
            handle.write(record[-1] + "\n")


def extract_exon_coord_splice_site(bam_path, cluster_dic):
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
            best_alignment_dic[query] = (rank, splice_site_li, exon_coord_li)

    bam.close()
    trans_splice_site_dic = {
        query: splice_site_li
        for query, (_, splice_site_li, _) in best_alignment_dic.items()
    }
    exon_coord_dic = {
        query: exon_coord_li
        for query, (_, _, exon_coord_li) in best_alignment_dic.items()
    }
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

def dedup_ss(splice_site_dic, trans_li):
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
    ss_dic = group_by_longest_superstring(trans_intron_coord_dic)
    return ss_dic

def get_last_trans(splice_site_dic, cluster_dic, gene_trans_dic):
    '''out gene transcript map'''
    # Abo6B402900.1 [(2229, 2339), (2404, 2486), (3760, 3842), (4548, 5809)]
    gene_ss_dic = {}
    for gene, gene_li in cluster_dic.items():
        trans_li = []
        for i in gene_li:
            if i in gene_trans_dic:
                trans_li += gene_trans_dic[i] # [.1, .2]
        gene_ss_dic[gene] = dedup_ss(splice_site_dic, trans_li)
    return gene_ss_dic

# JM471A001600 {'2811-2915:2957-3221:3289-3367:3410-3511:3543-3922:4093-4222:4400-4493:4691-4780:4989-5086': ['JM471A001600.1', 'ZM22Ctg10006500.1'], '2811-2915:2957-3221:3289-3367:3410-3511:3543-3922:4093-4222': ['Abo1A000300.1']}
def generate_gtf(gene_splice_site_dic, exon_coord_dic, gene_len_dic, trans_len_dic, gene_strand_dic, rename_map, gtf_path):
    """生成GTF文件，每个剪切位点模式选择最长的转录本作为代表"""
    with open(gtf_path, "w") as gtf_file:
        for gene_id, group_dic in gene_splice_site_dic.items():
            if gene_id not in gene_strand_dic or gene_id not in gene_len_dic:
                continue
            strand = gene_strand_dic[gene_id]
            gene_len = gene_len_dic[gene_id]
            output_gene_id = rename_map.get(gene_id, gene_id) if rename_map else gene_id

            for index, trans_li in enumerate(group_dic.values(), start=1):
                valid_trans = [
                    trans
                    for trans in trans_li
                    if trans in trans_len_dic and trans in exon_coord_dic and exon_coord_dic[trans]
                ]
                if not valid_trans:
                    continue

                rep_trans = max(valid_trans, key=lambda trans_id: trans_len_dic[trans_id])
                exon_coord_li = sorted(exon_coord_dic[rep_trans], key=lambda coord: coord[0])
                trans_id = f"{output_gene_id}.{index}"
                anno_col = f'transcript_id "{trans_id}"; gene_id "{output_gene_id}";'

                start_i = max(1, exon_coord_li[0][0])
                end_i = min(exon_coord_li[-1][1], gene_len)
                if end_i < start_i:
                    continue

                gtf_trans_li = [
                    gene_id,
                    "Pan",
                    "transcript",
                    str(start_i),
                    str(end_i),
                    ".",
                    strand,
                    ".",
                    anno_col,
                ]
                gtf_file.write("\t".join(gtf_trans_li) + "\n")
                for exon_start, exon_end in exon_coord_li:
                    exon_start = max(1, exon_start)
                    exon_end = min(exon_end, gene_len)
                    if exon_end < exon_start:
                        continue
                    gtf_exon_li = [
                        gene_id,
                        "Pan",
                        "exon",
                        str(exon_start),
                        str(exon_end),
                        ".",
                        strand,
                        ".",
                        anno_col,
                    ]
                    gtf_file.write("\t".join(gtf_exon_li) + "\n")

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
    gtf_path,
    seed_gtf_path=None,
):
    """基于剪切位点进行转录本去重并生成GTF文件"""
    trans_splice_site_dic, exon_coord_dic = extract_exon_coord_splice_site(bam_path, cluster_dic)

    if seed_gtf_path is not None:
        seed_models = load_gtf_transcript_models(seed_gtf_path)
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
                "Skipped %d historical GTF transcript%s from %d non-current "
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
            trans_splice_site_dic[transcript_id] = seed_models["splice_sites"][
                transcript_id
            ]
            exon_coord_dic[transcript_id] = seed_models["exon_coords"][
                transcript_id
            ]

    gene_trans_dic = get_gene_trans_dic(trans_splice_site_dic)
    gene_splice_site_dic = get_last_trans(trans_splice_site_dic, cluster_dic, gene_trans_dic)
    generate_gtf(gene_splice_site_dic, exon_coord_dic, gene_len_dic, trans_len_dic, gene_strand_dic, rename_map, gtf_path)

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


def get_cdna_from_gtf(fasta_path, gtf_path, out_path):
    """从GTF文件中提取 cDNA 序列并写入 FASTA。

    Args:
        fasta_path (str): 基因组 FASTA 文件路径
        gtf_path (str): 注释 GTF 文件路径
        out_path (str): 输出 cDNA FASTA 路径
    Returns:
        int: 写出的转录本条目数量
    """
    genome = Bio.SeqIO.to_dict(Bio.SeqIO.parse(fasta_path, "fasta"))

    from collections import defaultdict
    transcripts = defaultdict(list)  # transcript_id -> list[(chrom,start,end,strand)]

    with open(gtf_path) as gtf:
        for line in gtf:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            chrom, _, feature, start, end, _, strand, _, attrs = fields
            if feature != "exon":
                continue

            # 提取 transcript_id
            kv = {}
            for item in attrs.strip().split(";"):
                if not item.strip():
                    continue
                if " " in item:
                    k, v = item.strip().replace('"', "").split(" ", 1)
                    kv[k] = v
            tid = kv.get("transcript_id") or kv.get("transcriptId") or kv.get("transcript")
            if not tid:
                continue

            start_i, end_i = int(start), int(end)
            transcripts[tid].append((chrom, start_i, end_i, strand))

    records = []
    for tid, exons in transcripts.items():
        if not exons:
            continue
        strand = exons[0][3]
        exons_sorted = sorted(exons, key=lambda x: x[1])

        seq_parts = []
        for chrom, s, e, _ in exons_sorted:
            if chrom not in genome:
                continue  # 跳过缺失的染色体
            seq_parts.append(genome[chrom].seq[s - 1 : e])  # GTF 是 1-based 闭区间

        if not seq_parts:
            continue

        cdna = Seq("").join(seq_parts)
        if strand == "-":
            cdna = cdna.reverse_complement()
        records.append(SeqRecord(cdna, id=tid, description="cDNA"))

    Bio.SeqIO.write(records, out_path, "fasta")
    return len(records)
