import hashlib
import os
import pysam
import subprocess
import logging

logger = logging.getLogger(__name__)


def sequence_identity(sequence):
    """Return a stable identity for one biological sequence."""
    normalized = sequence.upper()
    return {
        "length": len(normalized),
        "sha256": hashlib.sha256(normalized.encode("ascii")).hexdigest(),
    }


def get_fasta_sequence_identities(fasta_path):
    """Return sequence identities keyed by the first token of each FASTA header."""
    identities = {}
    record_id = None
    seq_parts = []
    with open(fasta_path, "rt", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if record_id is not None:
                    identities[record_id] = sequence_identity("".join(seq_parts))
                record_id = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)
    if record_id is not None:
        identities[record_id] = sequence_identity("".join(seq_parts))
    return identities


def get_fasta_len(fasta_path):
    """读取FASTA文件，计算每个基因的总长度"""
    gene_len_dic = {}
    gene_n = None
    seq_parts = []
    with open(fasta_path, "rt") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                # 保存上一个基因的序列长度
                if gene_n is not None:
                    gene_len_dic[gene_n] = sum(len(s) for s in seq_parts)
                # 提取基因名（只取第一个空格前的部分）
                gene_n = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)
        # 保存最后一个基因
        if gene_n is not None:
            gene_len_dic[gene_n] = sum(len(s) for s in seq_parts)
    return gene_len_dic

def get_bed(bed_path):
    """读取BED文件，返回基因坐标字典"""
    bed_dic = {}
    gene_strand_dic = {}
    with open(bed_path, "rt") as f:
        for line in f:
            li = line.strip().split("\t")
            gene_strand_dic[li[3]] = li[5]
            chrom, start, end, gene = li[0], int(li[1]), int(li[2]), li[3]
            # bed_dic 始终保存基因组坐标，链方向由 gene_strand_dic 单独维护。
            bed_dic[gene] = [chrom, min(start, end), max(start, end)]
    return bed_dic, gene_strand_dic

def get_bed_rows(bed_path):
    """读取 BED 文件并保留原始顺序。"""
    rows = []
    with open(bed_path, "rt") as handle:
        for line in handle:
            li = line.strip().split("\t")
            if len(li) < 6:
                continue
            rows.append((li[0], int(li[1]), int(li[2]), li[3], li[5]))
    return rows

def cluster2dic(cluster_li):
    """将簇列表转换为基因-簇字典。

    支持两种簇表示形式：
    - List[str]：例如 ['gene1', 'gene2', 'gene3']
    - str：例如 'gene1\\tgene2\\tgene3'
    """
    cluster_dic = {}
    for line in cluster_li:
        if not line:
            continue
        # 如果是字符串，按制表符拆分；如果已经是列表，则直接使用
        if isinstance(line, str):
            nodes = line.rstrip().split("\t")
        else:
            nodes = list(line)
        if not nodes:
            continue
        hub = nodes[0]
        others = nodes[:]
        cluster_dic[hub] = others
    return cluster_dic

def extract_fasta_subset_by_names(
    fasta_file: str,
    gene_names: list,
    output_file: str,
    match_mode: str = "exact",         # "exact" 或 "token"
    case_sensitive: bool = False,
    header_style: str = "original",        # "gene", "original", "both"
    line_width: int = 60
):
    """
    根据基因名称列表从 fasta 中提取同名序列子集（不使用坐标）。
    如果 .fai 文件不存在，则自动创建。
    """

    # 检查并生成 fai 索引
    fai_file = fasta_file + ".fai"
    if not os.path.exists(fai_file):
        logger.info(f"未找到索引文件 {fai_file}，正在创建...")
        subprocess.run(["samtools", "faidx", fasta_file], check=True)
        logger.info("索引文件已生成。")

    # 读取基因名列表

    def norm(s: str) -> str:
        return s if case_sensitive else s.lower()

    def head_token(h: str) -> str:
        t = h.split(" ", 1)[0]
        t = t.split("|", 1)[0]
        return t

    fa = pysam.FastaFile(fasta_file)

    headers = list(fa.references)
    if match_mode == "exact":
        key_map = {norm(h): h for h in headers}
    elif match_mode == "token":
        key_map = {norm(head_token(h)): h for h in headers}
    else:
        fa.close()
        raise ValueError("match_mode 仅支持 'exact' 或 'token'。")

    with open(output_file, "w") as out:
        for g in gene_names:
            gkey = norm(g if match_mode == "exact" else head_token(g))
            if gkey in key_map:
                h = key_map[gkey]
                seq = fa.fetch(h)
                if header_style == "gene":
                    header = f">{g}"
                elif header_style == "original":
                    header = f">{h}"
                elif header_style == "both":
                    header = f">{g}|{h}"
                else:
                    fa.close()
                    raise ValueError("header_style 仅支持 'gene'、'original'、'both'。")
                out.write(header + "\n")
                for i in range(0, len(seq), line_width):
                    out.write(seq[i:i+line_width] + "\n")
    fa.close()

def concat_fasta_files(input_paths, output_path):
    """按顺序合并多个 FASTA 文件。"""
    with open(output_path, "w") as out_handle:
        for path in input_paths:
            with open(path, "rt") as in_handle:
                for line in in_handle:
                    out_handle.write(line)
    return output_path

def concat_text_files(input_paths, output_path):
    """按顺序合并多个文本文件。"""
    with open(output_path, "w") as out_handle:
        for path in input_paths:
            with open(path, "rt") as in_handle:
                for line in in_handle:
                    out_handle.write(line)
    return output_path

def extract_transcripts_by_gene_names(fasta_file, gene_names, output_file, line_width=60):
    """从 transcript FASTA 中提取指定 gene 的所有 transcript。"""
    gene_name_set = set(gene_names)
    kept = 0
    with open(fasta_file, "rt") as in_handle, open(output_file, "w") as out_handle:
        keep_record = False
        header = None
        seq_parts = []
        for raw_line in in_handle:
            line = raw_line.rstrip("\n")
            if line.startswith(">"):
                if header is not None and keep_record:
                    out_handle.write(header + "\n")
                    seq = "".join(seq_parts)
                    for i in range(0, len(seq), line_width):
                        out_handle.write(seq[i:i + line_width] + "\n")
                    kept += 1
                transcript_id = line[1:].split()[0]
                gene_id = ".".join(transcript_id.split(".")[:-1]) if "." in transcript_id else transcript_id
                keep_record = gene_id in gene_name_set
                header = line
                seq_parts = []
            else:
                if header is not None:
                    seq_parts.append(line.strip())
        if header is not None and keep_record:
            out_handle.write(header + "\n")
            seq = "".join(seq_parts)
            for i in range(0, len(seq), line_width):
                out_handle.write(seq[i:i + line_width] + "\n")
            kept += 1
    return kept

def extract_fasta_records_by_exact_names(fasta_file, record_names, output_file, line_width=60):
    """从 FASTA 中按精确 header token 提取记录，不依赖 faidx。"""
    record_name_set = set(record_names)
    kept = 0
    with open(fasta_file, "rt") as in_handle, open(output_file, "w") as out_handle:
        keep_record = False
        header = None
        seq_parts = []
        for raw_line in in_handle:
            line = raw_line.rstrip("\n")
            if line.startswith(">"):
                if header is not None and keep_record:
                    out_handle.write(header + "\n")
                    seq = "".join(seq_parts)
                    for i in range(0, len(seq), line_width):
                        out_handle.write(seq[i:i + line_width] + "\n")
                    kept += 1
                record_id = line[1:].split()[0]
                keep_record = record_id in record_name_set
                header = line
                seq_parts = []
            else:
                if header is not None:
                    seq_parts.append(line.strip())
        if header is not None and keep_record:
            out_handle.write(header + "\n")
            seq = "".join(seq_parts)
            for i in range(0, len(seq), line_width):
                out_handle.write(seq[i:i + line_width] + "\n")
            kept += 1
    return kept

def build_transcript_index(fasta_file):
    """从 transcript FASTA 建立 gene_id -> [(transcript_id, header, seq)] 索引。"""
    transcript_index = {}
    with open(fasta_file, "rt") as in_handle:
        header = None
        transcript_id = None
        seq_parts = []
        for raw_line in in_handle:
            line = raw_line.rstrip("\n")
            if line.startswith(">"):
                if header is not None and transcript_id is not None:
                    gene_id = ".".join(transcript_id.split(".")[:-1]) if "." in transcript_id else transcript_id
                    transcript_index.setdefault(gene_id, []).append((transcript_id, header, "".join(seq_parts)))
                header = line
                transcript_id = line[1:].split()[0]
                seq_parts = []
            else:
                if header is not None:
                    seq_parts.append(line.strip())
        if header is not None and transcript_id is not None:
            gene_id = ".".join(transcript_id.split(".")[:-1]) if "." in transcript_id else transcript_id
            transcript_index.setdefault(gene_id, []).append((transcript_id, header, "".join(seq_parts)))
    return transcript_index

def write_transcripts_from_index(transcript_index, gene_names, output_file, line_width=60):
    """根据 gene_id -> transcript 索引写出指定 gene 的所有 transcript。"""
    kept = 0
    with open(output_file, "w") as out_handle:
        for gene_id in gene_names:
            for transcript_id, header, seq in transcript_index.get(gene_id, []):
                out_handle.write(header + "\n")
                for i in range(0, len(seq), line_width):
                    out_handle.write(seq[i:i + line_width] + "\n")
                kept += 1
    return kept

def extract_fasta_records_by_exact_names_indexed(fasta_file, record_names, output_file, line_width=60):
    """用 faidx 按精确记录名提取 FASTA，避免全文件扫描。"""
    record_name_set = set(record_names)
    fai_file = fasta_file + ".fai"
    if not os.path.exists(fai_file):
        logger.info(f"未找到索引文件 {fai_file}，正在创建...")
        subprocess.run(["samtools", "faidx", fasta_file], check=True)
        logger.info("索引文件已生成。")
    fa = pysam.FastaFile(fasta_file)
    kept = 0
    with open(output_file, "w") as out_handle:
        for record_id in fa.references:
            if record_id not in record_name_set:
                continue
            seq = fa.fetch(record_id)
            out_handle.write(f">{record_id}\n")
            for i in range(0, len(seq), line_width):
                out_handle.write(seq[i:i + line_width] + "\n")
            kept += 1
    fa.close()
    return kept
