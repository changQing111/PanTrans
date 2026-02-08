import os
import sys
import pysam
import subprocess
import logging

logger = logging.getLogger(__name__)


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
            chrom, start, end, gene, strand = li[0], li[1], li[2], li[3], li[-1]
            if strand == "+":
                bed_dic[gene] = [chrom, start, end]
            else:
                bed_dic[gene] = [chrom, end, start] # 反向链时，起始和终止互换
    return bed_dic, gene_strand_dic

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

def get_longest_pep(fasta_file, output_file=None):
    """
    从FASTA文件中提取每个基因的最长转录本
    
    参数:
        fasta_file (str): 输入FASTA文件路径
        output_file (str): 可选，输出FASTA文件路径
    
    返回:
        dict: key 为基因 ID, value 为 (转录本ID, 序列)
    """
    longest = {}
    with open(fasta_file) as f:
        transcript_id = None
        seq_parts = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                # 保存上一个转录本
                if transcript_id:
                    seq = "".join(seq_parts)
                    gene_id = transcript_id.split('.')[0]
                    if gene_id not in longest or len(seq) > len(longest[gene_id][1]):
                        longest[gene_id] = (transcript_id, seq)
                # 新转录本开始
                transcript_id = line[1:].split()[0]  # 取ID
                seq_parts = []
            else:
                seq_parts.append(line)
        # 最后一条
        if transcript_id:
            seq = "".join(seq_parts)
            gene_id = transcript_id.split('.')[0]
            if gene_id not in longest or len(seq) > len(longest[gene_id][1]):
                longest[gene_id] = (transcript_id, seq)

    # 如果指定了输出文件路径，则写入结果
    if output_file:
        with open(output_file, "w") as out:
            for gene_id, (tid, seq) in longest.items():
                out.write(f">{tid}\n")
                # 按 60 字符换行，符合FASTA规范
                for i in range(0, len(seq), 60):
                    out.write(seq[i:i+60] + "\n")

    return longest


