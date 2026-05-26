import subprocess
import pysam

# 阈值可按需调整
COV_MIN = 0.80      # 覆盖度阈值
PID_MIN = 0.90      # 近似一致性阈值
SOFT_MAX = 0.10     # 软截断比例上限

def minimap2_map(cdna_path, gdna_path, threads, bam_path):
    # use minimap2 to align cdna to gdna and write an unsorted BAM
    cmd = f"minimap2 -ax splice:hq -uf -I 16G --secondary=yes -N 100 -t {threads} {gdna_path} {cdna_path} | samtools view -@{threads} -hb - > {bam_path}"
    subprocess.run(cmd, shell=True, check=True)

def minimap2_map_rescue(cdna_path, gdna_path, threads, bam_path):
    """Run a lighter splice-aware minimap2 pass for rescue alignments."""
    cmd = (
        f"minimap2 -ax splice:hq -uf --secondary=yes -N 100 -t 1 "
        f"{gdna_path} {cdna_path} | "
        f"samtools view -@1 -hb - > {bam_path}"
    )
    subprocess.run(cmd, shell=True, check=True)

def parse_cigar(cigartuples):
    """解析 CIGAR，返回 query 覆盖长度、对齐长度、软截断长度"""
    if not cigartuples:
        return 0, 0, 0
    mlen = 0
    aln_len = 0
    soft = 0
    for op, n in cigartuples:
        if op in (0, 7, 8):  # M, =, X：同时消耗 query 和 reference
            mlen += n
            aln_len += n
        elif op == 1:        # I：只消耗 query
            mlen += n
            aln_len += n
        elif op == 2:        # D：只消耗 reference，不应计入 query 覆盖长度
            aln_len += n
        elif op == 4:        # S
            soft += n
    return mlen, aln_len, soft

def get_tag(read, key, default=0):
    """安全获取 BAM 标签"""
    try:
        return read.get_tag(key)
    except (KeyError, ValueError):
        return default

def filter_bam(bam_path, output_bam, bed_dic):
    """
    过滤 BAM 文件并生成新的 BAM 文件。
    过滤逻辑：
    - 覆盖度 >= COV_MIN
    - PID >= PID_MIN
    - 软截断比例 <= SOFT_MAX

    参数：
        bam_path: 输入 BAM 文件路径
        output_bam: 输出 BAM 文件路径
    """
    bam_file = pysam.AlignmentFile(bam_path, "rb")
    out_bam = pysam.AlignmentFile(output_bam, "wb", header=bam_file.header)

    target_lengths = dict(zip(bam_file.references, bam_file.lengths))
    pair_gene_li = []
    # bed_dic[gene] = [chrom, start, end]
    for read in bam_file.fetch(until_eof=True):
        if read.is_unmapped or read.is_supplementary:
            continue

        refer_name = bam_file.get_reference_name(read.reference_id)
        if refer_name not in bed_dic:
            # Ignore alignments to references that are not present in BED annotation.
            continue
        query_trans = read.query_name
        query_parts = query_trans.split(".")
        query_trans = ".".join(query_parts[:-1]) if len(query_parts) > 1 else query_trans
        if query_trans == refer_name:
            out_bam.write(read)
            pair_gene_li.append((query_trans, refer_name))
            continue
        cigar = read.cigartuples
        seq = read.query_sequence
        mlen, aln_len, soft = parse_cigar(cigar)
        #if (qstart + mlen + soft) > target_lengths[refer_name]:
        #    continue

        qlen = len(seq) if seq else mlen
        cov = mlen / qlen if qlen > 0 else 0.0
        soft_ratio = soft / qlen if qlen > 0 else 0.0

        nm = get_tag(read, "NM", 0)
        pid = 1.0
        if aln_len > 0:
            pid = 1.0 - (float(nm) / float(aln_len))

        passed_criteria = sum([
            cov >= COV_MIN,
            pid >= PID_MIN,
            soft_ratio <= SOFT_MAX
        ])
        if passed_criteria == 3:
            out_bam.write(read)
            pair_gene_li.append((query_trans, refer_name))

    bam_file.close()
    out_bam.close()

    return pair_gene_li, target_lengths
