import re
import os
import subprocess
import pysam

# 阈值可按需调整
COV_MIN = 0.80      # 覆盖度阈值
PID_MIN = 0.90      # 近似一致性阈值
SOFT_MAX = 0.10     # 软截断比例上限

def diamond_makedb(pep_path, threads, db_n):
    cmd = f"diamond makedb --threads {threads} --in {pep_path} --db {db_n}"
    subprocess.run(cmd, shell=True, check=True)

def diamond_blastp(db_n, query_path, threads, blastp_res_path):
    # use diamond blastp self all.vs. self all
    cmd = f"diamond blastp --threads {threads} --db {db_n} --query {query_path} --max-target-seqs 100 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen > {blastp_res_path}"
    subprocess.run(cmd, shell=True, check=True)
    
def parse_blast_line(line):
    """
    解析 DIAMOND/BLAST outfmt 6 格式的一行结果
    """
    fields = line.strip().split("\t")
    if len(fields) < 14:
        raise ValueError(f"BLAST outfmt 6 格式至少需要14个字段，但当前行只有 {len(fields)} 个字段: {line[:100]}")
    return {
        "qseqid": fields[0],
        "sseqid": fields[1],
        "pident": float(fields[2]),
        "length": int(fields[3]),
        "qstart": int(fields[6]),
        "qend": int(fields[7]),
        "evalue": float(fields[10]),
        "bitscore": float(fields[11]),
        "qlen": int(fields[12]),
        "slen": int(fields[13])
    }

def filter_blastp(infile, outfile,
                evalue_thresh=1e-20,
                identity_thresh=70.0,
                qcov_thresh=0.7,
                scov_thresh=0.7,
                best_hit_only=True):
    """
    过滤 BLAST/DIAMOND 结果
    参数:
      evalue_thresh   - E-value 阈值
      identity_thresh - 百分比身份率阈值
      qcov_thresh     - query 覆盖度阈值 (alignment length / qlen)
      scov_thresh     - subject 覆盖度阈值 (alignment length / slen)
      best_hit_only   - 是否只保留每个 query 的最佳 hit
    """
    hits_by_query = {}

    with open(infile) as fin:
        for line in fin:
            if line.startswith("#") or not line.strip():
                continue
            try:
                hit = parse_blast_line(line)
            except (ValueError, IndexError) as e:
                # 跳过格式错误的行，记录警告
                import warnings
                warnings.warn(f"跳过格式错误的行: {e}")
                continue

            # 覆盖度计算
            qcov = hit["length"] / hit["qlen"] if hit["qlen"] > 0 else 0
            scov = hit["length"] / hit["slen"] if hit["slen"] > 0 else 0

            if (hit["evalue"] <= evalue_thresh and
                hit["pident"] >= identity_thresh and
                qcov >= qcov_thresh and
                scov >= scov_thresh):

                if best_hit_only:
                    # 按 bitscore 保留最佳 hit
                    if hit["qseqid"] not in hits_by_query or hit["bitscore"] > hits_by_query[hit["qseqid"]]["bitscore"]:
                        hits_by_query[hit["qseqid"]] = {"line": line, "bitscore": hit["bitscore"]}
                else:
                    if hit["qseqid"] not in hits_by_query:
                        hits_by_query[hit["qseqid"]] = []
                    hits_by_query[hit["qseqid"]].append({"line": line, "bitscore": hit["bitscore"]})

    # 输出结果
    with open(outfile, "w") as fout:
        if best_hit_only:
            for q in hits_by_query:
                fout.write(hits_by_query[q]["line"])
        else:
            for q in hits_by_query:
                for h in hits_by_query[q]:
                    fout.write(h["line"])


def minimap2_map(cdna_path, gdna_path, threads, bam_path):
    # use minimap2 to align cdna to gdna, bam sorted by coordinate
    cmd = f"minimap2 -ax splice:hq -uf -I 16G --secondary=yes -N 100 --MD --cs=long -t {threads} {gdna_path} {cdna_path} | samtools view -@{threads} -hb - | samtools sort -@{threads} > {bam_path}"
    subprocess.run(cmd, shell=True, check=True)

def bam_index(bam_path):
    cmd_index = f"samtools index {bam_path}"
    try:
        subprocess.run(cmd_index, shell=True, check=True)
        return
    except (KeyError, ValueError):
        cmd_index = f"samtools index -c {bam_path}"
        subprocess.run(cmd_index, shell=True, check=True)


def parse_cigar(cigartuples):
    """解析 CIGAR，返回匹配长度、对齐长度、软截断长度"""
    mlen = 0
    aln_len = 0
    soft = 0
    for op, n in cigartuples:
        if op in (0, 7, 8):  # M, =, X
            mlen += n
            aln_len += n
        elif op in (1, 2):   # I, D
            mlen += n
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
        if read.is_unmapped:
            continue

        refer_name = bam_file.get_reference_name(read.reference_id)
        refer_coord = bed_dic[refer_name][1:]
        refer_start, refer_end = min(refer_coord), max(refer_coord)
        query_trans = read.query_name
        query_trans = ".".join(query_trans.split(".")[:-1])
        if query_trans == refer_name:
            out_bam.write(read)
            pair_gene_li.append((query_trans, refer_name))
            continue
        qstart = read.reference_start
        cigar = read.cigartuples
        seq = read.query_sequence
        qstart + ((refer_end - refer_start) / 4 * 3)
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
