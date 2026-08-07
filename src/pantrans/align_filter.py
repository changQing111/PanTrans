import subprocess
import shutil
import os
import shlex
import tempfile
import pysam
from .version import __version__

MINIMAP2_OPTIONS = "-ax splice:hq -uf -I 16G --secondary=yes -N 100"
FILTER_LOGIC_ID = "pantrans-filter-v1"
FILTERED_BAM_MARKER = "pantrans-filtered=" + FILTER_LOGIC_ID

def merge_bams_with_full_header(source_bam_paths, full_header_bam_path, output_bam_path):
    """Merge BAM records while remapping references by name to one full header."""
    with pysam.AlignmentFile(full_header_bam_path, "rb") as header_bam:
        full_references = dict(zip(header_bam.references, header_bam.lengths))
        full_header_text = str(header_bam.header)
    for source_path in source_bam_paths:
        with pysam.AlignmentFile(source_path, "rb") as source_bam:
            source_references = dict(zip(source_bam.references, source_bam.lengths))
        missing_references = sorted(set(source_references) - set(full_references))
        mismatched_lengths = sorted(
            reference_name
            for reference_name, source_length in source_references.items()
            if reference_name in full_references
            and source_length != full_references[reference_name]
        )
        if missing_references or mismatched_lengths:
            details = []
            if missing_references:
                details.append("missing=" + ", ".join(missing_references[:10]))
            if mismatched_lengths:
                details.append("length_mismatch=" + ", ".join(mismatched_lengths[:10]))
            raise ValueError(
                "BAM references do not match full header: " + "; ".join(details)
            )

    output_bam_path = os.path.abspath(output_bam_path)
    os.makedirs(os.path.dirname(output_bam_path), exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="wt", encoding="utf-8", suffix=".sam", delete=False
    ) as header_file:
        header_file.write(full_header_text)
        header_path = header_file.name
    try:
        subprocess.run(
            ["samtools", "merge", "-f", "-h", header_path, output_bam_path]
            + list(source_bam_paths),
            check=True,
        )
    finally:
        if os.path.exists(header_path):
            os.remove(header_path)
    return output_bam_path

# 阈值可按需调整
COV_MIN = 0.80      # 覆盖度阈值
PID_MIN = 0.90      # 近似一致性阈值
SOFT_MAX = 0.10     # 软截断比例上限


def transcript_to_gene_id(transcript_id):
    """Remove the final transcript suffix used by PanTrans FASTA records."""
    parts = transcript_id.split(".")
    return ".".join(parts[:-1]) if len(parts) > 1 else transcript_id


def iter_filtered_bam_gene_edges(filtered_bam_path):
    """Stream directed gene edges from a BAM already filtered by PanTrans."""
    with pysam.AlignmentFile(filtered_bam_path, "rb") as bam_file:
        for read in bam_file.fetch(until_eof=True):
            if read.is_unmapped or read.is_supplementary:
                continue
            yield (
                transcript_to_gene_id(read.query_name),
                bam_file.get_reference_name(read.reference_id),
            )


def get_bam_target_lengths(bam_path):
    """Return target lengths from a BAM header, keyed by reference name."""
    with pysam.AlignmentFile(bam_path, "rb") as bam_file:
        return dict(zip(bam_file.references, bam_file.lengths))


def alignment_provenance():
    """Describe the executable and thresholds that define graph edges."""
    minimap2_path = shutil.which("minimap2")
    minimap2_version = None
    if minimap2_path:
        result = subprocess.run(
            [minimap2_path, "--version"],
            check=False,
            capture_output=True,
            text=True,
        )
        version_text = (result.stdout or result.stderr).strip()
        minimap2_version = version_text.splitlines()[0] if version_text else None
    return {
        "pantrans_version": __version__,
        "minimap2_path": minimap2_path,
        "minimap2_version": minimap2_version,
        "minimap2_options": MINIMAP2_OPTIONS,
        "filter_thresholds": {
            "coverage_min": COV_MIN,
            "identity_min": PID_MIN,
            "soft_clip_max": SOFT_MAX,
        },
        "filter_logic_id": FILTER_LOGIC_ID,
        "filter_thresholds_assumed": False,
    }


def _minimap2_options_from_command(command):
    tokens = shlex.split(command or "")
    if len(tokens) < 3 or os.path.basename(tokens[0]) != "minimap2":
        return None
    option_tokens = []
    index = 1
    option_end = len(tokens) - 2
    while index < option_end:
        token = tokens[index]
        if token == "-t":
            index += 2
            continue
        if token.startswith("-t") and token[2:].isdigit():
            index += 1
            continue
        option_tokens.append(token)
        index += 1
    return " ".join(option_tokens)


def bam_alignment_provenance(bam_path):
    """Recover edge-generation metadata from an input BAM header when possible."""
    runtime = alignment_provenance()
    with pysam.AlignmentFile(bam_path, "rb") as bam_file:
        header_dict = bam_file.header.to_dict()
        programs = header_dict.get("PG", [])
        comments = header_dict.get("CO", [])
    minimap2_programs = [program for program in programs if program.get("PN") == "minimap2"]
    provenance = dict(runtime)
    provenance.update(
        {
            "scope": "provided_bam",
            "graph_source": "provided_bam",
            "pantrans_version": None,
            "exported_by_pantrans_version": runtime["pantrans_version"],
            "source_bam_programs": programs,
            "filter_thresholds_assumed": True,
            "filter_logic_id": None,
            "pantrans_filtered_bam": FILTERED_BAM_MARKER in comments,
        }
    )
    if minimap2_programs:
        minimap2_program = minimap2_programs[-1]
        provenance["minimap2_version"] = minimap2_program.get("VN")
        provenance["minimap2_options"] = _minimap2_options_from_command(
            minimap2_program.get("CL")
        )
        provenance["minimap2_command"] = minimap2_program.get("CL")
    else:
        provenance["minimap2_version"] = None
        provenance["minimap2_options"] = None
        provenance["minimap2_command"] = None
    return provenance


def validate_resume_bam(bam_path, expected_query_ids, expected_target_lengths, label):
    """Validate a raw cross-alignment BAM before reusing it."""
    provenance = bam_alignment_provenance(bam_path)
    if provenance["pantrans_filtered_bam"]:
        raise ValueError(f"{label} is already PanTrans-filtered; provide the raw BAM")
    actual_target_lengths = get_bam_target_lengths(bam_path)
    if actual_target_lengths != dict(expected_target_lengths):
        raise ValueError(f"{label} target header does not match expected gDNA")
    expected_query_ids = set(expected_query_ids)
    actual_query_ids = set()
    with pysam.AlignmentFile(bam_path, "rb") as bam_file:
        for read in bam_file.fetch(until_eof=True):
            if read.is_secondary or read.is_supplementary:
                continue
            actual_query_ids.add(read.query_name)
            if read.query_name not in expected_query_ids:
                raise ValueError(
                    f"{label} contains query transcript outside expected cDNA: "
                    f"{read.query_name}"
                )
    missing_query_ids = expected_query_ids - actual_query_ids
    if missing_query_ids:
        missing = ", ".join(sorted(missing_query_ids))
        raise ValueError(f"{label} missing expected query transcripts: {missing}")
    return provenance

def minimap2_map(cdna_path, gdna_path, threads, bam_path):
    # use minimap2 to align cdna to gdna and write an unsorted BAM
    cmd = f"minimap2 {MINIMAP2_OPTIONS} -t {threads} {gdna_path} {cdna_path} | samtools view -@{threads} -hb - > {bam_path}"
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

def filter_bam(bam_path, output_bam, bed_dic, collect_edges=True):
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
    output_header = bam_file.header.to_dict()
    comments = list(output_header.get("CO", []))
    if FILTERED_BAM_MARKER not in comments:
        comments.append(FILTERED_BAM_MARKER)
    output_header["CO"] = comments
    out_bam = pysam.AlignmentFile(
        output_bam, "wb", header=pysam.AlignmentHeader.from_dict(output_header)
    )

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
        query_trans = transcript_to_gene_id(read.query_name)
        if query_trans == refer_name:
            out_bam.write(read)
            if collect_edges:
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
            if collect_edges:
                pair_gene_li.append((query_trans, refer_name))

    bam_file.close()
    out_bam.close()

    return pair_gene_li, target_lengths
