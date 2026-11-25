"""Predict transcript ORFs and project coding features onto gene-relative GFF3."""

import logging
import os
import shlex
import sys
from time import perf_counter
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
from urllib.parse import unquote

from Bio import SeqIO

from .gff3 import read_gff3, write_gff3, models_from_rows, make_row

logger = logging.getLogger(__name__)
CODING_FEATURES = {"CDS", "five_prime_UTR", "three_prime_UTR", "start_codon", "stop_codon"}


def resolve_transdecoder(directory=None):
    """Return classic two-step executables or the newer wrapper, without running them."""
    def find(name):
        return shutil.which(str(Path(directory).resolve() / name) if directory else name)

    longorfs, predict = find("TransDecoder.LongOrfs"), find("TransDecoder.Predict")
    if longorfs and predict:
        return [longorfs, predict]
    wrapper = find("TransDecoder")
    if wrapper:
        return [wrapper]
    raise FileNotFoundError(
        "TransDecoder not found: install TransDecoder on PATH or use --transdecoder-dir "
        "with the directory containing TransDecoder.LongOrfs/TransDecoder.Predict or TransDecoder."
    )


def _parallel_predict(executables, work, threads):
    """Adapt classic 5.7.1 selection in both modes; parallelize scoring optionally."""
    if len(executables) != 2:
        return None
    predict = Path(executables[1]).resolve()
    if not predict.is_file():
        return None
    source = predict.read_text()
    scoring = "$UTIL_DIR/score_CDS_likelihood_all_6_frames.pl $cds_file $hexamer_scores_file > $cds_scores_file"
    if ('my $VERSION = "5.7.1";' not in source or source.count(scoring) != 1
            or source.count('$FindBin::RealBin') != 2):
        return None
    scorer = predict.parent / "util" / "score_CDS_likelihood_all_6_frames.pl"
    if not scorer.is_file():
        return None
    selector_path = predict.parent / "util" / "select_best_ORFs_per_transcript.pl"
    if not selector_path.is_file():
        return None
    selector = selector_path.read_text()
    tie = '$b->{length} <=> $a->{length}'
    selection_call = '$UTIL_DIR/select_best_ORFs_per_transcript.pl'
    if selector.count(tie) != 2 or selector.count('$FindBin::Bin/../PerlLib') != 1 or source.count(selection_call) != 1:
        return None
    # The original second (stable) sort prioritizes homology and length,
    # retaining the first sort's score ordering for ties. Make that explicit.
    identity = '$a->{gene_obj}->{Model_feat_name} cmp $b->{gene_obj}->{Model_feat_name}'
    selector = selector.replace(tie, tie + '\n            || ' + identity, 1)
    start = selector.index('if ($SINGLE_BEST_ORF_FLAG)')
    selector = selector[:start] + selector[start:].replace(
        tie, tie + '\n                || $b->{cds_scores}->[0] <=> $a->{cds_scores}->[0]\n                || ' + identity, 1)
    selector = selector.replace('$FindBin::Bin/../PerlLib', '$ENV{PANTRANS_TD_ROOT}/PerlLib')
    private_selector = work / "select_best_ORFs_per_transcript.deterministic.pl"
    private_selector.write_text(selector)
    private_selector.chmod(0o755)
    source = source.replace(selection_call, '$ENV{PANTRANS_SELECT_COMMAND}')
    source = source.replace('$FindBin::RealBin', '$ENV{PANTRANS_TD_ROOT}')
    if threads > 1:
        source = source.replace(scoring, '$ENV{PANTRANS_SCORE_COMMAND} $cds_file $hexamer_scores_file > $cds_scores_file')
    adapter = work / ("TransDecoder.Predict.parallel" if threads > 1 else "TransDecoder.Predict.deterministic")
    adapter.write_text(source)
    adapter.chmod(0o755)
    env = os.environ.copy()
    env["PANTRANS_TD_ROOT"] = str(predict.parent)
    env["PANTRANS_SELECT_COMMAND"] = shlex.quote(str(private_selector))
    env["PANTRANS_SCORE_COMMAND"] = shlex.join([
        sys.executable, str(Path(__file__).with_name("transdecoder_parallel.py").resolve()),
        "--workers", str(threads), "--scorer", str(scorer),
    ])
    return str(adapter), env


def run_transdecoder(cdna_path, work_root, directory=None, min_orf_length=50, threads=1):
    """Use a fresh retained run directory, avoiding stale TransDecoder checkpoints."""
    if min_orf_length < 1:
        raise ValueError("Minimum ORF length must be positive")
    if threads < 1:
        raise ValueError("TransDecoder threads must be positive")
    executables = resolve_transdecoder(directory)
    root = Path(work_root).resolve()
    root.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix="run-", dir=str(root)))
    # A local, simple filename also avoids version-dependent output placement.
    fasta = work / "transcripts.fasta"
    shutil.copyfile(cdna_path, fasta)
    common = ["-t", fasta.name]
    if len(executables) == 2:
        commands = [
            [executables[0]] + common + ["-S", "-m", str(min_orf_length)],
            [executables[1]] + common + ["--single_best_only"],
        ]
    else:
        commands = [[executables[0]] + common + ["-S", "-m", str(min_orf_length), "--single_best_only"]]
    parallel = _parallel_predict(executables, work, threads)
    if parallel is None:
        logger.warning("This TransDecoder installation does not support the deterministic adapter; using original serial mode (ORF ties may vary)")
    else:
        logger.info("ORF tie-breaking: original priorities, then lexicographically smallest ORF ID (serial and parallel)")
    if parallel:
        commands[-1][0] = parallel[0]
    log_path = work / "transdecoder.log"
    logger.info("Predict coding regions; log: %s", log_path)
    with log_path.open("w") as log:
        for command in commands:
            log.write("Running: " + repr(command) + "\n")
            log.flush()
            try:
                started = perf_counter()
                options = {"env": parallel[1]} if parallel and command is commands[-1] else {}
                subprocess.run(command, cwd=str(work), stdout=log, stderr=subprocess.STDOUT, check=True, **options)
                elapsed = perf_counter() - started
                log.write("Finished %s in %.3f s\n" % (Path(command[0]).name, elapsed))
                log.flush()
                logger.info("TransDecoder %s finished in %.3f s", Path(command[0]).name, elapsed)
            except (OSError, subprocess.CalledProcessError) as exc:
                raise RuntimeError(f"TransDecoder failed; see {log_path}") from exc
    result = work / "transcripts.fasta.transdecoder.gff3"
    if not result.is_file():
        raise RuntimeError(f"TransDecoder did not produce {result}; see {log_path}")
    return result


def _read_orfs(path):
    models, coding = {}, {}
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError("Malformed TransDecoder GFF3 row")
            attrs = dict(unquote(item).split("=", 1) for item in fields[8].split(";") if "=" in item)
            if fields[2] == "mRNA":
                identifier = attrs["ID"]
                if identifier in models:
                    raise ValueError(f"Duplicate ORF ID: {identifier}")
                match = re.search(r"(?:^|[;\s])(?:ORF\s+)?type:([\w_]+)", attrs.get("Name", ""))
                models[identifier] = (unquote(fields[0]), match.group(1) if match else None)
            elif fields[2] == "CDS":
                if fields[6] != "+":
                    raise ValueError("Expected sense-strand ORFs from TransDecoder -S")
                coding.setdefault(attrs["Parent"], []).append((fields, attrs))
    result = {}
    for identifier, segments in coding.items():
        if identifier not in models or len(segments) != 1:
            raise ValueError(f"Expected one transcript-coordinate CDS for ORF {identifier}")
        tid, kind = models[identifier]
        fields, attrs = segments[0]
        if unquote(fields[0]) != tid or tid in result:
            raise ValueError(f"Inconsistent or multiple ORFs for {tid}; use --single_best_only")
        # Older releases encode completeness in Name; newer ones also use CDS flags.
        if kind not in {"complete", "5prime_partial", "3prime_partial", "internal"}:
            raise ValueError(f"Missing/unsupported TransDecoder ORF type for {identifier}: {kind}")
        complete5 = kind in {"complete", "3prime_partial"} and attrs.get("5_prime_partial") != "true"
        complete3 = kind in {"complete", "5prime_partial"} and attrs.get("3_prime_partial") != "true"
        phase = int(fields[7]) if fields[7] != "." else 0
        if phase not in (0, 1, 2):
            raise ValueError(f"Invalid CDS phase for {identifier}")
        result[tid] = (int(fields[3]), int(fields[4]), phase, complete5, complete3)
    return result


def _project(exons, strand, start, end):
    """Yield genomic intervals in transcription order; all coordinates are 1-based."""
    offset = 1
    for low, high in sorted(exons, reverse=strand == "-"):
        length = high - low + 1
        left, right = max(start, offset), min(end, offset + length - 1)
        if left <= right:
            if strand == "+":
                yield low + left - offset, low + right - offset
            else:
                yield high - right + offset, high - left + offset
        offset += length


def annotate_gff3(gff3_path, cdna_path, prediction_path):
    """Project TransDecoder CDS (including stop codon) onto gene-relative GFF3.

    All segments of a discontinuous CDS share one ID and their mRNA Parent.
    Only CDS has phase; UTR and optional codon features have phase '.'.
    """
    rows = [row for row in read_gff3(gff3_path) if row[2] not in CODING_FEATURES]
    models = models_from_rows(rows, require_unrenamed=True)
    sequences = {}
    for record in SeqIO.parse(str(cdna_path), "fasta"):
        if record.id in sequences:
            raise ValueError(f"Duplicate cDNA ID: {record.id}")
        sequences[record.id] = str(record.seq).upper()
    if set(sequences) != set(models["transcript_gene"]):
        raise ValueError("GFF3 and cDNA transcript IDs do not match")
    orfs = _read_orfs(prediction_path)
    if set(orfs) - set(sequences):
        raise ValueError("TransDecoder contains unknown transcript IDs")
    for tid, exons in models["exon_coords"].items():
        if sum(b - a + 1 for a, b in exons) != len(sequences[tid]):
            raise ValueError(f"GFF3/cDNA length mismatch for {tid}")
    for tid, (start, end, phase, complete5, complete3) in orfs.items():
        seq = sequences[tid]
        seqid, strand = models["transcript_context"][tid]
        if strand not in ("+", "-") or not 1 <= start <= end <= len(seq):
            raise ValueError(f"Invalid ORF coordinates or transcript strand for {tid}")
        if complete5 and (seq[start - 1:start + 2] != "ATG" or phase != 0):
            raise ValueError(f"Invalid complete start codon for {tid}")
        if complete3 and (seq[end - 3:end] not in {"TAA", "TAG", "TGA"} or (end - start + 1 - phase) % 3):
            raise ValueError(f"Invalid complete stop codon/frame for {tid}")
        if end - start + 1 < (6 if complete5 and complete3 else 3):
            raise ValueError(f"ORF is too short for {tid}")
        features = [("CDS", start, end, phase)]
        if complete5:
            features += [("five_prime_UTR", 1, start - 1, None), ("start_codon", start, start + 2, None)]
        if complete3:
            features += [("stop_codon", end - 2, end, None), ("three_prime_UTR", end + 1, len(seq), None)]
        for feature, left, right, frame in features:
            for index, (low, high) in enumerate(_project(models["exon_coords"][tid], strand, left, right), 1):
                identifier = f"{tid}.cds" if feature == "CDS" else f"{tid}.{feature}{index}"
                row = make_row(seqid, feature, low, high, strand, identifier, tid,
                               source="TransDecoder", phase="." if frame is None else frame)
                if feature == "CDS":
                    if not complete5:
                        row[8]["five_prime_partial"] = "true"
                    if not complete3:
                        row[8]["three_prime_partial"] = "true"
                rows.append(row)
                if frame is not None:
                    frame = (frame - (high - low + 1)) % 3
    write_gff3(gff3_path, rows)
    logger.info("Annotated %d/%d transcripts with predicted CDS", len(orfs), len(sequences))
    return len(orfs)
