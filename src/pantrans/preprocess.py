"""Prepare strand-oriented, coordinate-audited construct/append inputs."""
import csv
import gzip
import json
import logging
import os
from pathlib import Path
import re
import tempfile
from urllib.parse import unquote
from collections import Counter

from Bio.Seq import Seq
import pysam
from .cds_preparation import coding_model

LOG = logging.getLogger(__name__)
SAFE_ID = re.compile(r"^[A-Za-z0-9_][A-Za-z0-9_.:-]*$")
TRANSCRIPT_TYPES = {"mRNA", "transcript", "ncRNA", "lnc_RNA", "lncRNA", "tRNA", "rRNA",
                    "snRNA", "snoRNA", "miRNA", "primary_transcript", "pseudogenic_transcript"}


def read_list(path, paths=False):
    """One entry per line; relative paths are relative to the list file."""
    path = Path(path).resolve()
    values = []
    for line in path.read_text(encoding="utf-8-sig").splitlines():
        value = line.strip()
        if not value or value.startswith("#"):
            continue
        if paths:
            entry = Path(value).expanduser()
            value = str(entry if entry.is_absolute() else path.parent / entry)
        values.append(value)
    return values


def _annotation(path):
    """Read GFF3 hierarchy independent of feature order; never infer CDS as exon."""
    genes, candidates, exons, coding = {}, {}, {}, {}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8-sig") as handle:
        for number, line in enumerate(handle, 1):
            if line.startswith("##FASTA"):
                break
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"{path}:{number}: expected nine GFF3 columns")
            chrom, _, feature, start, end, _, strand, _, attr = fields
            try:
                start, end = int(start) - 1, int(end)
            except ValueError:
                raise ValueError(f"{path}:{number}: invalid coordinates") from None
            if start < 0 or end <= start:
                raise ValueError(f"{path}:{number}: invalid interval")
            attributes = dict(item.split("=", 1) for item in attr.split(";") if "=" in item)
            identifier = unquote(attributes.get("ID", ""))
            parents = [unquote(p) for p in attributes.get("Parent", "").split(",") if p]
            record = (unquote(chrom), start, end, strand)
            if feature in TRANSCRIPT_TYPES and (not identifier or len(parents) != 1):
                raise ValueError(f"{path}:{number}: transcript requires ID and one gene Parent")
            if feature in ("gene", "pseudogene"):
                if not identifier or identifier in genes:
                    raise ValueError(f"{path}:{number}: missing/duplicate gene ID {identifier!r}")
                genes[identifier] = record
            elif feature == "exon":
                if not parents:
                    raise ValueError(f"{path}:{number}: exon lacks Parent")
                for parent in parents:
                    exons.setdefault(parent, set()).add(record)
            elif feature == "CDS":
                if not parents:
                    raise ValueError(f"{path}:{number}: CDS lacks Parent")
                phase = int(fields[7]) if fields[7] in ("0", "1", "2") else -1
                for parent in parents:
                    coding.setdefault(parent, set()).add(record + (phase,))
            elif feature in TRANSCRIPT_TYPES:
                if identifier in candidates:
                    raise ValueError(f"{path}:{number}: duplicate feature ID {identifier}")
                candidates[identifier] = (parents, record)
    transcripts = {}
    for tid, (parents, record) in candidates.items():
        if parents[0] not in genes:
            raise ValueError(f"Transcript {tid} references missing gene {parents[0]}")
        if len(parents) == 1 and parents[0] in genes:
            gene = parents[0]
            g = genes[gene]
            if record[0] != g[0] or record[3] != g[3] or not g[1] <= record[1] < record[2] <= g[2]:
                raise ValueError(f"Transcript {tid} is inconsistent with gene {gene}")
            blocks = sorted(exons.get(tid, ()), key=lambda e: e[1])
            if not blocks:
                raise ValueError(f"Transcript {tid} has no exons; CDS-only annotations are not supported")
            for i, block in enumerate(blocks):
                if block[0] != g[0] or block[3] != g[3] or not record[1] <= block[1] < block[2] <= record[2]:
                    raise ValueError(f"Exon outside transcript or inconsistent strand: {tid}")
                if i and blocks[i - 1][2] > block[1]:
                    raise ValueError(f"Overlapping exons in transcript {tid}")
            transcripts[tid] = (gene, blocks)
    orphan = (set(exons) | set(coding)) - set(transcripts)
    if orphan:
        raise ValueError(f"Exons/CDS have unresolved transcript parents: {sorted(orphan)[:5]}")
    if not transcripts:
        raise ValueError(f"No gene/transcript/exon models found in GFF3: {path}")
    for tid, pieces in coding.items():
        blocks = transcripts[tid][1]
        pieces = sorted(pieces, key=lambda p: p[1])
        for i, (chrom, a, b, strand, phase) in enumerate(pieces):
            if not any(chrom == c and strand == s and x <= a < b <= y for c,x,y,s in blocks):
                raise ValueError(f"CDS outside transcript exons: {tid}")
            if i and pieces[i-1][2] > a:
                raise ValueError(f"Overlapping CDS segments: {tid}")
    return genes, transcripts, coding


def _fasta(handle, name, sequence, description=""):
    handle.write(f">{name}{' ' + description if description else ''}\n")
    for offset in range(0, len(sequence), 80):
        handle.write(sequence[offset:offset + 80] + "\n")


def _local(start, end, extracted_start, extracted_end, strand):
    if strand == "-":
        return extracted_end - end, extracted_end - start
    return start - extracted_start, end - extracted_start


def preprocess(genomes, annotations, varieties, output, prefix="Pan", flank=3000,
               namespace_ids=False):
    """Write combined FASTAs, BED and local-coordinate tables without overwriting."""
    if not genomes or len(genomes) != len(annotations) or len(genomes) != len(varieties):
        raise ValueError("Genome, annotation and variety lists must have equal nonzero lengths")
    if flank < 0:
        raise ValueError("flank must be nonnegative")
    if len(set(varieties)) != len(varieties):
        raise ValueError("Variety names must be unique")
    for name in [prefix] + list(varieties):
        if not SAFE_ID.fullmatch(name) or name in (".", "..", "0"):
            raise ValueError(f"Unsafe or reserved prefix/variety name: {name!r}")
    genomes = [Path(p).expanduser().resolve() for p in genomes]
    annotations = [Path(p).expanduser().resolve() for p in annotations]
    for path in genomes + annotations:
        if not path.is_file():
            raise FileNotFoundError(path)
    # Indexes stay in the temporary workspace, not beside potentially read-only inputs.
    if any(p.suffix.lower() in (".gz", ".bgz", ".bgzf") for p in genomes):
        raise ValueError("Genome FASTA must be uncompressed; GFF3 may be gzip-compressed")
    out = Path(output).resolve()
    out.mkdir(parents=True, exist_ok=True)
    names = {"gdna": f"{prefix}_gene{flank}bp_gdna.fa", "cdna": f"{prefix}_mrna.fa",
             "gene": f"{prefix}_gene.fa", "cds": f"{prefix}_cds.fa",
             "cds_valid": f"{prefix}_cds.valid.fa", "cds_regions": f"{prefix}_cds_regions.tsv",
             "cds_status": f"{prefix}_cds_status.tsv", "models": f"{prefix}_transcript_models.jsonl",
             "bed": f"{prefix}_genes.bed", "genes": f"{prefix}_gene_regions.tsv",
             "exons": f"{prefix}_exon_regions.tsv", "ids": f"{prefix}_id_map.tsv",
             "manifest": f"{prefix}_preprocess.json"}
    fasta_keys = ("gdna", "cdna", "gene", "cds", "cds_valid")
    final_names = list(names.values()) + [names[k] + ".fai" for k in fasta_keys]
    lock = out / f".{prefix}.preprocess.lock"
    with lock.open("x"):
        pass
    try:
        for name in final_names:
            if (out / name).exists():
                raise FileExistsError(f"Output already exists: {out / name}; use a new prefix/output")
        with tempfile.TemporaryDirectory(prefix=f".{prefix}-preprocess-", dir=out) as temp:
            work = Path(temp)
            from contextlib import ExitStack
            with ExitStack() as stack:
                handles = {k: stack.enter_context((work / name).open("w", encoding="utf-8", newline=""))
                           for k, name in names.items() if k != "manifest"}
                gw, ew, iw = [csv.writer(handles[k], delimiter="\t", lineterminator="\n")
                              for k in ("genes", "exons", "ids")]
                gw.writerow(["gene", "variety", "chrom", "gene_start0", "gene_end0", "strand",
                             "extract_start0", "extract_end0", "gdna_length", "body_start0",
                             "body_end0", "actual_upstream", "actual_downstream"])
                ew.writerow(["gene", "transcript", "exon_number", "chrom", "start0", "end0",
                             "gdna_start0", "gdna_end0", "gdna_strand"])
                iw.writerow(["variety", "feature", "original_id", "output_id", "output_gene"])
                cw = csv.writer(handles["cds_regions"], delimiter="\t", lineterminator="\n")
                cw.writerow(["gene", "transcript", "cds_number", "chrom", "start0", "end0",
                             "gdna_start0", "gdna_end0", "gene_start0", "gene_end0", "phase", "gdna_strand"])
                sw = csv.writer(handles["cds_status"], delimiter="\t", lineterminator="\n")
                sw.writerow(["variety", "gene", "transcript", "status", "cds_length"])
                used_genes, used_transcripts = set(), set()
                summaries = []
                for index, (genome, annotation, variety) in enumerate(zip(genomes, annotations, varieties)):
                    LOG.info("Preprocess %s: %s", variety, genome)
                    genes, transcripts, coding = _annotation(annotation)
                    cds_counts = Counter()
                    grouped = {}
                    for tid, (gene, blocks) in transcripts.items():
                        grouped.setdefault(gene, []).append((tid, blocks))
                    fai = work / f"input-{index}.fai"
                    contig_ids = set()
                    with genome.open("rb") as source:
                        for line in source:
                            if line.startswith(b">"):
                                words = line[1:].split()
                                if not words or words[0] in contig_ids:
                                    raise ValueError(f"Empty/duplicate FASTA sequence ID in {genome}")
                                contig_ids.add(words[0])
                    pysam.faidx("-o", str(fai), str(genome))
                    stats = dict(variety=variety, genes=0, transcripts=0, clipped_genes=0,
                                 skipped_genes_without_transcripts=len(genes) - len(grouped))
                    with pysam.FastaFile(str(genome), filepath_index=str(fai)) as reference:
                        ref_lengths = dict(zip(reference.references, reference.lengths))
                        if len(ref_lengths) != len(reference.references):
                            raise ValueError(f"Duplicate FASTA sequence names: {genome}")
                        for gene in sorted(grouped, key=lambda g: (genes[g][0], genes[g][1], g)):
                            chrom, start, end, strand = genes[gene]
                            if chrom not in ref_lengths or end > ref_lengths[chrom] or strand not in ("+", "-"):
                                raise ValueError(f"Gene {gene}: missing contig, out-of-bounds interval or unknown strand")
                            gid = f"{variety}__{gene}" if namespace_ids else gene
                            if not SAFE_ID.fullmatch(gid):
                                raise ValueError(f"Unsupported gene ID {gid!r}; use simple FASTA/GFF-safe IDs")
                            if gid in used_genes:
                                raise ValueError(f"Duplicate gene ID {gid}; use --namespace-ids for distinct varieties")
                            used_genes.add(gid)
                            left, right = max(0, start - flank), min(ref_lengths[chrom], end + flank)
                            sequence = reference.fetch(chrom, left, right).upper()
                            if strand == "-":
                                sequence = str(Seq(sequence).reverse_complement())
                            body_start, body_end = _local(start, end, left, right, strand)
                            upstream, downstream = body_start, len(sequence) - body_end
                            _fasta(handles["gdna"], gid, sequence,
                                   f"chrom={chrom}:{left + 1}-{right} strand={strand} "
                                   f"upstream_downstream={flank} actual_upstream={upstream} actual_downstream={downstream}")
                            _fasta(handles["gene"], gid, sequence[body_start:body_end],
                                   f"chrom={chrom}:{start + 1}-{end} strand={strand} flank=0")
                            handles["bed"].write(f"{chrom}\t{start}\t{end}\t{gid}\t{variety}\t{strand}\n")
                            gw.writerow([gid, variety, chrom, start, end, strand, left, right,
                                         len(sequence), body_start, body_end, upstream, downstream])
                            iw.writerow([variety, "gene", gene, gid, gid])
                            items = sorted(grouped[gene])
                            # Preserve compatible IDs; reserve them before assigning other isoforms.
                            reserved = {tid.rsplit(".", 1)[1] for tid, _ in items
                                        if "." in tid and tid.rsplit(".", 1)[0] == gene
                                        and tid.rsplit(".", 1)[1].isdigit()}
                            counter = 1
                            for tid, blocks in items:
                                if "." in tid and tid.rsplit(".", 1)[0] == gene and tid.rsplit(".", 1)[1].isdigit():
                                    suffix = tid.rsplit(".", 1)[1]
                                else:
                                    while str(counter) in reserved:
                                        counter += 1
                                    suffix = str(counter)
                                    reserved.add(suffix)
                                output_tid = f"{gid}.{suffix}"
                                if output_tid in used_transcripts:
                                    raise ValueError(f"Duplicate output transcript {output_tid}")
                                used_transcripts.add(output_tid)
                                ordered = sorted((_local(a, b, left, right, strand) for _, a, b, _ in blocks))
                                cdna = "".join(sequence[a:b] for a, b in ordered)
                                _fasta(handles["cdna"], output_tid, cdna,
                                       f"gene={gid} chrom={chrom}:{start + 1}-{end} strand={strand}")
                                coding_blocks = sorted((*_local(a,b,left,right,strand), phase)
                                                       for _,a,b,_,phase in coding.get(tid, ()))
                                cds_sequence, status, model = coding_model(sequence, ordered, coding_blocks)
                                model.update(gene=gid, transcript=output_tid, cds_status=status)
                                handles["models"].write(json.dumps(model, sort_keys=True) + "\n")
                                sw.writerow([variety, gid, output_tid, status, len(cds_sequence)])
                                cds_counts[status] += 1
                                if cds_sequence:
                                    _fasta(handles["cds"], output_tid, cds_sequence, f"gene={gid} cds_status={status}")
                                if status == "valid":
                                    _fasta(handles["cds_valid"], output_tid, cds_sequence, f"gene={gid}")
                                for ordinal, (a,b,phase) in enumerate(coding_blocks, 1):
                                    ga,gb = (right-b,right-a) if strand == "-" else (left+a,left+b)
                                    cw.writerow([gid,output_tid,ordinal,chrom,ga,gb,a,b,
                                                 a-body_start,b-body_start,phase,"+"])
                                for ordinal, (a, b) in enumerate(ordered, 1):
                                    genome_a, genome_b = ((right - b, right - a) if strand == "-"
                                                          else (left + a, left + b))
                                    ew.writerow([gid, output_tid, ordinal, chrom, genome_a, genome_b, a, b, "+"])
                                iw.writerow([variety, "transcript", tid, output_tid, gid])
                                stats["transcripts"] += 1
                            stats["genes"] += 1
                            stats["clipped_genes"] += upstream < flank or downstream < flank
                    for label, path in (("genome", genome), ("annotation", annotation)):
                        stat = path.stat()
                        stats[label] = dict(path=str(path), size=stat.st_size, mtime_ns=stat.st_mtime_ns)
                    summaries.append(stats)
                    stats["cds_status_counts"] = dict(cds_counts)
            empty_fastas = []
            for key in fasta_keys:
                if (work / names[key]).stat().st_size:
                    pysam.faidx(str(work / names[key]))
                else:
                    empty_fastas.append(key)
            manifest = dict(schema_version=2, flank=flank, namespace_ids=namespace_ids,
                            empty_fastas=empty_fastas, cds_validation="complete-ORF-phase-v1",
                            graph_query="cds_valid", graph_target="gene",
                            coordinate_system="TSV and BED: 0-based half-open; FASTA headers: 1-based inclusive",
                            gdna_orientation="gene strand; all emitted transcript exons are + relative to gdna",
                            files=names, varieties=summaries)
            (work / names["manifest"]).write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
            # Publish the manifest last; validation failures leave no named outputs.
            for name in final_names:
                if name != names["manifest"] and (work / name).exists():
                    os.replace(work / name, out / name)
            os.replace(work / names["manifest"], out / names["manifest"])
        return {key: str(out / name) for key, name in names.items()}
    finally:
        lock.unlink()
