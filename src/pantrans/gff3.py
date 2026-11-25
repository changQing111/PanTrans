"""GFF3 I/O and gene-relative transcript models used by PanTrans.

Rows use the nine GFF3 columns, with decoded seqid and a dict in column nine.
Parent values are lists, so escaped commas in IDs never become delimiters.
"""

import os
from pathlib import Path
import re
import tempfile
from urllib.parse import quote, unquote


def parse_attributes(text):
    attributes = {}
    if text == ".":
        return attributes
    for item in text.split(";"):
        if not item:
            continue
        if "=" not in item:
            raise ValueError("Expected GFF3 key=value attributes")
        key, value = item.split("=", 1)
        key = unquote(key)
        if key in attributes:
            raise ValueError(f"Duplicate GFF3 attribute: {key}")
        attributes[key] = [unquote(v) for v in value.split(",")] if key == "Parent" else unquote(value)
    return attributes


def format_attributes(attributes):
    def escape(value):
        return quote(str(value), safe=".:^*$@!+_?-|")
    return ";".join(escape(key) + "=" + (
        ",".join(escape(v) for v in value) if isinstance(value, list) else escape(value)
    ) for key, value in attributes.items()) or "."


def read_gff3(path):
    rows = []
    version = False
    with open(path, encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            line = line.rstrip("\r\n")
            if line == "##gff-version 3":
                version = True
            if line == "##FASTA":
                raise ValueError(f"{path}: embedded FASTA is not supported; supply a separate FASTA")
            if not line or line.startswith("#"):
                continue
            row = line.split("\t")
            if len(row) != 9:
                raise ValueError(f"{path}:{number}: expected nine GFF3 columns")
            if int(row[3]) < 1 or int(row[4]) < int(row[3]):
                raise ValueError(f"{path}:{number}: invalid coordinates")
            if row[6] not in {"+", "-", ".", "?"}:
                raise ValueError(f"{path}:{number}: invalid strand")
            if row[7] not in ({"0", "1", "2"} if row[2] == "CDS" else {"."}):
                raise ValueError(f"{path}:{number}: invalid phase for {row[2]}")
            row[0], row[8] = unquote(row[0]), parse_attributes(row[8])
            rows.append(row)
    if not version:
        raise ValueError(f"{path}: missing ##gff-version 3 header")
    return rows


def models_from_rows(rows, require_unrenamed=False):
    genes, transcripts, ids = {}, {}, {}
    for row in rows:
        if len(row) != 9 or int(row[3]) < 1 or int(row[4]) < int(row[3]):
            raise ValueError("Invalid GFF3 coordinates or column count")
        if row[6] not in {"+", "-", ".", "?"}:
            raise ValueError("Invalid GFF3 strand")
        if row[7] not in ({"0", "1", "2"} if row[2] == "CDS" else {"."}):
            raise ValueError(f"Invalid GFF3 phase for {row[2]}")
        identifier = row[8].get("ID")
        if not identifier:
            raise ValueError(f"{row[2]} is missing ID")
        if identifier in ids:
            old = ids[identifier]
            # A discontinuous CDS is a single feature with multiple segments.
            if row[2] != "CDS" or old[2] != "CDS" or old[0] != row[0] or old[6] != row[6] or old[8] != row[8]:
                raise ValueError(f"Duplicate/inconsistent GFF3 ID: {identifier}")
        ids[identifier] = row
        if row[2] == "gene":
            if row[8].get("Parent"):
                raise ValueError("Gene must not have a Parent")
            genes[identifier] = row
        elif row[2] in {"mRNA", "transcript"}:
            transcripts[identifier] = row
    transcript_gene, context, exons = {}, {}, {}
    for tid, row in transcripts.items():
        parents = row[8].get("Parent", [])
        if len(parents) != 1 or parents[0] not in genes:
            raise ValueError(f"Transcript {tid} has no unique gene Parent")
        gid = parents[0]
        gene = genes[gid]
        if row[0] != gene[0] or (gene[6] in {"+", "-"} and row[6] != gene[6]):
            raise ValueError(f"Transcript {tid} has mixed sequence names or strands")
        if not int(gene[3]) <= int(row[3]) <= int(row[4]) <= int(gene[4]):
            raise ValueError(f"Transcript {tid} lies outside gene")
        if require_unrenamed:
            if row[0] != gid:
                raise ValueError(f"Sequence name {row[0]} does not match gene ID {gid}")
            prefix, sep, suffix = tid.rpartition(".")
            if not sep or prefix != gid or not suffix.isdigit():
                raise ValueError(f"Transcript ID {tid} must use gene ID {gid} plus a numeric suffix")
        transcript_gene[tid], context[tid], exons[tid] = gid, (row[0], row[6]), []
    for row in rows:
        if row[2] in {"gene", "mRNA", "transcript"}:
            continue
        parents = row[8].get("Parent", [])
        if not parents or any(tid not in transcripts for tid in parents):
            raise ValueError(f"{row[2]} has no transcript Parent")
        for tid in parents:
            trans = transcripts[tid]
            if context[tid] != (row[0], row[6]):
                raise ValueError(f"Transcript {tid} has mixed sequence names or strands")
            if not int(trans[3]) <= int(row[3]) <= int(row[4]) <= int(trans[4]):
                raise ValueError(f"{row[2]} lies outside transcript {tid}")
            if row[2] == "exon":
                exons[tid].append((int(row[3]), int(row[4])))
    splice_sites = {}
    for tid, coords in exons.items():
        coords.sort()
        if not coords:
            raise ValueError(f"Transcript {tid} has no exon records")
        if any(b[0] <= a[1] + 1 for a, b in zip(coords, coords[1:])):
            raise ValueError(f"Transcript {tid} has overlapping or adjacent exons")
        splice_sites[tid] = [(a[1] + 1, b[0] - 1) for a, b in zip(coords, coords[1:])]
    for row in rows:
        if row[2] in {"CDS", "five_prime_UTR", "three_prime_UTR", "start_codon", "stop_codon"}:
            for tid in row[8]["Parent"]:
                if not any(a <= int(row[3]) <= int(row[4]) <= b for a, b in exons[tid]):
                    raise ValueError(f"{row[2]} segment lies outside exons of {tid}")
    return dict(transcript_gene=transcript_gene, exon_coords=exons, splice_sites=splice_sites,
                transcript_context=context, transcript_rows=transcripts, gene_rows=genes)


def load_gff3_transcript_models(path, require_unrenamed=True):
    return models_from_rows(read_gff3(path), require_unrenamed)


def load_history_transcript_models(path):
    """Read GFF3, or legacy PanTrans GTF for migration of historical runs only."""
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("##gff-version") or (line.strip() and not line.startswith("#")):
                if line.startswith("##gff-version") or "ID=" in line.split("\t")[-1]:
                    return load_gff3_transcript_models(path)
                break
    from .legacy_gtf import load_gtf_transcript_models
    return load_gtf_transcript_models(path)


def _natural_key(text):
    return tuple((1, int(v)) if v.isdigit() else (0, v) for v in re.split(r"(\d+)", text))


def write_gff3(path, rows):
    """Validate, sort parents before children, and atomically write one version header."""
    models = models_from_rows(rows)
    def key(row):
        attrs = row[8]
        if row[2] == "gene":
            gene, tid, rank = attrs["ID"], "", 0
        elif row[2] in {"mRNA", "transcript"}:
            gene, tid, rank = attrs["Parent"][0], attrs["ID"], 1
        else:
            tid = attrs["Parent"][0]
            gene, rank = models["transcript_gene"][tid], 2
        return (_natural_key(gene), _natural_key(tid), rank, int(row[3]), int(row[4]), row[2], attrs["ID"])
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=Path(path).resolve().parent,
                                         suffix=".gff3.tmp", delete=False) as out:
            temporary = out.name
            out.write("##gff-version 3\n")
            for row in sorted(rows, key=key):
                fields = row[:]
                fields[0] = quote(fields[0], safe=".:^*$@!+_?-|")
                fields[8] = format_attributes(fields[8])
                out.write("\t".join(fields) + "\n")
        os.replace(temporary, path)
    finally:
        if temporary and os.path.exists(temporary):
            os.unlink(temporary)


def make_row(seqid, feature, start, end, strand, identifier, parent=None, source="Pan", phase="."):
    attrs = {"ID": identifier}
    if parent is not None:
        attrs["Parent"] = [parent]
    return [seqid, source, feature, str(start), str(end), ".", strand, str(phase), attrs]


def rename_rows(rows, rename_map):
    models = models_from_rows(rows, require_unrenamed=True)
    mapping = {gid: rename_map[gid] for gid in models["gene_rows"]}
    for tid, gid in models["transcript_gene"].items():
        mapping[tid] = mapping[gid] + "." + tid.rsplit(".", 1)[1]
    for row in rows:
        identifier = row[8]["ID"]
        if identifier not in mapping:
            tid = row[8]["Parent"][0]
            # PanTrans children use transcript-derived IDs. For imported arbitrary
            # child IDs, namespace them under their renamed transcript.
            suffix = identifier[len(tid):] if identifier.startswith(tid + ".") else "." + identifier
            mapping[identifier] = mapping[tid] + suffix
    if len(set(mapping.values())) != len(mapping):
        raise ValueError("Renaming would produce duplicate GFF3 IDs")
    result = []
    for row in rows:
        new = row[:]
        new[8] = dict(row[8], ID=mapping[row[8]["ID"]])
        if "Parent" in row[8]:
            new[8]["Parent"] = [mapping[p] for p in row[8]["Parent"]]
        result.append(new)
    return result


def rename_gff3_ids(source_path, output_path, rename_map):
    write_gff3(output_path, rename_rows(read_gff3(source_path), rename_map))


def sort_gff3_by_gene_id(path):
    write_gff3(path, read_gff3(path))
