#!/usr/bin/env python3
"""Rebuild official GFF3/cDNA IDs from unrenamed outputs without realignment."""

import argparse
import json
import os
from pathlib import Path
import shutil
import tempfile

from pantrans.common_func import get_bed_varieties
from pantrans.pipeline import _build_gene_rename_map
from pantrans.gff3 import (
    read_gff3, load_gff3_transcript_models,
    rename_gff3_ids,
    sort_gff3_by_gene_id,
)


def repair(bed, output_dir, prefix):
    output_dir = Path(output_dir).resolve()
    varieties = get_bed_varieties(bed, required=True)
    with open(output_dir / f"{prefix}_last.tmp.cluster") as handle:
        clusters = {line.split()[0]: [] for line in handle if line.strip()}
    mapping = _build_gene_rename_map(
        clusters, bed, list(dict.fromkeys(varieties.values())), "Pan",
        gene_variety_dic=varieties,
    )
    gff3_name = f"{prefix}.gff3"
    fasta_name = f"{prefix}_cdna.refer.fasta"
    with tempfile.TemporaryDirectory(prefix=".repair-ids-", dir=output_dir) as tmp:
        stage = Path(tmp)
        rename_gff3_ids(output_dir / f"{prefix}_unrenamed.gff3", stage / gff3_name, mapping)
        sort_gff3_by_gene_id(stage / gff3_name)
        models = load_gff3_transcript_models(stage / gff3_name, require_unrenamed=False)
        transcript_ids = set(models["transcript_gene"])
        gene_ids = set(models["transcript_gene"].values())
        for tid, gene in models["transcript_gene"].items():
            seqid, _ = models["transcript_context"][tid]
            if mapping[seqid] != gene or not tid.startswith(gene + "."):
                raise ValueError("GFF3 IDs do not match the BED-based rename map")
        rows = len(read_gff3(stage / gff3_name))
        fasta_ids = set()
        # Preserve the unrenamed cDNA sequence bytes and descriptions exactly.
        with open(output_dir / f"{prefix}_unrenamed_cdna.refer.fasta", "rb") as source, open(stage / fasta_name, "wb") as dest:
            for line in source:
                if line.startswith(b">"):
                    original = line[1:].split(None, 1)[0]
                    gene, sep, suffix = original.decode().rpartition(".")
                    if not sep or not suffix.isdigit():
                        raise ValueError(f"Invalid unrenamed transcript ID: {original!r}")
                    renamed = mapping[gene] + "." + suffix
                    if renamed in fasta_ids:
                        raise ValueError(f"Duplicate transcript: {renamed}")
                    fasta_ids.add(renamed)
                    line = b">" + renamed.encode() + line[1 + len(original):]
                dest.write(line)
        if fasta_ids != transcript_ids:
            raise ValueError("Repaired GFF3 and cDNA transcript sets differ")
        report = {"genes": len(gene_ids), "transcripts": len(transcript_ids), "gff3_rows": rows}
        # Preserve both old outputs before replacing either one.
        backup = Path(tempfile.mkdtemp(prefix="official-id-backup-", dir=output_dir))
        for name in (gff3_name, fasta_name):
            shutil.copy2(output_dir / name, backup / name)
        report["backup"] = str(backup)
        try:
            for name in (gff3_name, fasta_name):
                os.replace(stage / name, output_dir / name)
        except Exception:
            for name in (gff3_name, fasta_name):
                shutil.copy2(backup / name, output_dir / name)
            raise
        with open(backup / "repair-report.json", "w") as handle:
            json.dump(report, handle, indent=2)
        return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bed", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--prefix", required=True)
    args = parser.parse_args()
    print(json.dumps(repair(args.bed, args.output_dir, args.prefix), indent=2))
