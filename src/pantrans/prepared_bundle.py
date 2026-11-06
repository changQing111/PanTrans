"""Consume schema-2 preprocessing outputs without reparsing original GFF3."""
import collections
import csv
from functools import lru_cache
import json
from pathlib import Path
import shutil

from Bio import SeqIO


@lru_cache(maxsize=16)
def manifest(prepared):
    paths = list(Path(prepared).glob("*_preprocess.json"))
    if len(paths) > 1:
        raise ValueError("Multiple preprocessing manifests; use a directory containing one dataset")
    return json.loads(paths[0].read_text()) if paths else {}


def prepared_file(prepared, key):
    names = manifest(str(Path(prepared).resolve())).get("files", {})
    legacy = dict(gdna="five_variety_gene3000bp_gdna.fa", cdna="five_variety_mrna.fa",
                  genes="five_variety_gene_regions.tsv", exons="five_variety_exon_regions.tsv",
                  ids="five_variety_id_map.tsv")
    return Path(prepared) / names.get(key, legacy.get(key, key))


def load_prepared_bundle(prepared, output, selected=None):
    regions = {}
    with prepared_file(prepared, "genes").open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if selected is None or row["gene"] in selected:
                regions[row["gene"]] = row
    if selected is not None and set(selected) != set(regions):
        raise ValueError("Selected genes are absent from preprocessing tables")
    models, gene_tids, counts = {}, collections.defaultdict(list), collections.Counter()
    with prepared_file(prepared, "models").open() as handle:
        for line in handle:
            model = json.loads(line)
            if model["gene"] not in regions:
                continue
            tid = model.pop("transcript")
            if tid in models:
                raise ValueError(f"Duplicate prepared transcript: {tid}")
            for key in ("cds", "exons", "qcds"):
                model[key] = [tuple(block) for block in model[key]]
            models[tid] = model
            gene_tids[model["gene"]].append(tid)
            counts[model["cds_status"]] += 1
    if counts["valid"] == 0:
        raise ValueError("No valid CDS for graph mapping; inspect the preprocessing CDS status report")
    for key, filename in (("gene", "body.fa"), ("cds_valid", "cds.fa")):
        source, destination = prepared_file(prepared, key), Path(output) / filename
        if selected is None:
            shutil.copyfile(source, destination)
        else:
            with source.open() as src, destination.open("w") as dst:
                for record in SeqIO.parse(src, "fasta"):
                    keep = record.id in regions if key == "gene" else record.id in models and models[record.id]["valid"]
                    if keep:
                        SeqIO.write(record, dst, "fasta")
    with (Path(output) / "cds_status.tsv").open("w") as handle:
        handle.write("transcript\tgene\tstatus\tcds_length\n")
        for tid, model in models.items():
            handle.write(f"{tid}\t{model['gene']}\t{model['cds_status']}\t{model['cds_len']}\n")
    (Path(output) / "cds_status.json").write_text(json.dumps(counts, indent=2))
    return regions, models, gene_tids
