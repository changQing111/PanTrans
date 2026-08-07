"""Reusable directed-graph metadata for incremental append runs."""

import json
import os
import subprocess
import tempfile


GRAPH_PACKAGE_VERSION = 1
EDGE_PROVENANCE_KEYS = (
    "minimap2_version",
    "minimap2_options",
    "filter_thresholds",
    "filter_logic_id",
)


def _absolute(path):
    return os.path.abspath(os.fspath(path))


def file_identity(path):
    """Return a cheap identity record for an input that must remain stable."""
    path = _absolute(path)
    stat_result = os.stat(path)
    return {
        "path": path,
        "size": stat_result.st_size,
        "mtime_ns": stat_result.st_mtime_ns,
    }


def edge_provenance_generations(provenance):
    """Normalize legacy single-generation metadata into a generation list."""
    provenance = dict(provenance or {})
    generations = provenance.get("edge_generations")
    if generations is None:
        return [provenance] if provenance else []
    if not isinstance(generations, list) or not all(
        isinstance(generation, dict) for generation in generations
    ):
        raise ValueError("Graph edge provenance generations are invalid")
    return [dict(generation) for generation in generations]


def validate_edge_provenance(historical_provenance, current_provenance):
    """Reject historical edge metadata incompatible with the current filter."""
    historical_generations = edge_provenance_generations(historical_provenance)
    if not historical_generations:
        raise ValueError("Graph package has no edge provenance")
    mismatches = []
    for generation_index, generation in enumerate(historical_generations):
        if generation.get("filter_thresholds_assumed", True):
            mismatches.append(
                f"generation {generation_index} has assumed filter thresholds"
            )
        if not generation.get("filter_logic_id"):
            mismatches.append(f"generation {generation_index} has no filter logic ID")
        for key in EDGE_PROVENANCE_KEYS:
            if generation.get(key) != current_provenance.get(key):
                mismatches.append(
                    f"generation {generation_index} {key}: "
                    f"{generation.get(key)!r} != {current_provenance.get(key)!r}"
                )
    if mismatches:
        raise ValueError(
            "incompatible edge provenance; historical graph cannot be safely "
            "combined with current alignments: " + "; ".join(mismatches)
        )
    return historical_generations


def _verify_identity(identity, label):
    path = identity["path"]
    if not os.path.isfile(path):
        raise ValueError(f"{label} is missing: {path}")
    current = file_identity(path)
    if (
        current["size"] != identity.get("size")
        or current["mtime_ns"] != identity.get("mtime_ns")
    ):
        raise ValueError(f"{label} source file changed: {path}")


def _sidecar_identity(manifest_path, sidecar_path):
    manifest_dir = os.path.dirname(_absolute(manifest_path))
    identity = file_identity(sidecar_path)
    identity["path"] = os.path.relpath(identity["path"], manifest_dir)
    return identity


def _resolve_sidecar(manifest_path, identity, label):
    manifest_dir = os.path.dirname(_absolute(manifest_path))
    path = _absolute(os.path.join(manifest_dir, identity["path"]))
    current = file_identity(path)
    if (
        current["size"] != identity.get("size")
        or current["mtime_ns"] != identity.get("mtime_ns")
    ):
        raise ValueError(f"{label} changed: {path}")
    return path


def _read_bed_rows(path):
    rows = []
    rows_by_gene = {}
    with open(path, "rt", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            if not raw_line.strip():
                continue
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) < 4:
                raise ValueError(f"BED row {line_number} in {path} has fewer than 4 columns")
            gene_id = fields[3]
            previous = rows_by_gene.get(gene_id)
            if previous is not None and previous != fields:
                raise ValueError(f"Conflicting BED rows for gene {gene_id} in {path}")
            if previous is None:
                rows_by_gene[gene_id] = fields
                rows.append(fields)
    return rows, rows_by_gene


def merge_history_and_query_bed(
    history_bed_path,
    append_bed_path,
    history_gene_ids,
    query_gene_ids,
    output_path,
):
    """Merge complete historical BED rows with representative-plus-query BED."""
    history_rows, history_by_gene = _read_bed_rows(history_bed_path)
    append_rows, append_by_gene = _read_bed_rows(append_bed_path)
    history_gene_ids = set(history_gene_ids)
    query_gene_ids = set(query_gene_ids)

    allowed_gene_ids = history_gene_ids | query_gene_ids
    unexpected = sorted(set(append_by_gene) - allowed_gene_ids)
    if unexpected:
        preview = ", ".join(unexpected[:10])
        suffix = "..." if len(unexpected) > 10 else ""
        raise ValueError(
            "Append BED contains genes not present in history or query: "
            f"{preview}{suffix}"
        )

    missing_query = sorted(query_gene_ids - set(append_by_gene))
    if missing_query:
        preview = ", ".join(missing_query[:10])
        suffix = "..." if len(missing_query) > 10 else ""
        raise ValueError(f"Append BED is missing query genes: {preview}{suffix}")

    merged_rows = list(history_rows)
    for gene_id, fields in append_by_gene.items():
        if gene_id in history_by_gene:
            if history_by_gene[gene_id] != fields:
                raise ValueError(f"Conflicting historical BED row for gene {gene_id}")
            continue
        merged_rows.append(fields)

    output_path = _absolute(output_path)
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "wt", encoding="utf-8") as handle:
        for fields in merged_rows:
            handle.write("\t".join(fields) + "\n")
    return output_path


def write_graph_package(
    manifest_path,
    edge_iter,
    gene_len_dic,
    bed_path,
    filtered_bam_path,
    cdna_paths,
    gdna_paths,
    variety_names,
    reference_name,
    main_chroms,
    provenance=None,
):
    """Write graph edges, node metadata, and a validated manifest."""
    manifest_path = _absolute(manifest_path)
    os.makedirs(os.path.dirname(manifest_path), exist_ok=True)
    stem, _ = os.path.splitext(manifest_path)
    edges_path = stem + ".edges.tsv"
    nodes_path = stem + ".nodes.tsv"

    bed_rows, bed_by_gene = _read_bed_rows(bed_path)
    missing_lengths = sorted({row[3] for row in bed_rows} - set(gene_len_dic))
    if missing_lengths:
        preview = ", ".join(missing_lengths[:10])
        suffix = "..." if len(missing_lengths) > 10 else ""
        raise ValueError(f"Graph package is missing gene lengths: {preview}{suffix}")
    with open(nodes_path, "wt", encoding="utf-8") as handle:
        for fields in bed_rows:
            chrom, start, end, gene_id = fields[:4]
            strand = fields[5] if len(fields) > 5 else "."
            handle.write(
                f"{gene_id}\t{gene_len_dic[gene_id]}\t{chrom}\t{start}\t{end}\t{strand}\n"
            )

    raw_handle = tempfile.NamedTemporaryFile(
        mode="wt",
        encoding="utf-8",
        prefix=os.path.basename(edges_path) + ".",
        suffix=".unsorted",
        dir=os.path.dirname(edges_path),
        delete=False,
    )
    raw_edges_path = raw_handle.name
    try:
        with raw_handle:
            for query, target in edge_iter:
                if query not in bed_by_gene or target not in bed_by_gene:
                    raise ValueError(
                        "Graph edge endpoint outside graph BED nodes: "
                        f"{query}\t{target}"
                    )
                raw_handle.write(f"{query}\t{target}\n")
        sort_env = os.environ.copy()
        sort_env["LC_ALL"] = "C"
        subprocess.run(
            ["sort", "-u", raw_edges_path, "-o", edges_path],
            check=True,
            env=sort_env,
        )
    finally:
        if os.path.exists(raw_edges_path):
            os.remove(raw_edges_path)

    manifest = {
        "format": "pantrans.graph",
        "version": GRAPH_PACKAGE_VERSION,
        "edges": _sidecar_identity(manifest_path, edges_path),
        "nodes": _sidecar_identity(manifest_path, nodes_path),
        "bed": file_identity(bed_path),
        "filtered_bam": file_identity(filtered_bam_path),
        "cdna": [file_identity(path) for path in cdna_paths],
        "gdna": [file_identity(path) for path in gdna_paths],
        "variety_names": list(variety_names),
        "reference_name": reference_name,
        "main_chroms": list(main_chroms) if main_chroms is not None else None,
        "provenance": dict(provenance or {}),
    }
    with open(manifest_path, "wt", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return manifest_path


def _read_nodes(path):
    gene_len_dic = {}
    with open(path, "rt", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) != 6:
                raise ValueError(f"Graph node row {line_number} in {path} is invalid")
            gene_len_dic[fields[0]] = int(fields[1])
    return gene_len_dic


def load_graph_package(manifest_path):
    """Load and validate a graph package, returning resolved paths and metadata."""
    manifest_path = _absolute(manifest_path)
    with open(manifest_path, "rt", encoding="utf-8") as handle:
        manifest = json.load(handle)
    if manifest.get("format") != "pantrans.graph" or manifest.get("version") != GRAPH_PACKAGE_VERSION:
        raise ValueError(f"Unsupported graph package format: {manifest_path}")

    edges_path = _resolve_sidecar(manifest_path, manifest["edges"], "Graph edge sidecar")
    nodes_path = _resolve_sidecar(manifest_path, manifest["nodes"], "Graph node sidecar")
    for label, identity in (
        ("Graph BED", manifest["bed"]),
        ("Graph filtered BAM", manifest["filtered_bam"]),
    ):
        _verify_identity(identity, label)
    for identity in manifest.get("cdna", []) + manifest.get("gdna", []):
        _verify_identity(identity, "Graph sequence")

    package = dict(manifest)
    package["manifest_path"] = manifest_path
    package["edges_path"] = edges_path
    package["nodes_path"] = nodes_path
    package["bed_path"] = manifest["bed"]["path"]
    package["filtered_bam_path"] = manifest["filtered_bam"]["path"]
    package["cdna_paths"] = [identity["path"] for identity in manifest.get("cdna", [])]
    package["gdna_paths"] = [identity["path"] for identity in manifest.get("gdna", [])]
    package["gene_len_dic"] = _read_nodes(nodes_path)
    package["history_gene_ids"] = set(package["gene_len_dic"])
    return package


def iter_graph_edges(edges_path):
    with open(edges_path, "rt", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) != 2:
                raise ValueError(f"Graph edge row {line_number} in {edges_path} is invalid")
            yield fields[0], fields[1]
