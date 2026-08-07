#!/usr/bin/env python3
"""Export a reusable PanTrans graph package from an existing filtered BAM."""

import argparse
import os

from pantrans.align_filter import (
    MINIMAP2_OPTIONS,
    bam_alignment_provenance,
    get_bam_target_lengths,
    iter_filtered_bam_gene_edges,
)
from pantrans.graph_package import write_graph_package


def _parse_names(value):
    if os.path.isfile(value):
        with open(value, "rt", encoding="utf-8") as handle:
            names = [line.strip() for line in handle if line.strip()]
    else:
        names = [
            item
            for item in value.replace(";", ",").replace(",", " ").split()
            if item
        ]
    if not names:
        raise ValueError("No variety names were provided")
    return names


def _read_lines(path):
    if not path:
        return None
    with open(path, "rt", encoding="utf-8") as handle:
        values = [line.strip() for line in handle if line.strip()]
    if not values:
        raise ValueError(f"No values were found in {path}")
    return values


def read_parameters(argv=None):
    parser = argparse.ArgumentParser(
        description="Export a PanTrans graph package from an existing filtered BAM"
    )
    parser.add_argument("--filtered-bam", required=True)
    parser.add_argument("--cdna", required=True, action="append")
    parser.add_argument("--gdna", required=True, action="append")
    parser.add_argument("--bed", required=True)
    parser.add_argument("--name", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--main-chroms")
    parser.add_argument("--coverage-min", required=True, type=float)
    parser.add_argument("--identity-min", required=True, type=float)
    parser.add_argument("--soft-clip-max", required=True, type=float)
    parser.add_argument("--filter-logic-id", required=True)
    parser.add_argument("--output", required=True, help="output .graph.json path")
    return parser.parse_args(argv)


def main(argv=None):
    args = read_parameters(argv)
    provenance = bam_alignment_provenance(args.filtered_bam)
    provenance["scope"] = "imported_historical_all_to_all"
    provenance["graph_source"] = "existing_pantrans_filtered_bam"
    thresholds = {
        "coverage_min": args.coverage_min,
        "identity_min": args.identity_min,
        "soft_clip_max": args.soft_clip_max,
    }
    if any(value < 0.0 or value > 1.0 for value in thresholds.values()):
        raise ValueError("Filter thresholds must be between 0 and 1")
    provenance["filter_thresholds"] = thresholds
    provenance["filter_thresholds_assumed"] = False
    provenance["filter_thresholds_source"] = "user_supplied_for_existing_filtered_bam"
    provenance["filter_logic_id"] = args.filter_logic_id
    if provenance["minimap2_options"] != MINIMAP2_OPTIONS:
        raise ValueError(
            "Historical BAM minimap2 options do not match this PanTrans version: "
            f"{provenance['minimap2_options']!r} != {MINIMAP2_OPTIONS!r}"
        )
    write_graph_package(
        manifest_path=args.output,
        edge_iter=iter_filtered_bam_gene_edges(args.filtered_bam),
        gene_len_dic=get_bam_target_lengths(args.filtered_bam),
        bed_path=args.bed,
        filtered_bam_path=args.filtered_bam,
        cdna_paths=args.cdna,
        gdna_paths=args.gdna,
        variety_names=_parse_names(args.name),
        reference_name=args.reference,
        main_chroms=_read_lines(args.main_chroms),
        provenance={"edge_generations": [provenance]},
    )
    print(os.path.abspath(args.output))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
