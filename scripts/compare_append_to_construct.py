#!/usr/bin/env python3
"""Compare an incremental append graph/clusters with a direct construct run."""

import argparse
import hashlib
import json
import os
import subprocess
import tempfile

from pantrans.align_filter import iter_filtered_bam_gene_edges


def _sort_unique(source_path, output_path):
    environment = os.environ.copy()
    environment["LC_ALL"] = "C"
    subprocess.run(
        ["sort", "-u", source_path, "-o", output_path],
        check=True,
        env=environment,
    )


def _next_line(handle):
    line = handle.readline()
    return line.rstrip("\n") if line else None


def _compare_sorted(left_path, right_path, left_only_path=None, right_only_path=None):
    left_only = right_only = common = 0
    left_output = open(left_only_path, "wt", encoding="utf-8") if left_only_path else None
    right_output = open(right_only_path, "wt", encoding="utf-8") if right_only_path else None
    try:
        with open(left_path, "rt", encoding="utf-8") as left, open(
            right_path, "rt", encoding="utf-8"
        ) as right:
            left_line = _next_line(left)
            right_line = _next_line(right)
            while left_line is not None or right_line is not None:
                if right_line is None or (
                    left_line is not None and left_line < right_line
                ):
                    left_only += 1
                    if left_output:
                        left_output.write(left_line + "\n")
                    left_line = _next_line(left)
                elif left_line is None or right_line < left_line:
                    right_only += 1
                    if right_output:
                        right_output.write(right_line + "\n")
                    right_line = _next_line(right)
                else:
                    common += 1
                    left_line = _next_line(left)
                    right_line = _next_line(right)
    finally:
        if left_output:
            left_output.close()
        if right_output:
            right_output.close()
    return left_only, right_only, common


def _set_metrics(left_path, right_path, left_only_path=None, right_only_path=None):
    left_only, right_only, common = _compare_sorted(
        left_path,
        right_path,
        left_only_path=left_only_path,
        right_only_path=right_only_path,
    )
    left_total = left_only + common
    right_total = right_only + common
    union = left_only + right_only + common
    return {
        "incremental_total": left_total,
        "construct_total": right_total,
        "common": common,
        "incremental_only": left_only,
        "construct_only": right_only,
        "precision": common / left_total if left_total else 0.0,
        "recall": common / right_total if right_total else 0.0,
        "jaccard": common / union if union else 0.0,
    }


def _write_construct_edges(filtered_bam_path, output_path):
    raw_path = output_path + ".unsorted"
    try:
        with open(raw_path, "wt", encoding="utf-8") as output:
            for query, target in iter_filtered_bam_gene_edges(filtered_bam_path):
                output.write(f"{query}\t{target}\n")
        _sort_unique(raw_path, output_path)
    finally:
        if os.path.exists(raw_path):
            os.remove(raw_path)


def _write_reciprocal_pairs(edges_path, output_path, temp_dir):
    reversed_path = os.path.join(temp_dir, os.path.basename(output_path) + ".reversed")
    reciprocal_directed_path = os.path.join(
        temp_dir, os.path.basename(output_path) + ".directed"
    )
    with open(edges_path, "rt", encoding="utf-8") as source, open(
        reversed_path, "wt", encoding="utf-8"
    ) as output:
        for line in source:
            query, target = line.rstrip("\n").split("\t")
            output.write(f"{target}\t{query}\n")
    _sort_unique(reversed_path, reversed_path)

    environment = os.environ.copy()
    environment["LC_ALL"] = "C"
    with open(reciprocal_directed_path, "wt", encoding="utf-8") as output:
        subprocess.run(
            ["comm", "-12", edges_path, reversed_path],
            stdout=output,
            check=True,
            env=environment,
        )
    with open(reciprocal_directed_path, "rt", encoding="utf-8") as source, open(
        output_path, "wt", encoding="utf-8"
    ) as output:
        for line in source:
            query, target = line.rstrip("\n").split("\t")
            if query < target:
                output.write(line)


def _write_cluster_signatures(cluster_path, signature_path):
    raw_path = signature_path + ".unsorted"
    try:
        with open(cluster_path, "rt", encoding="utf-8") as source, open(
            raw_path, "wt", encoding="utf-8"
        ) as output:
            for line_number, line in enumerate(source, start=1):
                members = line.rstrip("\n").split("\t")
                if not members or not members[0]:
                    raise ValueError(
                        f"Empty cluster at line {line_number} in {cluster_path}"
                    )
                representative = members[0]
                digest = hashlib.sha256(
                    "\0".join(sorted(members)).encode("utf-8")
                ).hexdigest()
                output.write(f"{representative}\t{digest}\n")
        _sort_unique(raw_path, signature_path)
    finally:
        if os.path.exists(raw_path):
            os.remove(raw_path)


def _compare_cluster_signatures(
    incremental_path,
    construct_path,
    incremental_only_path,
    construct_only_path,
):
    common_representatives = 0
    exact_members = 0
    incremental_total = 0
    construct_total = 0
    with open(incremental_path, "rt", encoding="utf-8") as incremental, open(
        construct_path, "rt", encoding="utf-8"
    ) as construct, open(
        incremental_only_path, "wt", encoding="utf-8"
    ) as incremental_only, open(
        construct_only_path, "wt", encoding="utf-8"
    ) as construct_only:
        incremental_line = _next_line(incremental)
        construct_line = _next_line(construct)
        while incremental_line is not None or construct_line is not None:
            if construct_line is None:
                incremental_rep = incremental_line.split("\t", 1)[0]
                incremental_total += 1
                incremental_only.write(incremental_rep + "\n")
                incremental_line = _next_line(incremental)
                continue
            if incremental_line is None:
                construct_rep = construct_line.split("\t", 1)[0]
                construct_total += 1
                construct_only.write(construct_rep + "\n")
                construct_line = _next_line(construct)
                continue
            incremental_rep, incremental_digest = incremental_line.split("\t", 1)
            construct_rep, construct_digest = construct_line.split("\t", 1)
            if incremental_rep < construct_rep:
                incremental_total += 1
                incremental_only.write(incremental_rep + "\n")
                incremental_line = _next_line(incremental)
            elif construct_rep < incremental_rep:
                construct_total += 1
                construct_only.write(construct_rep + "\n")
                construct_line = _next_line(construct)
            else:
                incremental_total += 1
                construct_total += 1
                common_representatives += 1
                exact_members += incremental_digest == construct_digest
                incremental_line = _next_line(incremental)
                construct_line = _next_line(construct)

    union = incremental_total + construct_total - common_representatives
    return {
        "incremental_total": incremental_total,
        "construct_total": construct_total,
        "common_representatives": common_representatives,
        "incremental_only_representatives": incremental_total - common_representatives,
        "construct_only_representatives": construct_total - common_representatives,
        "representative_jaccard": common_representatives / union if union else 0.0,
        "exact_members_for_common_representative": exact_members,
        "exact_member_fraction": (
            exact_members / common_representatives if common_representatives else 0.0
        ),
    }


def read_parameters(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--incremental-edges", required=True)
    parser.add_argument("--construct-bam", required=True)
    parser.add_argument("--incremental-pre", required=True)
    parser.add_argument("--construct-pre", required=True)
    parser.add_argument("--incremental-last", required=True)
    parser.add_argument("--construct-last", required=True)
    parser.add_argument("--output-dir", required=True)
    return parser.parse_args(argv)


def main(argv=None):
    args = read_parameters(argv)
    os.makedirs(args.output_dir, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="pantrans-compare-", dir=args.output_dir) as temp_dir:
        construct_edges = os.path.join(temp_dir, "construct.edges.tsv")
        incremental_reciprocal = os.path.join(temp_dir, "incremental.reciprocal.tsv")
        construct_reciprocal = os.path.join(temp_dir, "construct.reciprocal.tsv")
        _write_construct_edges(args.construct_bam, construct_edges)
        _write_reciprocal_pairs(args.incremental_edges, incremental_reciprocal, temp_dir)
        _write_reciprocal_pairs(construct_edges, construct_reciprocal, temp_dir)

        summary = {
            "directed_edges": _set_metrics(args.incremental_edges, construct_edges),
            "reciprocal_pairs": _set_metrics(
                incremental_reciprocal, construct_reciprocal
            ),
            "clusters": {},
        }
        for label, incremental_cluster, construct_cluster in (
            ("pre", args.incremental_pre, args.construct_pre),
            ("last", args.incremental_last, args.construct_last),
        ):
            incremental_signatures = os.path.join(
                temp_dir, f"incremental.{label}.signatures.tsv"
            )
            construct_signatures = os.path.join(
                temp_dir, f"construct.{label}.signatures.tsv"
            )
            _write_cluster_signatures(incremental_cluster, incremental_signatures)
            _write_cluster_signatures(construct_cluster, construct_signatures)
            summary["clusters"][label] = _compare_cluster_signatures(
                incremental_signatures,
                construct_signatures,
                os.path.join(
                    args.output_dir, f"{label}.incremental_only_representatives.txt"
                ),
                os.path.join(
                    args.output_dir, f"{label}.construct_only_representatives.txt"
                ),
            )

    summary_path = os.path.join(args.output_dir, "comparison_summary.json")
    with open(summary_path, "wt", encoding="utf-8") as output:
        json.dump(summary, output, indent=2, sort_keys=True)
        output.write("\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
