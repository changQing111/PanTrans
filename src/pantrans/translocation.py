"""Screen reciprocal pre-cluster homologs for candidate chromosome changes.

Reports describe positions relative to a reference variety, not evolutionary
direction or confirmed breakpoints. Only unambiguous reciprocal pairs anchor
collinear blocks; unplaced genes are excluded unless explicitly requested.
"""

import csv
import logging
import os
import shutil
import tempfile
from collections import Counter, defaultdict
from dataclasses import dataclass

from .cnv import load_cnv_bed


logger = logging.getLogger(__name__)

GENE_FIELDS = (
    "pre_cluster", "reference_variety", "query_variety", "reference_gene", "query_gene",
    "reference_chrom", "reference_start", "reference_end",
    "query_chrom", "query_start", "query_end", "reciprocal",
    "reference_match_count", "query_match_count",
    "reference_has_same_chrom_match", "query_has_same_chrom_match",
    "status", "reason", "block_id",
)
BLOCK_FIELDS = (
    "block_id", "reference_variety", "query_variety",
    "reference_chrom", "reference_start", "reference_end",
    "query_chrom", "query_start", "query_end", "orientation",
    "anchor_count", "pre_cluster_count", "status",
    "reference_genes", "query_genes", "pre_clusters",
)


@dataclass
class Pair:
    cluster: int
    reference: str
    query: str
    status: str
    reason: str
    reference_matches: int
    query_matches: int
    reference_same_chrom: bool
    query_same_chrom: bool
    block_id: str = "."


def load_pre_membership(path, positions):
    """Index genes by physical line number, rejecting overlapping pre clusters."""
    membership = {}
    count = 0
    with open(path, encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            count += 1
            for gene in line.split():
                if gene not in positions:
                    raise ValueError(f"{path}:{number}: gene absent from BED: {gene}")
                if gene in membership:
                    raise ValueError(f"{path}:{number}: duplicate pre-cluster gene: {gene}")
                membership[gene] = number
    return membership, count


def load_main_chroms(path, positions, reference):
    """Use exact BED chromosome labels; no gene-name or subgenome inference."""
    with open(path, encoding="utf-8") as handle:
        labels = [line.strip() for line in handle
                  if line.strip() and not line.lstrip().startswith("#")]
    if not labels or any(any(char.isspace() for char in label) for label in labels):
        raise ValueError("main-chroms must contain one chromosome name per line")
    main = set(labels)
    if not any(pos.variety == reference and pos.chrom in main for pos in positions.values()):
        raise ValueError("No reference BED chromosomes match --main-chroms (names must match exactly)")
    return main


def reciprocal_reference_pairs(edges_path, membership, positions, reference):
    """Stream edges, retaining only reference/query pairs inside the same cluster.

    Each pair uses two direction bits. Repeated edges never increase support.
    No complete graph or cross-product of cluster members is constructed.
    """
    directions = {}
    with open(edges_path, encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 2 or not all(fields):
                raise ValueError(f"{edges_path}:{number}: expected two tab-separated gene IDs")
            source, target = fields
            if source not in positions or target not in positions:
                raise ValueError(f"{edges_path}:{number}: edge endpoint absent from BED")
            if source == target or source not in membership or target not in membership:
                continue
            if membership[source] != membership[target]:
                continue
            source_ref = positions[source].variety == reference
            target_ref = positions[target].variety == reference
            if source_ref == target_ref:
                continue
            pair = (source, target) if source_ref else (target, source)
            directions[pair] = directions.get(pair, 0) | (1 if source_ref else 2)
    return [pair for pair, bits in directions.items() if bits == 3]


def classify_pairs(reciprocal, membership, positions, main_chroms):
    """Exclude normal same-chromosome pairs, retaining their ambiguity evidence."""
    reference_counts = Counter((ref, positions[query].variety) for ref, query in reciprocal)
    query_counts = Counter(query for _, query in reciprocal)
    reference_same, query_same = set(), set()
    for ref, query in reciprocal:
        r, q = positions[ref], positions[query]
        if r.chrom in main_chroms and q.chrom == r.chrom:
            reference_same.add((ref, q.variety))
            query_same.add(query)

    pairs = []
    normal = 0
    for ref, query in sorted(reciprocal, key=lambda pair: (
        membership[pair[0]], positions[pair[1]].variety, pair[0], pair[1]
    )):
        r, q = positions[ref], positions[query]
        r_count = reference_counts[(ref, q.variety)]
        q_count = query_counts[query]
        r_same, q_same = (ref, q.variety) in reference_same, query in query_same
        if r.chrom not in main_chroms or q.chrom not in main_chroms:
            status, reason = "unplaced", "reference_or_query_not_on_main_chromosome"
        elif r.chrom == q.chrom:
            normal += 1
            continue
        else:
            reasons = []
            if r_same or q_same:
                reasons.append("same_chromosome_alternative")
            if r_count > 1 or q_count > 1:
                reasons.append("multiple_reciprocal_matches")
            status = "ambiguous_copy" if reasons else "single_gene_candidate"
            reason = ";".join(reasons) if reasons else "insufficient_collinear_support"
        pairs.append(Pair(membership[ref], ref, query, status, reason,
                          r_count, q_count, r_same, q_same))
    return pairs, normal


def gene_ranks(positions):
    """Rank ALL BED genes, so intervening non-anchor genes count as skipped."""
    by_chrom = defaultdict(list)
    for gene, pos in positions.items():
        by_chrom[(pos.variety, pos.chrom)].append(gene)
    ranks = {}
    for genes in by_chrom.values():
        genes.sort(key=lambda gene: (positions[gene].start, positions[gene].end, gene))
        ranks.update((gene, rank) for rank, gene in enumerate(genes))
    return ranks


def _can_follow(left, right, positions, ranks, direction, max_gap, max_skipped):
    lr, rr = positions[left.reference], positions[right.reference]
    lq, rq = positions[left.query], positions[right.query]
    if rr.start <= lr.start or direction * (rq.start - lq.start) <= 0:
        return False
    if max(0, rr.start - lr.end) > max_gap:
        return False
    query_gap = rq.start - lq.end if direction == 1 else lq.start - rq.end
    if max(0, query_gap) > max_gap:
        return False
    return (
        ranks[right.reference] - ranks[left.reference] - 1 <= max_skipped
        and abs(ranks[right.query] - ranks[left.query]) - 1 <= max_skipped
    )


def collinear_chains(pairs, positions, ranks, max_gap, max_skipped, allow_reverse):
    """Deterministic disjoint chains, using bounded-window dynamic programming.

    Process longest paths first; stop extraction at an already assigned anchor.
    This is a conservative greedy selection, not exhaustive rearrangement inference.
    Equal start coordinates do not provide independent gene-order evidence.
    """
    pairs = sorted(pairs, key=lambda pair: (
        positions[pair.reference].start, positions[pair.reference].end,
        pair.reference, pair.query,
    ))
    proposals = []
    predecessors = {}
    for direction in ((1, -1) if allow_reverse else (1,)):
        lengths, previous = [1] * len(pairs), [-1] * len(pairs)
        for index, right in enumerate(pairs):
            for earlier in range(index - 1, -1, -1):
                left = pairs[earlier]
                if ranks[right.reference] - ranks[left.reference] - 1 > max_skipped:
                    break
                if _can_follow(left, right, positions, ranks, direction, max_gap, max_skipped):
                    if lengths[earlier] + 1 > lengths[index]:
                        lengths[index] = lengths[earlier] + 1
                        previous[index] = earlier
        predecessors[direction] = previous
        proposals.extend((-length, direction != 1, index, direction)
                         for index, length in enumerate(lengths))
    used = set()
    chains = []
    for _, _, endpoint, direction in sorted(proposals):
        if endpoint in used:
            continue
        path = []
        cursor = endpoint
        while cursor != -1 and cursor not in used:
            path.append(cursor)
            cursor = predecessors[direction][cursor]
        used.update(path)
        chains.append(([pairs[index] for index in reversed(path)], direction))
    return chains


def build_blocks(pairs, positions, reference, min_anchors, max_gap, max_skipped, allow_reverse):
    eligible = defaultdict(list)
    for pair in pairs:
        if pair.status == "single_gene_candidate":
            r, q = positions[pair.reference], positions[pair.query]
            eligible[(q.variety, r.chrom, q.chrom)].append(pair)
    ranks = gene_ranks(positions)
    blocks = []
    for (variety, ref_chrom, query_chrom), members in sorted(eligible.items()):
        chains = collinear_chains(members, positions, ranks, max_gap, max_skipped, allow_reverse)
        chains.sort(key=lambda item: (positions[item[0][0].reference].start, item[0][0].reference))
        for chain, direction in chains:
            if len(chain) < 2:
                continue
            clusters = sorted({pair.cluster for pair in chain})
            supported = len(clusters) >= min_anchors
            block_id = f"TR{len(blocks) + 1:06d}"
            for pair in chain:
                pair.block_id = block_id
                if supported:
                    pair.status, pair.reason = "supported_block", "independent_pre_clusters_in_collinear_block"
            refs = [positions[pair.reference] for pair in chain]
            queries = [positions[pair.query] for pair in chain]
            blocks.append({
                "block_id": block_id, "reference_variety": reference, "query_variety": variety,
                "reference_chrom": ref_chrom, "reference_start": min(pos.start for pos in refs),
                "reference_end": max(pos.end for pos in refs), "query_chrom": query_chrom,
                "query_start": min(pos.start for pos in queries), "query_end": max(pos.end for pos in queries),
                "orientation": "+" if direction == 1 else "-",
                "anchor_count": len(chain), "pre_cluster_count": len(clusters),
                "status": "supported_block" if supported else "low_support_block",
                "reference_genes": ";".join(pair.reference for pair in chain),
                "query_genes": ";".join(pair.query for pair in chain),
                "pre_clusters": ";".join(f"pre_{cluster:06d}" for cluster in clusters),
            })
    return blocks


def gene_report_rows(pairs, positions, reference):
    for pair in pairs:
        r, q = positions[pair.reference], positions[pair.query]
        yield {
            "pre_cluster": f"pre_{pair.cluster:06d}", "reference_variety": reference,
            "query_variety": q.variety, "reference_gene": pair.reference, "query_gene": pair.query,
            "reference_chrom": r.chrom, "reference_start": r.start, "reference_end": r.end,
            "query_chrom": q.chrom, "query_start": q.start, "query_end": q.end, "reciprocal": 1,
            "reference_match_count": pair.reference_matches, "query_match_count": pair.query_matches,
            "reference_has_same_chrom_match": int(pair.reference_same_chrom),
            "query_has_same_chrom_match": int(pair.query_same_chrom),
            "status": pair.status, "reason": pair.reason, "block_id": pair.block_id,
        }


def _write_reports(paths, pairs, blocks, positions, reference):
    """Stage both TSVs before publishing; restore prior files on publish failure."""
    directory = os.path.dirname(paths[0])
    os.makedirs(directory, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".pantrans-translocation-", dir=directory) as work:
        staged, backups = [], {}
        tables = ((GENE_FIELDS, gene_report_rows(pairs, positions, reference)), (BLOCK_FIELDS, blocks))
        for index, (fields, rows) in enumerate(tables):
            stage = os.path.join(work, f"report-{index}.tsv")
            with open(stage, "w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)
            staged.append(stage)
        for index, path in enumerate(paths):
            if os.path.exists(path):
                backup = os.path.join(work, f"backup-{index}.tsv")
                shutil.copy2(path, backup)
                backups[path] = backup
        published = []
        try:
            for stage, path in zip(staged, paths):
                os.replace(stage, path)
                published.append(path)
        except OSError:
            for path in published:
                if path in backups:
                    shutil.copy2(backups[path], path)
                else:
                    os.remove(path)
            raise


def detect_translocations(pre_cluster_path, bed_path, edges_path, main_chrom_path,
                          reference, output_prefix, min_anchor_genes=3, max_gap=1_000_000,
                          max_skipped_genes=5, allow_reverse=False, include_unplaced=False):
    if min_anchor_genes < 2 or max_gap < 0 or max_skipped_genes < 0:
        raise ValueError("Require min-anchor-genes >= 2, max-gap >= 0 and max-skipped-genes >= 0")
    prefix = os.path.abspath(output_prefix)
    paths = [prefix + ".translocation.genes.tsv", prefix + ".translocation.blocks.tsv"]
    sources = (pre_cluster_path, bed_path, edges_path, main_chrom_path)
    for path in paths:
        if os.path.islink(path):
            raise ValueError("Report destination must not be a symbolic link")
        for source in sources:
            if os.path.realpath(path) == os.path.realpath(source) or (
                os.path.exists(path) and os.path.samefile(path, source)
            ):
                raise ValueError("Translocation output must not overwrite an input file")
    positions = load_cnv_bed(bed_path)
    if reference not in {pos.variety for pos in positions.values()}:
        raise ValueError(f"Reference variety {reference!r} is absent from BED column 5")
    main_chroms = load_main_chroms(main_chrom_path, positions, reference)
    membership, cluster_count = load_pre_membership(pre_cluster_path, positions)
    reference_clusters = {membership[gene] for gene in membership if positions[gene].variety == reference}
    reciprocal = reciprocal_reference_pairs(edges_path, membership, positions, reference)
    all_reciprocal_count = len(reciprocal)
    if not include_unplaced:
        # Exclude both endpoints before counting alternative matches, so Ctg
        # genes do not affect placed-pair ambiguity or block eligibility.
        reciprocal = [
            (ref, query) for ref, query in reciprocal
            if positions[ref].chrom in main_chroms and positions[query].chrom in main_chroms
        ]
    pairs, normal = classify_pairs(reciprocal, membership, positions, main_chroms)
    blocks = build_blocks(pairs, positions, reference, min_anchor_genes, max_gap,
                          max_skipped_genes, allow_reverse)
    _write_reports(paths, pairs, blocks, positions, reference)
    counts = Counter(pair.status for pair in pairs)
    summary = {
        "parameters": {
            "reference": reference, "min_anchor_genes": min_anchor_genes,
            "max_gap": max_gap, "max_skipped_genes": max_skipped_genes,
            "allow_reverse": allow_reverse,
            "include_unplaced": include_unplaced,
        },
        "pre_clusters": cluster_count,
        "clusters_without_reference": cluster_count - len(reference_clusters),
        "reciprocal_pairs": len(reciprocal), "same_chromosome_pairs": normal,
        "excluded_unplaced_pairs": all_reciprocal_count - len(reciprocal),
        "reported_pairs": len(pairs), "pair_statuses": dict(counts),
        "supported_blocks": sum(block["status"] == "supported_block" for block in blocks),
        "low_support_blocks": sum(block["status"] == "low_support_block" for block in blocks),
        "genes_path": paths[0], "blocks_path": paths[1],
    }
    logger.info("Translocation screening: %s", summary)
    return summary
