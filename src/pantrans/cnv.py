"""Screen pre clusters for nearby same-variety gene copies."""

import logging
import os
import tempfile
from collections import defaultdict
from typing import NamedTuple


logger = logging.getLogger(__name__)


class GenePosition(NamedTuple):
    chrom: str
    start: int
    end: int
    variety: str


def load_cnv_bed(path):
    """Read gene intervals and variety labels without changing BED coordinates."""
    positions = {}
    with open(path, encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 5:
                raise ValueError(f"{path}:{number}: expected at least five BED columns")
            chrom, start_text, end_text, gene, variety = fields[:5]
            if not chrom or not gene or any(c.isspace() for c in chrom + gene):
                raise ValueError(f"{path}:{number}: invalid chromosome or gene ID")
            variety = variety.strip()
            if variety in {"", ".", "0"}:
                raise ValueError(f"{path}:{number}: missing variety in BED column 5")
            try:
                start, end = int(start_text), int(end_text)
            except ValueError as exc:
                raise ValueError(f"{path}:{number}: invalid BED coordinates") from exc
            if start < 0 or end <= start:
                raise ValueError(f"{path}:{number}: require 0 <= start < end")
            position = GenePosition(chrom, start, end, variety)
            if gene in positions and positions[gene] != position:
                raise ValueError(f"{path}:{number}: conflicting BED records for {gene}")
            positions[gene] = position
    return positions


def first_proximity_group(genes, positions, max_gap=1_000_000):
    """Return the first qualifying interval chain, stopping before other groups.

    Groups are visited in their first occurrence order in the pre cluster;
    intervals within a (variety, chromosome) group are sorted by coordinates.
    Running maximum end coordinates handle overlapping/nested gene intervals.
    """
    if max_gap < 0:
        raise ValueError("max_gap must be nonnegative")
    groups = defaultdict(list)
    for gene in genes:
        position = positions[gene]
        groups[(position.variety, position.chrom)].append(gene)
    for members in groups.values():
        if len(members) < 2:
            continue
        members.sort(key=lambda gene: (positions[gene].start, positions[gene].end, gene))
        chain = []
        right = None
        for gene in members:
            position = positions[gene]
            if right is not None and position.start - right > max_gap:
                if len(chain) >= 2:
                    return chain
                chain = []
                right = None
            chain.append(gene)
            right = position.end if right is None else max(right, position.end)
        if len(chain) >= 2:
            return chain
    return []


def select_cnv_representative(genes, positions, reference):
    """Choose across the entire hit cluster: reference, BED length, then ID."""
    return min(
        genes,
        key=lambda gene: (
            positions[gene].variety != reference,
            -(positions[gene].end - positions[gene].start),
            gene,
        ),
    )


def detect_cnv(pre_cluster_path, bed_path, reference, output_path, max_gap=1_000_000):
    """Write one original representative gene ID per hit pre cluster, atomically."""
    if max_gap < 0:
        raise ValueError("max_gap must be nonnegative")
    output_path = os.path.abspath(output_path)
    for source in (pre_cluster_path, bed_path):
        if os.path.realpath(output_path) == os.path.realpath(source) or (
            os.path.exists(output_path) and os.path.samefile(output_path, source)
        ):
            raise ValueError("CNV output must not overwrite an input file")
    positions = load_cnv_bed(bed_path)
    if reference not in {position.variety for position in positions.values()}:
        raise ValueError(f"Reference variety {reference!r} is absent from BED column 5")
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    temporary_path = None
    checked = hits = 0
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", dir=os.path.dirname(output_path),
            prefix=".pantrans-cnv-", delete=False,
        ) as output:
            temporary_path = output.name
            with open(pre_cluster_path, encoding="utf-8") as clusters:
                for number, line in enumerate(clusters, 1):
                    if not line.strip() or line.lstrip().startswith("#"):
                        continue
                    genes = line.split()
                    if len(set(genes)) != len(genes):
                        raise ValueError(f"{pre_cluster_path}:{number}: duplicate gene in cluster")
                    missing = [gene for gene in genes if gene not in positions]
                    if missing:
                        raise ValueError(
                            f"{pre_cluster_path}:{number}: genes absent from BED: "
                            + ", ".join(missing[:10])
                        )
                    checked += 1
                    if first_proximity_group(genes, positions, max_gap):
                        output.write(select_cnv_representative(genes, positions, reference) + "\n")
                        hits += 1
        os.replace(temporary_path, output_path)
        temporary_path = None
    finally:
        if temporary_path is not None:
            os.remove(temporary_path)
    logger.info("CNV screening: %d/%d pre clusters matched; wrote %s", hits, checked, output_path)
    return hits
