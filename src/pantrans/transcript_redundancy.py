"""Conservative structural deduplication with explicit missing-end evidence.

All models must already be on one representative's forward reference coordinates.
Missing sequence is not evidence for a splice event or for joining two isoforms.
"""

def _blocks(model, name):
    return tuple(sorted(tuple(pair) for pair in model[name]))


def _key(model):
    exons, cds = _blocks(model, "exons"), _blocks(model, "cds")
    junctions = tuple((left[1], right[0]) for left, right in zip(exons, exons[1:]))
    return junctions, cds, () if cds or junctions else exons


def _rank(model):
    return (model["complete"], sum(b-a for a, b in model["exons"]),
            model["pid"], model["source"])


def evidence(model):
    """Alignment completeness is independent of missing UTR annotation."""
    exons, cds = _blocks(model, "exons"), _blocks(model, "cds")
    cds_only = bool(cds) and exons == cds
    unknown5 = bool(model.get("missing_5prime", False)) or cds_only
    unknown3 = bool(model.get("missing_3prime", False)) or cds_only
    return dict(cds_only=cds_only, missing_5prime=unknown5, missing_3prime=unknown3,
                kind="cds_only" if cds_only else "partial" if not model["complete"] else "observed")


def compatible(short, full):
    """Require exact CDS and no disagreement anywhere actually observed.

    Full exons clipped to the short model's observed span must equal that model.
    This rejects intron retention, new internal junctions and crossing an intron.
    Any extension beyond that span requires explicit missing-end evidence.
    """
    if not full["complete"] or evidence(full)["cds_only"]:
        return False
    if not short["cds"] or _blocks(short, "cds") != _blocks(full, "cds"):
        return False
    s, f = _blocks(short, "exons"), _blocks(full, "exons")
    start, end = s[0][0], s[-1][1]
    if f[0][0] > start or f[-1][1] < end:
        return False
    ev = evidence(short)
    if f[0][0] < start and not ev["missing_5prime"]:
        return False
    if f[-1][1] > end and not ev["missing_3prime"]:
        return False
    clipped = tuple((max(a, start), min(b, end)) for a, b in f
                    if max(a, start) < min(b, end))
    return clipped == s


def deduplicate_models(candidates):
    """Fold incomplete evidence only into a unique observed structure group.

    Targets are frozen before folding: incomplete candidates never bridge groups.
    Exact full-structure duplicates retain the previous deterministic ranking.
    Input dictionaries are not modified. Source provenance survives every fold.
    """
    anchors, incomplete = {}, []
    for candidate in candidates:
        ev = evidence(candidate)
        if candidate["complete"] and not ev["cds_only"]:
            anchors.setdefault(_key(candidate), []).append(candidate)
        else:
            incomplete.append(candidate)
    entries = {}
    for key, members in anchors.items():
        entries[("full", key)] = dict(members=list(members), selected=max(members, key=_rank),
                                     status="observed", compatible_sources=[], provenance=[])
    for short in incomplete:
        matches = [key for key, members in anchors.items()
                   if any(compatible(short, full) for full in members)]
        targets = sorted({m["source"] for key in matches for m in anchors[key]})
        ev = evidence(short)
        if len(matches) == 1:
            entry = entries[("full", matches[0])]
            entry["members"].append(short)
            entry["provenance"].append(dict(source=short["source"],
                relation="unique_partial_compatible", evidence=ev["kind"]))
            # The observed full model stays selected even if the CDS-only seed
            # was previously labeled complete or had a better alignment score.
        else:
            status = "ambiguous" if matches else "unresolved"
            # Distinct observed ends/evidence are not collapsed merely because
            # a partial chain shares the same internal junctions.
            key = ("incomplete", _key(short), _blocks(short, "exons"),
                   ev["kind"], ev["missing_5prime"], ev["missing_3prime"], tuple(sorted(matches)))
            if key not in entries:
                entries[key] = dict(members=[], selected=short, status=status,
                                    compatible_sources=targets, provenance=[])
            entry = entries[key]
            entry["members"].append(short)
            entry["selected"] = max(entry["members"], key=_rank)
    result = []
    for entry in entries.values():
        model = entry["selected"].copy()
        model["sources"] = sorted({m["source"] for m in entry["members"]})
        model["utr_evidence"] = evidence(model)["kind"]
        model["redundancy_status"] = entry["status"]
        model["compatible_sources"] = entry["compatible_sources"]
        folded = {p["source"]: p for p in entry["provenance"]}
        model["source_relations"] = [folded.get(m["source"], dict(
            source=m["source"], relation="representative" if m["source"] == model["source"]
            else "exact_structure_duplicate", evidence=evidence(m)["kind"]))
            for m in sorted(entry["members"], key=lambda m: m["source"])]
        result.append(model)
    return sorted(result, key=lambda m: (_blocks(m, "exons"), _blocks(m, "cds"), m["source"]))
