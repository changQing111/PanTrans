"""Original-annotation CDS extraction; internal phase is validation, not trimming."""
from Bio.Seq import Seq


def coding_model(sequence, exons, cds):
    exons, cds = sorted(exons), sorted(cds)
    coding = "".join(sequence[a:b] for a, b, phase in cds)
    query_blocks, offset = [], 0
    for a, b in exons:
        for x, y, phase in cds:
            if max(a, x) < min(b, y):
                query_blocks.append((offset + max(a, x)-a, offset + min(b, y)-a))
        offset += b-a
    reason, cumulative = "valid", 0
    for a, b, phase in cds:
        if phase != (3-cumulative % 3) % 3:
            reason = "invalid_phase"
        cumulative += b-a
    if not cds:
        reason = "no_CDS"
    elif any(c not in "ACGT" for c in coding):
        reason = "ambiguous_CDS"
    elif not coding.startswith("ATG") or coding[-3:] not in ("TAA", "TAG", "TGA") or len(coding) % 3:
        reason = "partial_or_invalid_ORF"
    elif "*" in str(Seq(coding).translate())[:-1]:
        reason = "internal_stop"
    elif (sum(b-a for a, b in query_blocks) != len(coding)
          or query_blocks[-1][1]-query_blocks[0][0] != len(coding)):
        reason = "CDS_not_contiguous_in_cDNA"
    return coding, reason, dict(exons=exons, cds=[(a,b) for a,b,_ in cds],
                                qcds=query_blocks, valid=reason == "valid",
                                cds_len=len(coding), cdna_len=offset)
