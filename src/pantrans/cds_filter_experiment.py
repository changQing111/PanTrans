"""Opt-in structural sensitivity checks; production annotation rules are unchanged."""


def compensated_indels(read, query_blocks, target_blocks, window=12):
    """Allow a locally restored frame only inside one target CDS block.

    Query and reference spans from the first disruptive indel must both fit the
    window. Splice/clip boundaries cannot bridge an unresolved frameshift.
    This tests alignment compatibility, not translated protein equivalence.
    """
    if window < 0:
        raise ValueError('window must be nonnegative')
    q, r = 0, read.reference_start
    opened = None
    delta = 0
    for op,n in read.cigartuples or []:
        q0,r0=q,r
        coding=False
        change=0
        if op in (0,7,8):q+=n;r+=n
        elif op==1:
            coding=any(max(q,a)<min(q+n,b) for a,b in query_blocks)
            q+=n;change=n
        elif op==2:
            coding=any(max(r,a)<min(r+n,b) for a,b in target_blocks)
            r+=n;change=-n
        elif op==3:
            if opened is not None:return False
            r+=n
        elif op in (4,5):
            if opened is not None:return False
            q+=n
        if coding and (opened is not None or n%3):
            if opened is None:
                opened=(q0,r0)
            delta+=change
            if q-opened[0]>window or r-opened[1]>window:return False
            if not any(a<=opened[1] and r<=b for a,b in target_blocks):return False
            if not any(a<=opened[0] and q<=b for a,b in query_blocks):return False
            if delta%3==0:
                opened=None;delta=0
        elif opened is not None and (q-opened[0]>window or r-opened[1]>window):
            return False
    return opened is None


def compatible_junctions(observed, target, start, end, tolerance=3):
    """Match coding introns one-to-one with equal donor/acceptor displacement.

    No intron creation/deletion or retention is allowed. Equal displacement
    preserves intron length. Multiple possible matches are rejected.
    """
    if tolerance < 0:raise ValueError('tolerance must be nonnegative')
    used=set()
    for left,right in sorted(observed):
        candidates=[pair for pair in target if abs(left-pair[0])<=tolerance
                    and left-pair[0]==right-pair[1]]
        if len(candidates)!=1 or candidates[0] in used:return False
        used.add(candidates[0])
    return all(pair in used for pair in target if start<pair[0] and pair[1]<end)
