"""Deterministic graph-supported assignment within fixed per-variety copy slots."""


def maximum_assignment(scores):
    """Rectangular Hungarian matching maximizing lexicographic score sums."""
    n = len(scores)
    if not n:
        return []
    m = len(scores[0])
    if n > m:
        raise ValueError('More genes than available variety slots')
    bound = max(abs(x) for row in scores for cell in row for x in cell)
    base = 2 * n * bound + 1
    weights = []
    for row in scores:
        encoded = []
        for cell in row:
            value = 0
            for x in cell:
                value = value * base + x
            encoded.append(value)
        weights.append(encoded)
    u, v, p, way = [0]*(n+1), [0]*(m+1), [0]*(m+1), [0]*(m+1)
    for i in range(1, n+1):
        p[0] = i
        j0 = 0
        mins, used = [None]*(m+1), [False]*(m+1)
        while True:
            used[j0] = True
            i0, delta, j1 = p[j0], None, 0
            for j in range(1, m+1):
                if used[j]:
                    continue
                cost = -weights[i0-1][j-1] - u[i0] - v[j]
                if mins[j] is None or cost < mins[j]:
                    mins[j], way[j] = cost, j0
                if delta is None or mins[j] < delta:
                    delta, j1 = mins[j], j
            for j in range(m+1):
                if used[j]:
                    u[p[j]] += delta
                    v[j] -= delta
                elif mins[j] is not None:
                    mins[j] -= delta
            j0 = j1
            if p[j0] == 0:
                break
        while j0:
            j1 = way[j0]
            p[j0] = p[j1]
            j0 = j1
    result = [None]*n
    for j in range(1, m+1):
        if p[j]:
            result[p[j]-1] = j-1
    return result


def assign_main_by_graph(groups, by_variety, seed, varieties, bed, graph):
    """Joint matching and strict-improvement swaps; seed anchors stay fixed.

    Objective: complete reciprocal stars, best-hub coverage, reciprocal pairs,
    then chromosome agreement. Eight deterministic rounds bound local repair;
    this is not a claim of globally optimal multi-variety grouping.
    """
    nodes = {gene for genes in by_variety.values() for gene in genes}
    reciprocal = {a: {b for b in graph.successors(a)
                       if b != a and b in nodes and graph.has_edge(b, a)}
                  if a in graph else set() for a in nodes}

    def score(group):
        members = set(group['genes'])
        if not members:
            return (0, 0, 0, 0)
        support = [len(reciprocal[g] & members) for g in members]
        covered = max(support)+1
        return (int(covered == len(members)), covered, sum(support)//2,
                sum(bed[g][0] == group['main_chrom'] for g in members))

    def total(rows):
        values = [score(group) for group in rows]
        return tuple(sum(x[i] for x in values) for i in range(4))

    def place(rows, genes):
        previous = [score(group) for group in rows]
        matrix = []
        for gene in genes:
            matrix.append([tuple(a-b for a,b in zip(
                score(dict(genes=group['genes']+[gene], main_chrom=group['main_chrom'])), old))
                for group,old in zip(rows,previous)])
        for gene,index in zip(genes,maximum_assignment(matrix)):
            rows[index]['genes'].append(gene)

    for variety in varieties:
        if variety != seed:
            place(groups,by_variety[variety])
    for _ in range(8):
        improved = False
        for variety in varieties:
            if variety == seed or not by_variety[variety]:
                continue
            genes = set(by_variety[variety])
            candidate = [dict(genes=[g for g in group['genes'] if g not in genes],
                              main_chrom=group['main_chrom']) for group in groups]
            place(candidate,by_variety[variety])
            if total(candidate) > total(groups):
                groups[:] = candidate
                improved = True
        if not improved:
            break
