# Direct hub support for contig attachment

Contigs may join a cluster only when they have direct edges in both directions
with its representative (the first member). A reciprocal edge to another member
does not qualify. An attached contig cannot recruit additional genes.

This rule applies to the remaining-contig branch of `unit_recursion`, shared by
construct and append, and to final singleton attachment. The remaining-contig
branch keeps the first-compatible-cluster ordering; unmatched contigs establish
new singleton hubs. Final singleton attachment retains the existing preference
for same-variety support, then smaller destination clusters, then cluster order,
but support must now come from the destination hub. A singleton that moves loses
its former hub status. Main-chromosome hub selection is unchanged.

Existing result files do not change. Recompute clusters from the stored graph to
evaluate the new rule, then regenerate downstream annotations if those new
clusters are to be used. Representative selection, alignment inputs, and
historical-edge retention still differ between construct and append; this change
does not make those workflows equivalent or establish biological orthology.

## Validation on 2026-10-28

Tests cover contig-only chains, chromosome-hub chains, one-way connections,
singletons that move, deterministic ordering, unique gene membership, and direct
reciprocal hub support. Legacy-equivalence testing is restricted to paths without
the intentionally changed contig attachment behavior.

The stored induced subgraphs for the `CSChr1DG103800` pre-cluster case were
replayed with the old and new code. The old code reproduced each original
cluster exactly. The new code produced:

| Input case | Original members | New clusters | Hub cluster members | Largest cluster |
|---|---:|---:|---:|---:|
| 5+5 append | 1,812 | 799 | 47 | 47 |
| 10-variety construct | 40 | 6 | 35 | 35 |

All input genes were retained exactly once, and every non-hub member had direct
reciprocal support to its hub. These are local induced-subgraph replays, not a
full graph/SCC rerun or a new GFF3 comparison.

Remote validation and pre-change backups are stored under
`/data/changq/Project/pantrans_test/contig_hub_fix_20261028/` on PKU-IAAS.
The focused clustering, append/construct flow, graph append, and gene-renaming
suite passed 49 tests and 804 subtests using the original job's cactus tool PATH.
