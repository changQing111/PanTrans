# Append clustering independent of variety input order

## Cause and fix

`assign_gene_by_chrom` used the incoming variety list for seed-variety tie
breaking, variety assignment order, and contig traversal. Construct derives that
list from BED order; append concatenates historical and query varieties. Identical
variety sets consequently received different grouping priorities.

The function now canonicalizes unique variety names by `(name.casefold(), name)`.
The second field provides a deterministic tie-break for names differing only by
case. An explicit `Refer` alias retains priority when reference prefixes are used.
Caller lists are not mutated. This is a global deterministic policy, not a list
of priorities fitted to specific wheat varieties. It affects construct and append.

## Validation on 2026-11-02

Runs used the completed original-cDNA append
`20261102T112706_append5plus5_2c05bf21` and direct construct
`20261102T100522_construct10_6b3e42ee`. Node attributes were confirmed identical.
Only Last derivation was repeated using stored pre clusters, unchanged graph
edges, and canonical variety order. Existing result files and GFF3s were retained.

| Metric vs direct construct | Original append | Canonical-order append |
|---|---:|---:|
| Last representative Jaccard |82.20%|90.75%|
| Exact member-set Jaccard |61.28%|83.00%|
| Co-cluster gene-pair Jaccard |81.47%|92.93%|
| Last cluster count |304,532|304,069|

The corrected append has 289,809 shared representatives and a union of 319,361.
All 1,577,859 input genes occur exactly once. This experiment uses the real append
graph, without construct edge substitution. It does not include regenerated final
annotations or a new complete-pipeline runtime measurement.

For the A+B diagnostic graph, the same order correction raises Last representative
Jaccard from 84.90% to 94.16776%. Its representative set and complete cluster
membership sets exactly match the earlier 94.17% controlled experiment. This
resolves the apparent disagreement between those experiments.

An initial case-sensitive sort placed AMN before Abo, unlike the construct order;
it yielded 87.63% on actual append. The final case-insensitive canonical policy
matches the direct construct order in this dataset. Its benefit on other datasets
must be measured; exceeding 90% here is not a universal guarantee.

Regression checks exercise all permutations of three varieties across 40 generated
graphs, covering seed ties, main chromosomes, contigs, unique member preservation,
explicit reference aliases and input immutability. Related suites: 48 tests and
804 subtests passed. Pre-change live source is backed up remotely.

Remote evidence:
`/data/changq/Project/pantrans_test/variety_order_20261102/`.
Final result folders: `append_canonical/`, `AB_canonical/`.
Jobs: 18700936 and 18700937, respectively.
