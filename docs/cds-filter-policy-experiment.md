# CDS filtering experiment: 90% coverage and bounded structural tolerance

This is an opt-in sensitivity experiment, not a change to production defaults.
The driver is `scripts/test_cds_filter_policies.py`; constrained structural
checks live in `src/pantrans/cds_filter_experiment.py`.

The experiment reuses the original five-variety `cds.raw.bam` and schema-2
models. It scans all records, evaluates every relevant valid target isoform,
and tests these scenarios:

| Scenario | PID | Query/target CDS coverage | Structural rule |
|---|---:|---:|---|
| strict95 | .95 | .95 / .95 | Original strict checks |
| strict90 | .95 | .90 / .90 | Original strict checks |
| indel12_cov90 | .95 | .90 / .90 | Local frame restoration allowed |
| junction3_cov90 | .95 | .90 / .90 | Equal donor/acceptor displacement allowed |
| combined_cov90 | .95 | .90 / .90 | Both bounded tolerances |

Indel tolerance: disruptive coding indels can be compensated only when query
and reference spans from the first disruptive event both fit 12 bp, the events
stay within one target CDS block and one query coding interval, and the signed
insertion-minus-deletion length returns to a multiple of three. Unresolved
frameshifts cannot cross an intron, clipping boundary or alignment end. Net
frame restoration does not prove protein identity or biological correctness.

Splice tolerance: each observed coding junction maps to exactly one annotated
target junction; donor and acceptor must move by the same signed amount, at most
3 bp. This preserves intron length. Ambiguous matches, missing target introns,
and novel observed introns remain rejected. This is not an intron-retention
allowance and does not replace canonical-splice or sequence-level validation.

Baseline retained-record count must equal the historical graph summary, and
sampled isoform evaluations are checked against production `coding_support`.
Candidate edges emitted by the sweep may include a pair already supported by
another strict transcript. The graph replay deduplicates these against original
edges before reporting true added edges.

Reclustering is run for all scenarios with the same current SCC→pre code.
Last results are reported separately for (a) the corrected graph-supported
assignment and (b) the user's deterministic ABD / same-variety-same-chromosome
copy separation rule. This avoids conflating a filter change with a last-policy
change. No cDNA transfer/GFF3 annotation is performed. No existing source data,
BAM or published results are overwritten.

Because the raw BAM is fixed, this test cannot recover minimap2 candidates that
were never emitted. A recovered reciprocal edge does not guarantee a singleton
will be absorbed by greedy SCC→pre assignment. Counts of genes moving into a
multi-member pre are measured by replay, rather than assumed from edge counts.

Remote frozen experiment:
`/data/changq/Project/pantrans_test/cds_filter_policy_20261110/`

Validation: 233 tests and 832 subtests passed before the full-data sweep.
Structural tests cover short compensated indels, long/unresolved indels,
intron/block crossings, equal and unequal junction shifts, missing/novel
introns and ambiguous/multiple junction matching.

A real junction-tolerance example recovered by the sweep is
`ZM366Chr1AG197000.1 -> ZM16Chr1AG206300` (target transcript
`ZM16Chr1AG206300.1`): identity 0.999604, both CDS coverages 0.999208.
The observed intron in gene-body 0-based coordinates is `[2230,2321)`;
the target annotation is `[2228,2319)`. Both donor and acceptor shift +2 bp,
with unchanged 91 bp intron length. This supports testing a narrowly bounded
junction tolerance; it does not by itself validate all tolerated edges.

## Full five-variety result

Slurm job 18815356 completed successfully. All 13,440,669 raw records were
scanned; baseline retained count exactly reproduced 10,059,920 and 15,338
sampled isoform decisions matched the original function. Strict95 pre clusters
also exactly matched the original representative/member signatures.

| Scenario | Graph edges | pre clusters | pre singleton | Last (corrected graph policy) | Last (ABD policy) | Original self-only genes now in multi-member pre |
|---|---:|---:|---:|---:|---:|---:|
| strict95 |8382924|196539|73226|247263|243599|0|
| strict90 |8560471|192451|69493|244194|240350|3040|
| indel12_cov90 |8595550|192019|69388|244082|240231|3143|
| junction3_cov90 |8561380|192394|69433|244148|240304|3041|
| combined_cov90 |8596459|191962|69328|244036|240185|3144|

All scenarios had largest SCC=730 and largest pre=117. These two maximum-size
checks do not establish biological correctness or exclude inappropriate small
merges. Coverage relaxation is the dominant contributor in this experiment:
strict90 reduces pre by 4088 and pre singleton by 3733. The combined structural
tolerances add 489 further pre reductions and 165 further singleton reductions,
and recover only 104 additional genes from the original 53,569 self-only set.

The combined scenario leaves 50,425 of those genes as pre singleton, including
49,281 that remain self-loop-only. The difference consists of genes that gained
external edges but still ended up alone after SCC/pre processing. Recovered
gene count is not the number of removed clusters.

Recommended interpretation: 90% bidirectional CDS coverage with 95% PID is a
reasonable next graph-policy candidate to validate on ten varieties. Narrow
structural tolerances recover demonstrably near-matching records but have modest
additional impact. Do not wholesale disable splice/indel checks or silently
reuse relaxed graph rules in strict cDNA annotation transfer. Production defaults
and existing construct/append packages remain unchanged in this experiment.
