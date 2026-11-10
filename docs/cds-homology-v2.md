> Since 0.8.0, use `pantrans construct` and `pantrans append`. Implementations are in the installed `pantrans` package; experimental script paths are compatibility wrappers. See [CLI migration](cli-migration.md). Historical benchmark sections below retain the settings used in those runs.

# CDS homology graph policy v2

The default experimental CDS construct/append graph filter is now
`cds-homology-v2`: PID >= 0.90, query CDS coverage >= 0.90 and target CDS coverage
>= 0.90, without coding-intron compatibility or frame-preserving-indel checks.
Unmapped, supplementary, reverse and missing-NM alignments remain rejected.
Queries still come from preprocess `cds.valid.fa`; at least one valid target
isoform must meet the numeric criteria.

Coverage remains the established CDS-overlap metric: aligned query coding bases
projected inside target annotated CDS, divided by full query/target CDS lengths.
Indels still affect NM-based PID and/or coverage. Allowing an indel structurally
does not remove its contribution to these metrics.

Example: ZM366Chr1BG241400.1 -> CSChr1BG266700, 948M, NM=2,
query CDS overlap 906/948, target CDS overlap 906/906, now passes graph filtering.
The target's 42 bp annotated intron no longer vetoes homology evidence.

The graph package policy is versioned and records structure_checks=false. Old
v1 packages are rejected by normal append rather than silently mixing policies.
Recompute/refilter historical raw alignments using the new policy in a new run.
Previous experimental outputs are not automatically recomputed.

This change applies to graph homology only. cDNA annotation transfer still uses
strict coding_support and existing ORF checks; structural deduplication still
preserves observed splice differences. Therefore graph membership alone does
not guarantee all splice variants can be transferred onto a different hub.
Relaxing annotation transfer to represent retained-intron/alternative-CDS models
requires a separate change and validation. No existing GFF3 or BAM is rewritten.
