# Construct/append CLI migration in 0.8.0

The official installed commands now use the CDS workflow:

| Previous entry | 0.8 command |
|---|---|
| `python scripts/experimental_cds_construct.py --prepared P --output O` | `pantrans construct --prepared P --output O` |
| `python scripts/experimental_cds_append.py ...` | `pantrans append ...` |
| `pantrans construct --cdna C --gdna G --bed B ...` | `pantrans legacy-construct --cdna C --gdna G --bed B ...` |
| `pantrans append --history-graph ... --cdna ...` | `pantrans legacy-append --history-graph ... --cdna ...` |

This is an explicit CLI contract change. Legacy arguments are not accepted or
silently routed by the new default commands. Legacy function implementations
remain in `pipeline.py`; the new implementations are `pantrans.cds_construct`
and `pantrans.cds_append`. Experimental script paths are forwarding wrappers,
not separate algorithm implementations. Running immutable historical source
snapshots continues to use their original code.

`preprocess`, `cnv` and `translocation` retain their interfaces. Legacy experiment
runners/validation scripts now call the explicit legacy subcommands. Full new
constructs publish `cds.graph.json`; append checks schema, hashes and policy and
can chain through its `merged_prepared/` directory.

The CDS v2 graph uses 90% identity and both CDS-overlap coverages, without coding
splice/indel vetoes. Only last is annotated by default. UTR-aware structural
deduplication, persistent readers, bounded scheduling and transfer-evidence caches
are shared between new construct and append. Annotation transfer itself retains
strict coding-structure and ORF checks; the relaxed graph filter is not reused as
an annotation acceptance rule.

The current reference/chromosome policy remains wheat-specific (CS, Chr1A–Chr7D).
It was not generalized as part of this CLI migration. The independent ABD-only
partition experiment is not enabled by these commands.

Compatibility script options remain supported, including optional schema-1
`--annotations`, `--genes`, `--annotate-pre`, and annotation-only `--reuse-graph`.
Schema-2 full runs are the supported source of graph packages for append.
Packaged graph restarts validate current graph/input policy; unpackaged restarts
record unverified graph thresholds rather than attributing the current defaults
to an old graph. An append graph consists of multiple generations' alignments,
so the older construct raw-BAM restart interface is not a general append resume.

Validation includes installed-wheel console entrypoints outside the repository,
real minimap2 construct/append runs, repeat append, source/installed output
comparison and legacy parser regressions. Benchmark datasets and generated
experiment output are not shipped in the Python package.


Validation for this migration: 252 tests and 837 subtests passed on Linux with
real minimap2/samtools integration enabled. A built 0.8.0 wheel was installed in
an isolated virtual environment and its console commands were exercised from a
temporary directory with PYTHONPATH removed. Installed preprocess, construct,
append and a direct combined construct produced matching graph/cluster/GFF3/cDNA
and accepted alignment outputs on the smoke dataset. Console --version and
explicit legacy help were checked. This verifies package/CLI integration, not
universal full-cohort equivalence of incremental and direct alignment.
