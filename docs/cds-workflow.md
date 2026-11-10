> Since 0.8.0, use `pantrans construct` and `pantrans append`. Implementations are in the installed `pantrans` package; experimental script paths are compatibility wrappers. See [CLI migration](cli-migration.md). Historical benchmark sections below retain the settings used in those runs.

# CDS graph and UTR-aware annotation development branch

This branch contains the current implementation and experimental workflow:

1. `pantrans preprocess` writes complete cDNA, gene-body gDNA, flanked gDNA,
   all annotated CDS, a CDS subset passing complete-ORF/phase checks, coordinate
   tables, and a schema-2 manifest. These are input checks, not graph validation.
2. `scripts/experimental_cds_construct.py` maps CDS against unflanked gene bodies,
   builds the gene graph and pre/last clusters, then maps full cDNA to cluster
   representatives. It preserves original representative annotations and marks
   partial transferred models explicitly.
3. `scripts/restart_cds_annotation.py` reuses completed graph/cluster files and
   reruns annotation with persistent per-worker FASTA handles, a disk cache shared
   across pre/last demands, bounded scheduling and node-local temporary files.
4. `transcript_redundancy.py` folds CDS-only/partial evidence into a unique compatible
   complete structure; ambiguous cases and observed splice differences remain.

The CDS workflow is now the default implementation of `pantrans construct` and
`pantrans append`; the old cDNA workflow uses explicit legacy subcommands. Existing graph packages cannot be
assumed compatible with the experimental CDS graph. Source snapshots and completed
outputs are immutable; code updates do not retroactively change prior experiments.

## Annotation-only restart

```bash
python scripts/restart_cds_annotation.py \
  --source /path/to/PanTrans \
  --old-results /path/to/completed-cds-graph-results \
  --prepared /path/to/preprocessed-inputs \
  --annotations /path/to/original-annotation-list.txt \
  --output /path/to/new-annotation-results \
  --workers 8
```

`--source` may be a repository checkout or a frozen snapshot with `src/` and
`cds_pipeline.py`. Prefer a frozen snapshot for reproducibility. The existing raw
and filtered CDS BAMs, graph edges, pre/last clusters and graph summary must exist.
The output directory must be new. The restart script regenerates annotation
metadata but does not call CDS mapping or clustering. Its manifest records the
reused graph inputs; it does not imply full job-level checkpoint/resume support.

The original annotation list is retained by this restart interface, while schema-2
preprocessing supplies the CDS models directly. No CDS or UTR sequence is invented
to fill clipped ends. `PANTRANS_LOCAL_SCRATCH` may select node-local scratch;
otherwise temporary alignment directories use `/tmp`.

## Validation and limitations

The source/test suite passed 190 tests before publication preparation. The 8-gene
pilot forms the expected 5+3 coding groups and retains a supported UTR splice
variant. A 48-task annotation benchmark produced identical normalized models,
statuses and accepted SAM records, and the restart writer produced identical
pre/last GFF3 output before the subsequent D/E rule change. Those equivalence checks
apply to matching biological rule versions, not across the D/E change.

The 48-task benchmark was about 10x faster; this is not a promised whole-pipeline
speedup. Original representative annotation seeds make missing-model rates
incomparable with older alignment-only, no-rescue pipelines. Partial output models
and failed cross-variety transfers must be evaluated separately. The strict CDS
input policy, graph thresholds and short-ORF behavior still require broader
biological validation. No experiment data or generated output is included here.

## Last-only annotation (2026-11-09)

`experimental_cds_construct.py` now annotates only `last.cluster` by default.
Both cluster tables are still produced during graph construction. Add
`--annotate-pre` only when pre-level annotation is needed for comparison.
The annotation-only restart script uses the same last-only default.

Workers keep separate indexed FASTA readers open for the stage, rather than
reloading the complete indexes per cluster. Scratch files use the node-local
Python temporary directory (normally `/tmp`), or `PANTRANS_LOCAL_SCRATCH` when
set to an existing directory. A bounded rolling queue preserves cluster/output
order and removes the old whole-batch barrier. CDS-block indexing restricts
partial-model compatibility checks to complete models with identical CDS blocks.
Biological filtering and deterministic model selection rules are unchanged.

`progress.json` and stdout report completed/total clusters at stage start,
every 1,000 clusters, and completion. `--workers` controls annotation concurrency;
`--threads` controls the initial CDS alignment. Annotation still invokes minimap2
per nontrivial cluster; it does not perform one unrestricted all-cDNA alignment.

To reuse graph results from the same prepared dataset, choose a NEW output directory:

```bash
python scripts/experimental_cds_construct.py \
  --prepared /path/to/five_variety_prepared \
  --reuse-graph /path/to/five_variety_construct_out \
  --output /path/to/five_variety_last_optimized_out \
  --workers 8
```

Reuse checks that both cluster files partition the prepared gene IDs and records
source paths, sizes and mtimes in `reused_graph.json`. It does not prove that
same-named genes have unchanged sequences: use the original prepared dataset.
Existing graph/BAM/annotation outputs are not overwritten. The new directory
contains copied cluster tables and graph summary; large original graph/BAM files
remain at the paths recorded in `reused_graph.json`.

Validation: 24 real last clusters (126 genes/transcripts), four workers, produced
identical GFF3, cDNA, provenance/status tables, summary, BAM headers and alignment
records. Baseline annotation took 40.09 s versus 8.97 s with persistent readers
and rolling scheduling (4.47x). This small-sample result is not a full-run ETA.

## CDS append (2026-11-09)

The experimental incremental entry is now `scripts/experimental_cds_append.py`.
It reuses historical CDS graph edges, maps new CDS to all gene bodies and old
CDS to new gene bodies, recomputes clusters, and annotates only last clusters.
Fresh full schema-2 constructs emit a separate `cds.graph.json` with input and
edge hashes. See [CDS append](cds-append.md) for input contracts, old-result import,
chaining, validation and all-to-all equivalence limitations.

## Transfer evidence cache and stage timings

See [performance optimization](cds-performance.md) for automatic evidence caching,
per-stage timings, optional node-local reference staging, and the graph FASTA
symlink behavior. Existing historical results have no cache until a new annotation
run creates one; final GFF3 models alone do not reconstruct all original evidence.

## Default homology graph policy v2 (2026-11-10)

The graph now uses 90% PID and both CDS-overlap coverages, without coding splice
or indel structure vetoes. Annotation transfer remains separately validated.
This supersedes earlier graph-threshold descriptions above. See
[CDS homology v2](cds-homology-v2.md) for exact metric definitions and graph-package
compatibility. Existing experimental results retain their original policy.
