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

The CDS workflow remains an experimental script, **not the default implementation
of `pantrans construct` or `pantrans append`**. Existing graph packages cannot be
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
