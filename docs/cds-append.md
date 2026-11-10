> Since 0.8.0, use `pantrans construct` and `pantrans append`. Implementations are in the installed `pantrans` package; experimental script paths are compatibility wrappers. See [CLI migration](cli-migration.md). Historical benchmark sections below retain the settings used in those runs.

# CDS append: experimental incremental workflow

`experimental_cds_append.py` extends the CDS gene graph produced by
`experimental_cds_construct.py`. It does not accept the legacy cDNA-based
`pantrans append` graph package. The production CLI now exposes this workflow as `pantrans append`.

## Inputs and execution

```bash
PYTHONPATH=src python scripts/experimental_cds_append.py \
  --history-results /path/to/history_construct \
  --history-prepared /path/to/history_prepared \
  --query-prepared /path/to/new_varieties_prepared \
  --output /path/to/new_append_results \
  --threads 32 --workers 8
```

Both prepared directories must use schema 2 and compatible preprocessing policies.
Variety names, gene IDs and transcript IDs must be disjoint. Use all original
historical genes and transcripts, not only the representative transcriptome.
Output must be a new directory.

1. Validate the history CDS package, original prepared input SHA256 digests,
   graph edge digest, minimap2 version/options, and CDS filtering policy.
2. Merge complete prepared datasets into `merged_prepared/`; check identifiers
   against gene/model tables and create indexes for combined FASTAs.
3. Reuse historical directed edges and retain all gene nodes, including isolates.
4. Align new valid CDS against all gene bodies, including new-to-new mappings.
5. Align historical valid CDS against new gene bodies.
6. Apply the same CDS filtering function as construct; merge directed edges and
   recompute SCCs, representative selection, pre and last clusters globally.
7. Recompute all last-cluster cDNA transfers and UTR-aware deduplication from
   original transcript models using the optimized annotation runtime.

The first version intentionally does not reuse historical GFF3/cDNA annotations.
Representative changes can alter target coordinates, and cluster membership can
change even when a representative is unchanged. Recomputing last annotations
avoids carrying stale models forward. Historical input/output files are read-only.

## Outputs and chained append

- `cds.graph.json`: separate versioned CDS graph manifest, input hashes, mapping
  and filter policy, parent provenance. Emitted after the graph is complete;
  `summary.json` is the marker that annotation also completed successfully.
- `graph.edges.tsv`, `pre.cluster`, `last.cluster`, `graph.summary.json`: merged
  graph and recomputed clusters. Alignment counters in graph.summary describe
  the two new alignment batches; `historical_edges` counts reused edges.
- `new_to_all.raw.bam` / `.filtered.bam`, `history_to_new.raw.bam` / `.filtered.bam`:
  cross-alignment evidence. A batch with no query CDS is skipped and recorded.
- `last.*`, `progress.json`, `summary.json`: final annotation, source/audit tables,
  progress and execution statistics, with the same meanings as CDS construct.
- `merged_prepared/`: all historical and newly added original models/sequences.
- `body.fa`, `cds.fa`, `cds_status.*`: combined graph input sequences and CDS audit.

For another append, use this result directory as `--history-results` and its
`merged_prepared/` as `--history-prepared`. Parent provenance survives the chain.
There is no single `cds.raw.bam` spanning all generations; use the recorded parent
history plus each generation's cross-alignment BAMs.

Fresh schema-2 full constructs now emit `cds.graph.json`. Subset constructs,
schema-1 runs and `--reuse-graph` annotation restarts do not create an attested
full-dataset package. Running snapshots created before this feature are unchanged.

## Existing five-variety results

Old CDS results have no `cds.graph.json`. Use explicit `--import-legacy-history`
with the ORIGINAL graph directory containing `body.fa`, `cds.fa`, `cds.raw.bam`.
For the current experiment that directory is `five_variety_construct_out/`,
not the later `five_variety_last_optimized_out/` annotation restart directory.

Legacy import checks exact graph-query/target FASTA hashes and BAM reference
IDs/lengths, and refilters the historical raw BAM with the current CDS rules.
It writes `history.refiltered.bam`. It cannot prove the old mapping options or
that annotation coordinates are unchanged; this assumption is explicitly
recorded and retained in parent provenance. It does not rewrite the old result
or silently label its provenance as verified. Subsequent appends use the new
package normally.

## Validation and limits

Integration tests invoke actual minimap2 and samtools on miniature prepared
sequence datasets. They compare append versus direct construct edges, pre/last
clusters, GFF3, cDNA, source/status tables, summary counts and BAM records. Cases
cover changed hub selection, a second append generation, legacy import, invalid
new CDS, duplicate varieties/IDs, missing cDNA records, and tampered graph/input
or policy rejection. Existing annotation/filtering tests remain applicable.

Full wheat five-plus-five versus ten-variety equivalence has not yet been
validated. All-to-all candidate competition and `-N 100` truncation can change
which alignments minimap2 emits when the target set changes; an incremental
union is not guaranteed to equal a direct full construct. No old edges are
silently replaced with edges from the full-run comparator. Input hashing and
merging add sequential I/O; this version optimizes mapping reuse and annotation,
not disk capacity. The current wheat reference (`CS`) and main chromosome list
(`Chr1A` through `Chr7D`) follow the experimental construct policy.

## Transfer evidence caching and timing

The construct/append CLIs now record per-stage `timings.json` and persist
pre-deduplication query-to-hub evidence in `last.evidence.sqlite`. Append reuses
matching historical evidence and recomputes deduplication from all current
candidates. Additional caches can be supplied with `--history-evidence-cache`;
`--no-evidence-cache` disables caching. `--local-fasta` optionally stages annotation
references on node-local scratch. Graph-input FASTAs created by these CLIs are
now links to prepared sources. See [performance notes](cds-performance.md) for
invalidation, source retention, tests and limitations.
