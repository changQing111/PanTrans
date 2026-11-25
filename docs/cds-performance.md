> Since 0.8.0, use `pantrans construct` and `pantrans append`. Implementations are in the installed `pantrans` package; experimental script paths are compatibility wrappers. See [CLI migration](cli-migration.md). Historical benchmark sections below retain the settings used in those runs.

# CDS performance optimization: evidence cache and measured stage costs

## Implemented behavior

- Construct and append CLIs write `timings.json` with independently measured stages
  (input validation/loading, merge, mapping, filtering, clustering, packaging,
  annotation). Entries also identify a failed stage. Progress counters remain in
  `progress.json`.
- Both CLIs now build `last.evidence.sqlite` by default. This stores each query-to-hub
  transfer before structural deduplication: accepted candidate, accepted SAM record,
  or explicit no-accepted-model status. Original representative seeds are recreated.
- Append automatically reads the historical result's `last.evidence.sqlite` when
  present. `--history-evidence-cache PATH [PATH ...]` supplies additional completed
  caches. `--no-evidence-cache` disables reads and writes for benchmarking.
- Keys bind query ID, sequence and model; target ID, sequence and ordered original
  models; transfer code; minimap2 executable SHA256/version; Biopython and pysam
  versions. Changes produce misses. Query order is reconstructed deterministically.
- Deduplication runs again on the full current cluster, even if all transfers were
  cached. Newly added full models can change historical partial/ambiguous outcomes.
- Parent databases are copied to node-local scratch and opened read-only. Writes
  use a private local database and publish atomically only on successful stage
  exit. The resulting cache includes hits as well as new evidence so it can be
  used without its parents on the next append. `last.cache_stats.json` records
  hits, misses and computed pairs. A stopped/failed stage does not publish a cache.
- Merged FASTAs are copied without constructing SeqRecord objects or rewrapping
  sequence lines. IDs are still checked against metadata; hashes are calculated
  while writing merged files and used immediately in graph-package creation.
  Subsequent runs still rehash historical inputs to detect changes.
- CLI graph inputs `body.fa` and `cds.fa` are relative symlinks to their prepared
  sources, avoiding redundant sequence copies. Keep the prepared source files;
  these outputs are no longer independent standalone FASTA copies. Programmatic
  `load_prepared_bundle` retains its default copy behavior.
- CDS BAM filtering lazily caches unique target CDS structures in gene-body
  coordinates. Equivalent isoforms are checked once; filtering thresholds and
  coding support rules are unchanged.

## Optional node-local reference staging

`--local-fasta` copies cDNA and flanked gDNA plus their indexes to node-local
scratch once per annotation stage. This can reduce shared-filesystem random
fetch latency but adds a large sequential copy and needs sufficient local space.
It remains opt-in pending full-scale timing. `PANTRANS_LOCAL_SCRATCH` selects an
existing scratch directory; otherwise Python's default temporary directory is
used. Files are removed on context exit. `last.local_fasta.json` reports copy
bytes/time and original source paths.

## Existing results and cache warm-up

Previously completed five/ten-variety runs do not contain pre-dedup evidence
caches. They therefore cannot yield immediate hits merely by pointing append at
an old GFF3. A normal new construct or append creates a cache for future reuse;
a construct annotation-only rerun can also create one using `--reuse-graph`.
Do not reinterpret a final GFF3 as the missing original evidence.

The measured five-to-ten experiment has 1,648,751 valid query-to-hub transfers.
402,492 (24.41%) have the same historical hub and member relationship; 62,091
clusters have unchanged hub and gene membership. These are reuse opportunities,
not a promised speedup or an attestation of matching sequence/model fingerprints.
The cache performs the actual content checks.

## Validation

Tests cover cold/partial/full cache reuse; changed query model, target sequence
and algorithm keys; negative evidence; failed-stage nonpublication; read-only
missing-parent handling; local-reference staging; streamed merge hashes; symlink
content; chained append and existing biological behavior. Real-data benchmarks
compare GFF3, cDNA, source/status tables, counts, BAM headers and SAM records and
compare filtered alignments and graph edges before/after CDS-index caching.

Small benchmarks include FASTA-index startup and cannot predict whole-run speed.
Cache hit fraction is not equivalent to wall-time improvement. Reusing cached
transfers does not reuse completed GFF3 rows or skip deduplication; whole-cluster
final-output reuse and local graph updates remain separate future changes.

## Measured regression benchmark (2026-11-09)

Final candidate tested on PKU-IAAS, Slurm job 18812890, node tcu145, 8 annotation
workers. Stratified sample: 757 last clusters, 552 historical clusters, 7,061
query/hub pairs. The sample deliberately covers large/small clusters, unchanged
hubs, new hubs and singletons; it is not a random estimate of whole-run wall time.

| Annotation mode | Seconds | Cache hits / pairs |
|---|---:|---:|
| Previous deployed baseline | 33.887 | — |
| New code, cold cache | 32.749 | 0 / 7,061 |
| Reuse five-variety historical evidence | 32.348 | 3,324 / 7,061 |
| Fully warm evidence | 31.506 | 7,061 / 7,061 |
| Cold cache plus local FASTA staging | 74.593 | 0 / 7,061 |

Partial historical reuse was 4.5% faster in this single run; fully warm evidence
was 7.0% faster. These small differences include noise and startup costs and do
not establish a full-dataset speedup. Local staging copied 18.94 GB in 17.37 s
and was slower end-to-end, so it remains disabled by default and is not currently
recommended for this environment based on this sample.

For 10,000 real raw CDS BAM records, filtering took 2.874 s before and 2.608 s
after target-CDS indexing (9.3% reduction in this sample). Both retained exactly
6,559 records; graph edges and accepted SAM records were identical.

For all annotation modes, GFF3, cDNA, source/status tables, summary counts, BAM
headers and accepted SAM records matched the baseline. The full suite passed
222 tests plus 832 subtests. No full five-plus-five rerun was submitted for this
optimization round. Raw reports are under
`/data/changq/Project/pantrans_test/cds_speed_v2_20261109/`.
