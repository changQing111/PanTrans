# PanTrans

[中文](README.zh-CN.md)

PanTrans **0.8.0** constructs and extends a CDS homology graph and produces a
UTR-aware pan-transcriptome. `pantrans construct` and `pantrans append` now run
the CDS workflow directly from the installed Python package. The former
experimental script names remain compatibility wrappers.

## Install

```bash
pip install .
pantrans --version
```

Python 3.8+, pysam, NetworkX and Biopython are required. Install minimap2 and
samtools separately and put them on PATH. TransDecoder is not required by the
CDS workflow; it is used only by optional legacy prediction.

## Prepare and construct

```bash
pantrans preprocess \
  --genome-list genomes.txt --annotation-list annotations.txt \
  --variety-list varieties.txt --flank 3000 \
  --prefix sample --output prepared

pantrans construct --prepared prepared --output construct_results \
  --threads 32 --workers 8
```

The parallel lists contain one genome FASTA, GFF3 path and variety name per line
in matching order. Use uncompressed genome FASTAs. Preprocessing emits gene-body
DNA, flanked gene DNA, complete spliced cDNA, valid CDS, coordinate/model tables,
and a schema-2 manifest. Keep original globally unique IDs, or consistently use
`--namespace-ids` for all batches. See [input details](docs/preprocess-cds-inputs.md).

Construct maps **all valid CDS against all unflanked gene bodies**. Graph edges
require identity and both CDS-overlap coverages >=90%; coding indels and splice
structure differences do not veto graph edges. Coverage counts aligned bases
inside a target isoform's annotated CDS, not arbitrary aligned gene-body bases.
Unmapped, reverse, supplementary and missing-NM records remain excluded.

SCCs are split into pre clusters, then the graph-supported copy assignment
produces last clusters. The current policy prefers reference variety **CS** on
ties and recognizes wheat **Chr1A–Chr7D** (plus aliases) as main chromosomes.
This policy is fixed in this release; the standalone ABD repartition experiment
is not the default. See [assignment rules](docs/last-assignment-graph-support.md).

Only **last** is annotated by default (`--annotate-pre` additionally annotates
pre). Full cDNA is transferred to each hub's flanked DNA, combined with original
hub models and deduplicated by CDS/exon/splice evidence. Unique compatible partial
models can fold into full structures; ambiguous/unresolved evidence remains
explicit. **Annotation transfer still checks coding structure and ORF validity
more strictly than graph construction**. Graph homology alone does not guarantee
that every isoform transfers. See [graph policy](docs/cds-homology-v2.md).

## Append prepared varieties

```bash
pantrans append \
  --history-results construct_results --history-prepared prepared \
  --query-prepared new_prepared --output append_results \
  --threads 32 --workers 8
```

Append reuses the complete historical graph, aligns new CDS to all gene bodies
and historical CDS to new gene bodies, then recomputes clusters and last
annotation. It requires disjoint new varieties and gene/transcript IDs and
compatible schema-2 preparation. Original historical inputs are required;
representative-only sequences are insufficient.

For another generation:

```bash
pantrans append \
  --history-results append_results \
  --history-prepared append_results/merged_prepared \
  --query-prepared next_prepared --output next_append_results
```

The CDS graph package checks input/edge hashes and mapping/filter policy. Old
cDNA `.graph.json` packages and CDS v1 packages are not interchangeable with v2.
Unpackaged historical CDS raw-BAM results may be explicitly imported using
`--import-legacy-history`; missing historical coordinate/mapping provenance is
recorded as assumed. See [append contracts](docs/cds-append.md).

Incremental versus direct construction is **not guaranteed identical** because
minimap2's emitted candidates depend on the target set. Compare biological
outputs when validating a new dataset or policy.

## Outputs and execution

| Output | Meaning |
|---|---|
| `pre.cluster`, `last.cluster` | One cluster per line; first gene is the hub |
| `graph.edges.tsv`, `graph.summary.json` | Directed homology graph and counts |
| `cds.graph.json` | CDS graph/input policy and integrity manifest |
| `last.gff3`, `last.cdna.fa` | Final annotations and reconstructed cDNA |
| `last.annotation.bam` | Accepted transfer evidence before deduplication |
| `last.model_sources.tsv`, `last.alignment_status.tsv` | Source/merge evidence and transfer status |
| `last.evidence.sqlite`, `last.cache_stats.json` | Content-addressed pre-dedup transfer cache and counters |
| `last.summary.json`, `summary.json` | Counts and completion summary |
| `progress.json`, `timings.json` | Annotation progress and stage times |

GFF3 coordinates are **relative to the hub's gene-oriented flanked DNA**, not
original chromosome coordinates. Transcript IDs are locally assigned `hub.N`;
use the source table for provenance. Partial models are retained explicitly.
Graph-input `body.fa` and `cds.fa` are links to prepared data; keep those sources.
Choose a new output directory for each run.

`--threads` controls graph mapping; `--workers` controls annotation workers.
Append automatically reuses a historical `last.evidence.sqlite` when available;
additional caches use `--history-evidence-cache`. `--no-evidence-cache` disables
caching. Deduplication is recomputed with all current candidates. `--local-fasta`
is optional local-reference staging and is not enabled by default. See
[performance and validation](docs/cds-performance.md).

`construct --reuse-graph PATH` performs annotation-only restart. It preserves
original graph provenance; unpackaged history has unverified policy and is
reported with unknown graph thresholds. Subset `--genes` runs and annotation-only
restarts do not produce a new full-dataset graph package for append.

## Migration and legacy commands

Old cDNA arguments (`--cdna`, `--gdna`, `--bed`, `--bam`, `--annotate-cds`) belong to
`pantrans legacy-construct` / `pantrans legacy-append`. They are not silently
interpreted as CDS inputs. See [migration](docs/cli-migration.md) and
[legacy reference](LEGACY.md). `preprocess`, `cnv` and `translocation` remain
available. Scientific visualization scripts remain standalone utilities.

```bash
pantrans construct --help
pantrans append --help
pantrans legacy-construct --help
```
