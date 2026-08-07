# PanTrans

PanTrans is a command-line toolkit for constructing a pan-genome and pan-transcriptome from multiple varieties in polyploid organisms.

The current implementation builds gene-level relationship graphs from cDNA-to-gDNA alignments, derives `pre` and `last` gene clusters, performs transcript deduplication from filtered BAM alignments, and writes both intermediate and final reference sequence sets.

## Features

- Build pan-gene and pan-transcript clusters from multiple varieties with `construct`
- Reuse an existing BAM file with `construct --bam`
- Control the definition of main chromosomes with `construct --main-chroms`
- Persist a reusable complete directed graph package from `construct`
- Produce both `pre` and `last` outputs:
  `pre.tmp.cluster`, `last.tmp.cluster`, `pre.gtf`, `gtf`, reference FASTA files, and reference BED files
- Append new varieties to an existing reference set with `append`
- Retain both original representative IDs and official `Pan...` IDs in final GTF/cDNA outputs

## Installation

Install in editable mode from the project root:

```bash
git clone https://github.com/changQing111/PanTrans.git
cd PanTrans
pip install -e .
```

If your environment is offline and already has the required build tools installed, you may need:

```bash
pip install -e . --no-build-isolation
```

## Requirements

PanTrans depends on:

- Python 3.8+
- `pysam`
- `networkx`
- `biopython`
- `minimap2`
- `samtools`

## Input files

PanTrans expects the following input types:

- `cdna` FASTA: transcript sequences
- `gdna` FASTA: genomic gene sequences used as alignment targets
- `bed` file: gene coordinates and strands
- variety names: either a text file, a comma-separated string, or a space-separated string

The BED file must include gene IDs in column 4 and chromosome or contig names in column 1.

## Construct

Use `construct` to build a new pan-genome and pan-transcriptome from multiple varieties.

### Command

```bash
pantrans construct \
  --name <variety_file_or_list> \
  --cdna <all_cdna.fasta> \
  --gdna <all_gdna.fasta> \
  --bed <all.bed> \
  --reference <reference_variety_name> \
  --prefix <prefix> \
  --output <output_dir>
```

### Optional arguments

- `--bam`: use an existing cDNA-to-gDNA BAM file and skip the `minimap2` alignment step
- `--main-chroms`: text file containing main chromosome names, one per line
- `--threads`: number of threads for `minimap2`

### Example

```bash
pantrans construct \
  --name test/tmp_variety.txt \
  --cdna test/tmp_cdna.fasta \
  --gdna test/tmp_gdna.fasta \
  --bed test/tmp.bed \
  --reference CS \
  --prefix test_pantrans \
  --threads 32 \
  --output test/pantrans_construct_out
```

### Construct outputs

Given `--prefix test_pantrans`, PanTrans writes:

- `test_pantrans_pre.tmp.cluster`
- `test_pantrans_last.tmp.cluster`
- `test_pantrans_pre.gtf`
- `test_pantrans_pre_cdna.refer.fasta`
- `test_pantrans_pre_gdna.refer.fasta`
- `test_pantrans_pre.refer.bed`
- `test_pantrans_unrenamed.gtf`
- `test_pantrans_unrenamed_cdna.refer.fasta`
- `test_pantrans.gtf`
- `test_pantrans_cdna.refer.fasta`
- `test_pantrans_gdna.refer.fasta`
- `test_pantrans.refer.bed`
- `test_pantrans.graph.json`
- `test_pantrans.graph.edges.tsv`
- `test_pantrans.graph.nodes.tsv`
- alignment BAM files

### Notes on naming

- `pre.gtf` preserves the original cluster hub gene IDs and transcript IDs.
- `<prefix>_unrenamed.gtf` and `<prefix>_unrenamed_cdna.refer.fasta`
  retain the final `last` representative IDs used by the graph, BED, and cluster files.
- `<prefix>.gtf` and `<prefix>_cdna.refer.fasta` contain the same final
  transcript models renamed to official `Pan...` IDs such as `Pan1A000001`.

### Exporting a graph from an existing construct result

Construct runs created before graph packages were introduced do not need to repeat
their historical minimap2 alignment. Export the graph directly from the existing
PanTrans-filtered BAM and the original complete inputs:

```bash
PYTHONPATH=src python scripts/export_graph_package.py \
  --filtered-bam test/pantrans_construct_out/test_pantrans_cdna_align_gdna.filtered.bam \
  --cdna test/tmp_cdna.fasta \
  --gdna test/tmp_gdna.fasta \
  --bed test/tmp.bed \
  --name test/tmp_variety.txt \
  --reference CS \
  --main-chroms test/tmp_chrom.txt \
  --coverage-min 0.80 \
  --identity-min 0.90 \
  --soft-clip-max 0.10 \
  --filter-logic-id pantrans-filter-v1 \
  --output test/pantrans_construct_out/test_pantrans.graph.json
```

The BAM must already be the filtered BAM that generated the historical clusters.
Its reference names must cover every gene in the complete historical BED. The
three filtering thresholds and `--filter-logic-id` must describe the logic used
to create that filtered BAM; BAM headers do not retain PanTrans filtering
settings.
The exporter streams BAM records and records the source BAM program header in
the manifest; it does not rerun alignment or clustering.

## Append

Use `append` to reuse a historical construct graph and add only the new-to-history alignment blocks. The merged graph then follows the same clustering and transcript processing flow as `construct`.

### Command

```bash
pantrans append \
  --name <new_variety_file_or_list> \
  --cdna <previous_unrenamed_final_cdna.fasta> \
  --history-gtf <previous_unrenamed_final.gtf> \
  --query-cdna <new_variety_cdna.fasta> \
  --gdna <query_gdna.fasta> \
  --bed <combined_representatives_and_query.bed> \
  --history-graph <historical_prefix.graph.json> \
  --prefix <prefix> \
  --output <output_dir>
```

The three transcript inputs have distinct roles:

- `--cdna` is the previous run's `<prefix>_unrenamed_cdna.refer.fasta`.
- `--history-gtf` is its matching `<prefix>_unrenamed.gtf`.
- `--query-cdna` is the new variety's complete cDNA FASTA.

The historical GTF and cDNA must contain exactly the same transcript IDs. Their
gene IDs must be unrenamed IDs present in the historical graph package; passing
the official `Pan...` GTF is rejected. `--bed` is one combined BED containing
the first gene ID from every historical `pre.tmp.cluster` plus every gene from
the new variety. The graph package supplies complete historical gDNA, BED,
filtered BAM, graph edges, node lengths, variety order, reference, and alignment
provenance.

For restartability, append also accepts `--query-to-all-bam` and
`--history-to-query-bam` to reuse either completed cross-alignment BAM.

### Example

```bash
pantrans append \
  --name JM22 \
  --cdna test/pantrans_construct_out/test_pantrans_unrenamed_cdna.refer.fasta \
  --history-gtf test/pantrans_construct_out/test_pantrans_unrenamed.gtf \
  --query-cdna test/JM22_cdna.fasta \
  --gdna test/JM22.gdna.fasta \
  --bed test/test_pantrans_pre.refer_append_JM22.bed \
  --history-graph test/pantrans_construct_out/test_pantrans.graph.json \
  --prefix test_append_JM22 \
  --threads 32 \
  --output test/pantrans_append_out
```

### Append workflow

The `append` workflow performs the following steps:

1. Validate the historical graph package and its source-file identities
2. Build a working cDNA from historical nonredundant plus query transcripts,
   and build full historical-plus-query gDNA and BED files
3. Align new cDNA to all historical-plus-new gDNA
4. Align historical nonredundant cDNA to new gDNA
5. Reuse historical graph edges, merge the two cross-edge blocks, and merge their filtered BAMs with the historical filtered BAM
6. Rebuild graph-based `pre` and `last` clusters using exactly the construct clustering flow
7. Seed final transcript deduplication with historical GTF models whose original
   representative gene ID is still a current `last` representative
8. Write official and unrenamed final GTF/cDNA outputs plus a graph package for
   a later append

All sequence-backed historical and new genes are eligible cluster members. The append BED is only the user-facing representative-plus-query BED contract; the historical graph package supplies the complete historical BED and gDNA sequence set used internally. Historical GTF coordinates are reused only when their original representative gene ID remains a current `last` cluster representative; transcript models are not projected onto a different representative coordinate system.

The two cross-alignment blocks avoid recomputing historical cDNA against historical gDNA. Since minimap2 retains a bounded number of secondary alignments, the incremental graph is expected to be highly similar to, but not byte-identical with, a monolithic construct graph.

To append another variety, use the current append run's graph package together
with its `_unrenamed.gtf` and `_unrenamed_cdna.refer.fasta` outputs as the next
run's historical inputs.

If an append job fails after either cross-alignment finishes, rerun it with the
corresponding `--query-to-all-bam` and/or `--history-to-query-bam` path. The BAM
must contain the minimap2 `@PG` record with the same edge-generating options and
version; PanTrans re-filters reused BAMs with the current thresholds.

## Output interpretation

- `pre.tmp.cluster`: graph-derived intermediate clustering before final recursive refinement
- `last.tmp.cluster`: final cluster list after recursive splitting and chromosome-aware assignment
- `pre.gtf`: transcript deduplication result at the `pre` cluster level using original hub gene IDs
- `_unrenamed.gtf`: final `last` transcript models using original representative IDs
- `gtf`: transcript deduplication result at the `last` cluster level using `Pan...` gene naming
- `*.refer.bed`: BED entries for the cluster hub genes used as the reference set
- `*_gdna.refer.fasta`: genomic sequences of cluster hub genes
- `*_cdna.refer.fasta`: cDNA sequences reconstructed from the generated GTF
- `*_unrenamed_cdna.refer.fasta`: final cDNA sequences keyed by original representative transcript IDs
- `*.graph.json`: graph package manifest with source identities and alignment provenance
- `*.graph.edges.tsv`: unique directed graph edges
- `*.graph.nodes.tsv`: complete graph node lengths and BED metadata

## Current behavior and caveats

- Main chromosome detection is based on BED column 1
- Any chromosome name not recognized as a main chromosome is treated as a contig
- For `construct`, `pre` outputs keep original IDs while final `last` outputs use `Pan...` names
- For `append`, historical and new varieties are assigned as separate varieties, matching construct; no representative-only `Refer` block is used
- Append graph packages retain the source minimap2 version/options and filtering thresholds so incompatible historical graphs can be identified before reuse
- Graph manifests validate source size and modification time; moving a manifest with its sidecars is supported, but changing a recorded BED, BAM, cDNA, or gDNA source invalidates it

## Development notes

This repository is still under active development. Some workflows have been implemented incrementally, and behavior may continue to evolve as clustering, transcript recovery, and append logic are refined.
