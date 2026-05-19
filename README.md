# PanTrans

PanTrans is a command-line toolkit for constructing a pan-genome and pan-transcriptome from multiple varieties in polyploid organisms.

The current implementation builds gene-level relationship graphs from cDNA-to-gDNA alignments, derives `pre` and `last` gene clusters, performs transcript deduplication from filtered BAM alignments, and writes both intermediate and final reference sequence sets.

## Features

- Build pan-gene and pan-transcript clusters from multiple varieties with `construct`
- Reuse an existing BAM file with `construct --bam`
- Control the definition of main chromosomes with `construct --main-chroms`
- Produce both `pre` and `last` outputs:
  `pre.tmp.cluster`, `last.tmp.cluster`, `pre.gtf`, `gtf`, reference FASTA files, and reference BED files
- Append new varieties to an existing reference set with `append`
- Rename final `last`-level genes in GTF output using `Pan...` identifiers while preserving original IDs in `pre` outputs

## Installation

Install in editable mode from the project root:

```bash
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
- `test_pantrans.gtf`
- `test_pantrans_cdna.refer.fasta`
- `test_pantrans_gdna.refer.fasta`
- `test_pantrans.refer.bed`
- alignment BAM files

### Notes on naming

- `pre.gtf` preserves the original cluster hub gene IDs and transcript IDs
- final `last`-level GTF output renames genes to `Pan...` identifiers such as `Pan1A000001`

## Append

Use `append` to merge a new query set with an existing reference set and re-run graph construction and transcript processing on the merged data.

### Command

```bash
pantrans append \
  --name <new_variety_file_or_list> \
  --cdna <query_cdna.fasta> \
  --gdna <query_gdna.fasta> \
  --bed <query.bed> \
  --refer_cdna <refer_cdna.fasta> \
  --refer_gdna <refer_gdna.fasta> \
  --refer_bed <refer.bed> \
  --prefix <prefix> \
  --output <output_dir>
```

### Example

```bash
pantrans append \
  --name JM22 \
  --cdna test/JM22_cdna.fasta \
  --gdna test/JM22.gdna.fasta \
  --bed test/JM22_new_dedup.bed \
  --refer_cdna test/pantrans_construct_out/test_pantrans_pre_cdna.refer.fasta \
  --refer_gdna test/pantrans_construct_out/test_pantrans_pre_gdna.refer.fasta \
  --refer_bed test/pantrans_construct_out/test_pantrans_pre.refer.bed \
  --prefix test_append_JM22 \
  --threads 32 \
  --output test/pantrans_append_out
```

### Append workflow

The current `append` implementation performs the following steps:

1. Merge query and reference `cdna`, `gdna`, and `bed` inputs
2. Align merged `cdna` to merged `gdna`
3. Filter the merged BAM
4. Rebuild graph-based `pre` and `last` clusters from the merged alignment graph
5. Write `pre` and `last` GTF/FASTA/BED outputs for the merged result

Important:
`append` currently follows a re-clustering strategy. It does not preserve the original reference `pre` clusters as fixed units. New query sequences can bridge previously separate reference clusters, which may reduce the number of `pre` clusters after append.

## Output interpretation

- `pre.tmp.cluster`: graph-derived intermediate clustering before final recursive refinement
- `last.tmp.cluster`: final cluster list after recursive splitting and chromosome-aware assignment
- `pre.gtf`: transcript deduplication result at the `pre` cluster level using original hub gene IDs
- `gtf`: transcript deduplication result at the `last` cluster level using `Pan...` gene naming
- `*.refer.bed`: BED entries for the cluster hub genes used as the reference set
- `*_gdna.refer.fasta`: genomic sequences of cluster hub genes
- `*_cdna.refer.fasta`: cDNA sequences reconstructed from the generated GTF

## Current behavior and caveats

- Main chromosome detection is based on BED column 1
- Any chromosome name not recognized as a main chromosome is treated as a contig
- For `construct`, `pre` outputs keep original IDs while final `last` outputs use `Pan...` names
- For `append`, the reference set is treated as a merged reference block during cluster assignment rather than as separately fixed reference varieties

## Development notes

This repository is still under active development. Some workflows have been implemented incrementally, and behavior may continue to evolve as clustering, transcript recovery, and append logic are refined.