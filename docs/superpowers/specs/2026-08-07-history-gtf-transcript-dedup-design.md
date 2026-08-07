# Historical GTF-Aware Append Transcript Deduplication Design

## Goal

Preserve the representative-ID agreement already achieved by incremental append while making final transcript deduplication reuse the nonredundant transcript models produced by the previous construct or append run.

## User-Facing Contract

Construct keeps its existing inputs and produces two final transcriptome variants:

- `<prefix>.gtf` and `<prefix>_cdna.refer.fasta`: the existing official outputs with `Pan...` gene and transcript IDs.
- `<prefix>_unrenamed.gtf` and `<prefix>_unrenamed_cdna.refer.fasta`: the same transcript models with original last-cluster representative gene IDs and `<original_gene_id>.<index>` transcript IDs.

Append separates the historical transcriptome from the new-variety cDNA:

```text
--cdna         previous run's <prefix>_unrenamed_cdna.refer.fasta
--history-gtf  previous run's <prefix>_unrenamed.gtf
--query-cdna   new-variety cDNA FASTA
--gdna         new-variety gDNA FASTA
```

The existing `--history-graph`, `--bed`, `--name`, BAM resume options, output prefix, and thread options remain unchanged. This is an intentional append CLI change: `--cdna` now means historical deduplicated cDNA, while `--query-cdna` means the new variety.

## Construct Data Flow

For the final last-cluster output, transcript deduplication first writes an unrenamed GTF. Rescue records are added to this unrenamed GTF using original representative IDs. PanTrans then creates the official GTF by deterministically mapping each original representative gene ID through the existing `last_rename_map` and preserving the transcript's numeric suffix.

Both GTF files are converted independently to cDNA FASTA, so the sequence content and transcript structures are the same while their FASTA IDs follow their respective GTF IDs. Pre-cluster output remains unrenamed and unchanged.

## Append Alignment and Graph Data Flow

Append continues to reuse all historical graph edges and the historical filtered BAM from `--history-graph`. It creates only the two cross-alignment blocks:

1. new-variety `--query-cdna` against merged historical-plus-new gDNA;
2. historical nonredundant `--cdna` against new-variety gDNA.

The working cDNA used for transcript length lookup and rescue is the concatenation of historical nonredundant cDNA and new-variety cDNA. The complete historical raw cDNA recorded in older graph packages is no longer used as the history-to-query alignment input.

## Historical GTF Candidate Merge

Final append transcript deduplication reads transcript models from two sources:

- filtered merged BAM alignments for new and cross-aligned transcript candidates;
- the previous unrenamed final GTF for historical candidates already selected by the prior run.

The GTF parser reconstructs exon coordinates and splice intervals for each transcript. Historical GTF candidates override a BAM candidate with the same transcript ID because the historical GTF is the authoritative previously selected model.

A historical GTF model is directly reusable only when its sequence/representative gene ID is also a key in the current last-cluster dictionary. Coordinates in a gene-relative GTF cannot be transferred to a different representative sequence. Historical representatives that changed in the new clustering are skipped with a count in the log; BAM candidates may still recover them when an alignment to the new representative exists.

The combined candidates then pass through the existing splice-pattern grouping and longest-transcript selection. The final result is first written with original representative IDs and then renamed to the official `Pan...` IDs.

## Input Validation

Append validates the historical GTF and cDNA before running minimap2:

- every GTF transcript feature has `gene_id` and `transcript_id`;
- GTF sequence names and `gene_id` attributes agree and belong to the historical graph gene set;
- every cDNA FASTA transcript ID occurs in the GTF;
- every GTF transcript ID has a cDNA FASTA record;
- transcript IDs map back to their original gene IDs by the normal final-suffix rule;
- new `--query-cdna` transcript gene IDs belong only to the new variety and occur in new gDNA.

These checks reject an official renamed `Pan...` GTF/cDNA pair when an unrenamed pair is required.

## Outputs and Chaining

Append writes the same two final variants as construct. The append `_unrenamed.gtf` and `_unrenamed_cdna.refer.fasta` become the `--history-gtf` and `--cdna` inputs of a later append, allowing repeated incremental additions without restoring complete historical cDNA.

## Error Handling

Missing files, mismatched GTF/cDNA transcript sets, renamed historical IDs, and gene IDs outside the historical graph fail before alignment. Historical GTF representatives that are valid historical genes but are no longer current last-cluster representatives are not fatal; they are counted and skipped because their coordinates are not valid on the new representative sequence.

## Testing

Tests cover:

- parsing a historical GTF into exon and splice models;
- preserving a historical transcript absent from the BAM;
- combining historical and new BAM candidates before splice-pattern deduplication;
- deterministic unrenamed-to-renamed GTF conversion;
- construct production of both GTF/cDNA variants;
- append CLI and function signature changes;
- historical GTF/cDNA validation before alignment;
- use of historical nonredundant cDNA for history-to-query alignment;
- use of the historical GTF only for final last-cluster deduplication;
- existing graph, cluster, provenance, and validation-script tests.

## Non-Goals

- Making incremental clusters byte-identical to a monolithic construct run.
- Re-aligning all original historical cDNA.
- Projecting historical GTF coordinates onto a different representative gene without sequence alignment.
- Changing pre/last graph construction or representative selection.
