> Historical cDNA workflow reference. PanTrans 0.8 uses CDS by default; see [current CLI](README.md). The historical development notes below are not the current CDS command contract.

# PanTrans

[中文版](README.zh-CN.md)

PanTrans is a command-line toolkit for constructing a pan-genome and pan-transcriptome from multiple varieties in polyploid organisms.

The current implementation builds gene-level relationship graphs from cDNA-to-gDNA alignments, derives `pre` and `last` gene clusters, performs transcript deduplication from filtered BAM alignments, and writes both intermediate and final reference sequence sets.

## Features

Version **0.7.0** switches all new annotations to GFF3 and supports TransDecoder CDS/UTR prediction.

- Screen pre clusters for nearby gene copies with `cnv` and candidate chromosome changes with `translocation`

- Build pan-gene and pan-transcript clusters from multiple varieties with `legacy-construct`
- Reuse an existing BAM file with `construct --bam`
- Control the definition of main chromosomes with `construct --main-chroms`
- Persist a reusable complete directed graph package from `legacy-construct`
- Produce both `pre` and `last` outputs:
  `pre.tmp.cluster`, `last.tmp.cluster`, `pre.gff3`, `gff3`, reference FASTA files, and reference BED files
- Append new varieties to an existing reference set with `legacy-append`
- Retain both original representative IDs and official `Pan...` IDs in final GFF3/cDNA outputs

## Translocation visualization

The standalone `scripts/plot_translocations.py` draws chromosome overviews and local
gene correspondences with cubic Bezier connections. It supports `--flank`, optional
ambiguous matches, and PDF/PNG/SVG export. Install Matplotlib separately:

```bash
pip install matplotlib
python scripts/plot_translocations.py --help
```

See [plotting usage and input formats](scripts/plot_translocations.md). This script
does not rerun detection; plotted spans are candidate anchor spans, not breakpoints.

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

### Optional CDS / UTR annotation

All generated annotations, including pre, both final variants and rescue intermediates,
now use `.gff3`. Files contain `##gff-version 3` and a gene → mRNA → exon/CDS/UTR
hierarchy with escaped `ID`/`Parent` attributes. Official outputs rename IDs and
Parent references at every level while retaining the original representative
sequence name in column 1 and the same gene-relative coordinates.
Append uses `--history-gff3`, which also reads legacy historical GTF;
`--history-gtf` remains a compatibility alias. Existing files are not automatically
converted or removed; new runs no longer generate GTF.
`scripts/repair_official_ids.py` now repairs the new `.gff3` output pairs.

Both `legacy-construct` and `legacy-append` accept `--annotate-cds`. After final transcript
deduplication and rescue, TransDecoder predicts one primary ORF per transcript.
PanTrans adds `CDS`, `five_prime_UTR`, `three_prime_UTR`, and supported
`start_codon` / `stop_codon` features to both `<prefix>_unrenamed.gff3` and
`<prefix>.gff3`. Pre-cluster outputs are not annotated.

Add these options to an existing command:

```bash
--annotate-cds --transdecoder-dir /path/to/TransDecoder --min-orf-length 50
```

The executable directory is optional when TransDecoder is on `PATH`; specifying
it also enables annotation. Classic `TransDecoder.LongOrfs` / `TransDecoder.Predict`
and the newer `TransDecoder` wrapper are supported, with their dependencies in a
working Linux/WSL environment. Annotation is disabled by default. The minimum
candidate ORF length defaults to 50 amino acids.

TransDecoder scoring remains serial by default, regardless of `-t`. Add
`--transdecoder-parallel -t 8` to use up to eight scoring processes (`-t 1`
keeps the original serial path). Parallel scoring currently supports classic
TransDecoder 5.7.1 with the expected scoring command; other installations log a
warning and fall back to serial execution. All chunks use the same trained model
and original scorer, and scores are merged in input order. LongOrfs, training,
and downstream prediction stay serial. The installed TransDecoder is not modified;
a private Predict adapter and per-chunk diagnostics are retained in the run folder.


Prediction runs once on final unrenamed cDNA, using sense-strand mode and
`--single_best_only`. Coordinates are projected through exons onto the existing
representative gene sequences, then official IDs are generated from the same GFF3.
GFF3 CDS includes the stop codon, matching TransDecoder; the 3' UTR begins
after it. Spliced CDS segments share one ID and have appropriate phase.
Optional start/stop codon features overlap CDS and have phase `.`. Partial ORFs only
receive UTRs at supported boundaries. Transcripts without predicted ORFs retain
their transcript/exon models and are not automatically classified as noncoding.

Logs and raw predictions (GFF3, peptide and CDS FASTA with original transcript/ORF
IDs) are retained in a fresh `<prefix>_transdecoder/run-*/` directory on every run.
Failures or invalid predictions raise errors instead of silently omitting annotation.
Append recomputes predictions when enabled; historical CDS/UTR features are not
inherited, so enable the option for each append that needs coding annotation.
Transcript strands in local GFF3 come from the selected BAM alignment relative to
the input gDNA, not chromosome BED column 6. cDNA extraction joins exons in coordinate
order and reverse-complements only for GFF3 `-`; unknown transcript strands fail.
Already oriented gDNA is therefore not reverse-complemented twice. Opposite-strand
models are deduplicated separately; a gene containing both strands has gene strand
`.` and explicitly stranded transcripts/children. Previously generated incorrect
strands are not automatically repaired: reuse the BAM to regenerate annotations
and cDNA, then rerun TransDecoder. Repair historical annotations before using append.

These are ORF-based predictions using the standard genetic code, not experimentally
confirmed translation boundaries.

PanTrans expects the following input types:

- `cdna` FASTA: transcript sequences
- `gdna` FASTA: genomic gene sequences used as alignment targets
- `bed` file: gene coordinates and strands
- BED column 5: the variety name for each gene
- variety names are inferred from BED column 5

The BED file must include gene IDs in column 4, variety names in column 5, and chromosome or contig names in column 1. Column 6 may provide the strand (or is treated as `.` when omitted).

The column-5 label, rather than a gene-ID prefix, defines the variety for every
construct and append operation. Both commands derive their variety list from
the BED and do not accept a separate variety-name argument.

### Required BED format

The BED file must be tab-delimited and contain at least five columns on every
data row. Coordinates use the standard BED convention: 0-based, half-open
intervals (`start` inclusive, `end` exclusive). Do not add a header row. Blank
lines and lines beginning with `#` are ignored.

| Column | Required value |
| --- | --- |
| 1 | Chromosome or contig name |
| 2 | Integer start coordinate |
| 3 | Integer end coordinate |
| 4 | Unique gene ID, matching a gDNA FASTA record and the gene part of cDNA transcript IDs |
| 5 | Non-empty variety name for that gene; do not use `0`, `.`, or an empty value |
| 6 | Strand, `+` or `-`; optional, defaults to `.` when omitted |

Example:

```text
A1	277424819	277434365	gene_0001	CS	+
A1	154737879	154739482	gene_0002	JM22	-
```

For `legacy-construct`, the BED contains every input gene. For `legacy-append`, the combined
BED contains the historical representative genes required by the graph package
plus every query gene; column 5 must label both the historical and query rows.

## Construct

Use `legacy-construct` to build a new pan-genome and pan-transcriptome from multiple varieties.

### Command

```bash
pantrans legacy-construct \
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
- `--threads`: thread budget shared by minimap2 and samtools
- `--rescue`: explicitly enable rescue of cluster representatives missing from the generated GFF3; disabled by default

Rescue first reuses available pre-cluster GFF3 models, then runs one pooled
cDNA-to-gDNA alignment over all members of the remaining missing clusters per
output stage. Only alignments within each cluster are used to recover models.
The original mapping parameters are retained; some genes may remain missing.

### Example

```bash
pantrans legacy-construct \
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
- `test_pantrans_pre.gff3`
- `test_pantrans_pre_cdna.refer.fasta`
- `test_pantrans_pre_gdna.refer.fasta`
- `test_pantrans_pre.refer.bed`
- `test_pantrans_unrenamed.gff3`
- `test_pantrans_unrenamed_cdna.refer.fasta`
- `test_pantrans.gff3`
- `test_pantrans_cdna.refer.fasta`
- `test_pantrans_gdna.refer.fasta`
- `test_pantrans.refer.bed`
- `test_pantrans.graph.json`
- `test_pantrans.graph.edges.tsv`
- `test_pantrans.graph.nodes.tsv`
- alignment BAM files

### Notes on naming

- `pre.gff3` preserves the original cluster hub gene IDs and transcript IDs.
- `<prefix>_unrenamed.gff3` and `<prefix>_unrenamed_cdna.refer.fasta`
  retain the final `last` representative IDs used by the graph, BED, and cluster files.
- `<prefix>.gff3` and `<prefix>_cdna.refer.fasta` contain the same final
  transcript models renamed to official `Pan...` IDs such as `PanChr1A000001`.
  IDs preserve the representative gene's exact BED chromosome name after `Pan`.

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

Use `legacy-append` to reuse a historical construct graph and add only the new-to-history alignment blocks. The merged graph then follows the same clustering and transcript processing flow as `legacy-construct`.

### Command

```bash
pantrans legacy-append \
  --cdna <previous_unrenamed_final_cdna.fasta> \
  --history-gff3 <previous_unrenamed_final.gff3> \
  --query-cdna <new_variety_cdna.fasta> \
  --gdna <query_gdna.fasta> \
  --bed <combined_representatives_and_query.bed> \
  --history-graph <historical_prefix.graph.json> \
  --prefix <prefix> \
  --output <output_dir>
```

The three transcript inputs have distinct roles:

- `--cdna` is the previous run's `<prefix>_unrenamed_cdna.refer.fasta`.
- `--history-gff3` is its matching `<prefix>_unrenamed.gff3`.
- `--query-cdna` is the new variety's complete cDNA FASTA.

The historical GFF3 and cDNA must contain exactly the same transcript IDs. Their
gene IDs must be unrenamed IDs present in the historical graph package; passing
the official `Pan...` GFF3 is rejected. `--bed` is one combined BED containing
the first gene ID from every historical `pre.tmp.cluster` plus every gene from
the new variety. The graph package supplies complete historical gDNA, BED,
filtered BAM, graph edges, node lengths, variety order, reference, and alignment
provenance.

For restartability, append also accepts `--query-to-all-bam` and
`--history-to-query-bam` to reuse either completed cross-alignment BAM.
Pass `--rescue` to enable rescue of cluster representatives missing from the
generated GFF3; rescue is disabled by default.

### Example

```bash
pantrans legacy-append \
  --cdna test/pantrans_construct_out/test_pantrans_unrenamed_cdna.refer.fasta \
  --history-gff3 test/pantrans_construct_out/test_pantrans_unrenamed.gff3 \
  --query-cdna test/JM22_cdna.fasta \
  --gdna test/JM22.gdna.fasta \
  --bed test/test_pantrans_pre.refer_append_JM22.bed \
  --history-graph test/pantrans_construct_out/test_pantrans.graph.json \
  --prefix test_append_JM22 \
  --threads 32 \
  --output test/pantrans_append_out
```

### Append workflow

The `legacy-append` workflow performs the following steps:

1. Validate the historical graph package and its source-file identities
2. Build a working cDNA from historical nonredundant plus query transcripts,
   and build full historical-plus-query gDNA and BED files
3. Align new cDNA to all historical-plus-new gDNA
4. Align historical nonredundant cDNA to new gDNA
5. Reuse historical graph edges, merge the two cross-edge blocks, and merge their filtered BAMs with the historical filtered BAM
6. Rebuild graph-based `pre` and `last` clusters using exactly the construct clustering flow
7. Seed final transcript deduplication with historical GFF3 models whose original
   representative gene ID is still a current `last` representative
8. Write official and unrenamed final GFF3/cDNA outputs plus a graph package for
   a later append

All sequence-backed historical and new genes are eligible cluster members. The append BED is only the user-facing representative-plus-query BED contract; the historical graph package supplies the complete historical BED and gDNA sequence set used internally. Historical GFF3 coordinates are reused only when their original representative gene ID remains a current `last` cluster representative; transcript models are not projected onto a different representative coordinate system.

The two cross-alignment blocks avoid recomputing historical cDNA against historical gDNA. Since minimap2 retains a bounded number of secondary alignments, the incremental graph is expected to be highly similar to, but not byte-identical with, a monolithic construct graph.

To append another variety, use the current append run's graph package together
with its `_unrenamed.gff3` and `_unrenamed_cdna.refer.fasta` outputs as the next
run's historical inputs.

If an append job fails after either cross-alignment finishes, rerun it with the
corresponding `--query-to-all-bam` and/or `--history-to-query-bam` path. The BAM
must contain the minimap2 `@PG` record with the same edge-generating options and
version. PanTrans also verifies the primary query IDs and sequence identities
against the current cDNA FASTA before re-filtering the reused BAM with the
current thresholds.

## CNV screening

Screen existing pre clusters without rerunning alignment or clustering:

```bash
pantrans cnv \
  --pre-cluster five_variety_construct_out/five_variety_pre.tmp.cluster \
  --bed five_variety_genes.bed \
  --reference CS \
  --output five_variety.cnv.genes.txt
```

A cluster qualifies when two or more distinct genes from the same BED column-5
variety lie on the same BED chromosome with interval gaps at most 1,000,000 bp.
The gap is `max(0, next_start - current_group_max_end)` for sorted BED intervals;
overlaps have zero gap. Chains are transitive: A–B and B–C within the threshold
belong to one group even if A–C exceeds it. Strand is not a restriction.
Detection stops at the first qualifying group in each cluster; it does not
require different copy counts between varieties.

The output is a plain text file with **one original gene ID per hit pre cluster**,
without a header, in input cluster order. Representatives are chosen from the
**entire cluster**, preferring the `-r` variety, then the longest BED interval
(`end-start`), then the lexicographically smallest gene ID. The representative
need not belong to the first detected group. If the reference variety is absent
from a particular cluster, the longest member wins. The reference must exist in
the BED overall. No hits produces an empty file.

`--max-gap` can override the default in bp. The BED must cover every cluster
member and use 0-based, half-open coordinates and valid variety labels. Distinct
contigs must have distinct BED chromosome names; a shared `contig` label cannot
establish proximity on separate unplaced sequences. This command applies the
specified proximity rule, not a between-variety copy-number comparison.

## Candidate translocations

Compare each non-reference variety with `-r` using existing pre clusters, the
complete BED, and the directed graph edge sidecar from the same run:

```bash
pantrans translocation \
  --pre-cluster five_variety_construct_out/five_variety_pre.tmp.cluster \
  --bed five_variety_genes.bed \
  --edges five_variety_construct_out/five_variety.graph.edges.tsv \
  --main-chroms main_chroms \
  --reference CS \
  --output-prefix five_variety
```

No alignment, BAM loading, or re-clustering is performed. BED column 5 defines
varieties. All cluster and edge genes must exist in the complete BED; pre-cluster
membership must be unique. The edge TSV has two columns (`source_gene`,
`target_gene`), no header. Both directed edges must exist **within the same pre
cluster**; self edges, indirect paths and one-way matches do not support a pair.
Duplicate edges do not increase support. Clusters without reference genes are
counted in the log and skipped.

`--main-chroms` is required: use one exact BED chromosome label per line. Names
must be consistent between varieties; `Chr1A` and `1A` are not normalized by this
command. By default, pairs with either gene outside this list are excluded from
reports AND reciprocal-partner/ambiguity counts. Use `--include-unplaced` to
include these pairs as `unplaced`; they never anchor translocation blocks.
Same-chromosome reciprocal pairs are
considered normal and omitted from the gene report. A cross-chromosome pair is
`ambiguous_copy` if either endpoint also has a same-chromosome reciprocal match
or multiple reciprocal partners in this reference/query-variety comparison.
Only with `--include-unplaced` does an alternative unplaced partner count toward copy ambiguity. Ambiguous
pairs are retained for inspection but cannot anchor blocks.

Only one-to-one reciprocal cross-chromosome pairs anchor blocks. Anchors must
share the query variety and both chromosome names. Their start coordinates must
increase on both sides (or decrease on the query side with `--allow-reverse`).
BED intervals may overlap, but tied starts provide no gene-order evidence.

| Option | Default | Meaning |
| --- | --- | --- |
| `--min-anchor-genes` | `3` | Minimum **distinct pre clusters** supporting a block; minimum allowed value is 2 |
| `--max-gap` | `1000000` | Maximum interval gap in bp between successive anchors on **each** side; overlaps have gap 0 |
| `--max-skipped-genes` | `5` | Maximum intervening genes on each side, counted from the **complete BED**, including non-anchors |
| `--allow-reverse` | off | Also accept reversed query gene order; orientation is based on gene order, not BED strand |
| `--include-unplaced` | off | Include unplaced pairs in the gene report and ambiguity counts, but never as block anchors |

Chaining uses deterministic longest-path proposals within these limits, then
greedily selects disjoint chains. An anchor belongs to at most one block. This
does not enumerate all possible rearrangements. Gaps constrain adjacent anchors,
not total block length. Many pairs from one family/pre cluster cannot by
themselves produce a supported block.

Outputs (the `-o` / `--output-prefix` argument is a **prefix**, not a directory):

- `five_variety.translocation.genes.tsv`: original gene IDs, both BED intervals,
  varieties, reciprocal partner counts, same-chromosome alternatives, status,
  reason, and block ID. `pre_000123` refers to physical line 123 in the input
  cluster file. Status is `supported_block`, `single_gene_candidate`,
  `ambiguous_copy`, or `unplaced`.
- `five_variety.translocation.blocks.tsv`: candidate spans, orientation, anchor
  and independent-cluster counts, and ordered supporting genes. Chains with at
  least two pairs but too few distinct clusters are `low_support_block`; their
  gene rows retain `single_gene_candidate`. Isolated pairs have no block (`.`).

All output coordinates are BED 0-based half-open. Block spans cover supporting
genes and are **not breakpoint estimates**. Empty results still contain TSV
headers. Both reports are staged before publication; detected publication
failures restore prior outputs. The log records thresholds and screening counts.
`supported_block` denotes collinear evidence for a candidate position difference
relative to the reference, not proof of a translocation or its evolutionary
direction. Collinearity does not resolve all copy, assembly or annotation errors.

## Output interpretation

- `pre.tmp.cluster`: graph-derived intermediate clustering before final recursive refinement
- `last.tmp.cluster`: final cluster list after recursive splitting and chromosome-aware assignment
- `pre.gff3`: transcript deduplication result at the `pre` cluster level using original hub gene IDs
- `_unrenamed.gff3`: final `last` transcript models using original representative IDs
- `gff3`: transcript deduplication result at the `last` cluster level using `Pan...` gene naming
- `*.refer.bed`: BED entries for the cluster hub genes used as the reference set
- `*_gdna.refer.fasta`: genomic sequences of cluster hub genes
- `*_cdna.refer.fasta`: cDNA sequences reconstructed from the generated GFF3
- `*_unrenamed_cdna.refer.fasta`: final cDNA sequences keyed by original representative transcript IDs
- `*.graph.json`: graph package manifest with source identities and alignment provenance
- `*.graph.edges.tsv`: unique directed graph edges
- `*.graph.nodes.tsv`: complete graph node lengths and BED metadata

## Current behavior and caveats

- Main chromosome detection is based on BED column 1
- Any chromosome name not recognized as a main chromosome is treated as a contig
- For `legacy-construct`, `pre` outputs keep original IDs while final `last` outputs use `Pan...` names
- For `legacy-append`, historical and new varieties are assigned as separate varieties, matching construct; no representative-only `Refer` block is used
- Append graph packages retain the source minimap2 version/options and filtering thresholds so incompatible historical graphs can be identified before reuse
- Graph manifests validate source size and modification time; moving a manifest with its sidecars is supported, but changing a recorded BED, BAM, cDNA, or gDNA source invalidates it

## Development notes

This repository is still under active development. Some workflows have been implemented incrementally, and behavior may continue to evolve as clustering, transcript recovery, and append logic are refined.

For supported TransDecoder 5.7.1, both serial and parallel modes use a private deterministic ORF selector: original priorities are preserved, with the lexicographically smallest ORF ID winning complete ties. This can change previously arbitrary tied predictions. Unsupported installations warn that deterministic tie-breaking is unavailable.

## Prepare construct/append inputs

`pantrans preprocess` accepts parallel lists of uncompressed genome FASTAs, GFF3
annotations (optionally gzip-compressed), and unique variety names:

```bash
pantrans preprocess \
  --genomes /data/CS.fa /data/ZM366.fa \
  --annotations /data/CS.gff3 /data/ZM366.gff3 \
  --varieties CS ZM366 \
  --flank 3000 --prefix two_variety --output prepared
```

For long lists, use `--genome-list genomes.txt --annotation-list annotations.txt
--variety-list varieties.txt` instead. Each file contains one entry per line in
matching order; blank/comment lines are ignored. Relative paths in a list file
are resolved relative to that file, not the process working directory.

Outputs include `two_variety_gene3000bp_gdna.fa`, `two_variety_mrna.fa`,
`two_variety_genes.bed`, FASTA indexes, `_gene_regions.tsv`, `_exon_regions.tsv`,
`_id_map.tsv`, and `_preprocess.json`. BED column 5 is the supplied variety name.
Genes without transcript models are skipped and counted. Transcript models must
have explicit exons; CDS-only annotations and GTF are not supported by this
subcommand. Exons include UTR sequence. Missing contigs, invalid coordinates,
overlapping exons and unresolved parents are errors, not silently clipped models.

Flanks are clipped to chromosome/contig bounds, including contigs shorter than
6 kb; **gene/exon coordinates themselves are never clipped**. No padding is added.
Both gDNA and spliced cDNA are oriented in the gene's 5′→3′ direction. Tables use
0-based half-open coordinates; FASTA headers use 1-based inclusive coordinates.
The gene table records actual upstream/downstream lengths and the gene body in
local gDNA coordinates; the exon table records each transcript's local blocks.
These tables support subsequent position-aware masking/filtering; preprocessing
itself does not mask sequences or change the current BAM filter.

Gene IDs are preserved and must be globally unique. If different varieties reuse
IDs, pass `--namespace-ids` consistently to all related preprocessing runs; gene
IDs become `VARIETY__original`. Compatible transcript IDs (`GENE.number`) retain
their suffix; other transcript IDs are deterministically normalized to that form.
Original IDs are always recorded in `_id_map.tsv`. Output files are not overwritten;
use a fresh prefix or output directory. Input genome indexes are built privately
and do not require writing beside the genome FASTA.

Use the three combined files directly for `legacy-construct`. For `legacy-append`, preprocess
only the new varieties, pass their `_mrna.fa` as `--query-cdna` and gDNA as `--gdna`,
and concatenate the historical BED with the new BED for `--bed`. Keep original
historical IDs and variety labels; the historical final cDNA/GFF3 and graph package
are still supplied through the existing append options. Keep main chromosome
names explicit through `construct --main-chroms`; preprocessing preserves contig
names and does not guess which sequences are main chromosomes.

## Exon-aware mapping and filtering

Pass preprocessing tables to `legacy-construct` or `legacy-append`:

```bash
pantrans legacy-construct --cdna prepared/five_variety_mrna.fa \
  --gdna prepared/five_variety_gene3000bp_gdna.fa \
  --bed prepared/five_variety_genes.bed \
  --gene-regions prepared/five_variety_gene_regions.tsv \
  --exon-regions prepared/five_variety_exon_regions.tsv \
  --exon-overlap 0.90 --exon-query-overlap 0.80 \
  --reference CS --prefix five_variety --output results --threads 32
```

These options activate coordinate-preserving N masking outside each gene body
for fresh minimap2 mappings, followed by same-transcript exon overlap and strand
filtering. Introns remain unmasked. Query lengths come from input cDNA, including
secondary records without SEQ. Self-alignments follow the same quality/position
rules. Original gDNA remains the source for output sequences. Each filtered BAM
gets a `.filter.json` count report. Without region tables, legacy behavior remains.

`construct --bam` can apply the new filter to old raw BAMs, but cannot undo the
original candidate truncation. Such graphs are marked `provided-raw-bam` and are
not interchangeable with masked-reference graphs. Exon-aware append currently
requires fresh cross-alignments and tables covering all historical and new genes;
historical graph policy, thresholds, masking mode and historical per-gene region
digests must match. Rescue uses the same policy. Coordinates are in preprocessing
reference orientation; never pass genomic exon coordinates directly.

### CDS and gene-body preprocessing outputs

`preprocess` also emits `<prefix>_gene.fa` (gene body with introns, no flanks),
`<prefix>_cds.fa` (all annotated CDS), `<prefix>_cds.valid.fa` (CDS passing
complete-ORF/phase/sequence checks), CDS coordinates, per-transcript status and
`_transcript_models.jsonl`. These are input checks, **not alignment or graph
validation**. Unqualified annotated CDS remain in `_cds.fa`; genes lacking CDS
remain in the cDNA/gDNA/BED outputs. Empty CDS files do not receive a FASTA index.

Schema-2 preprocessing directories can be read directly by
`experimental_cds_construct.py --prepared prepared --output results` without
re-reading original annotations. Filenames are resolved from the manifest and
may use any preprocessing prefix. Legacy schema-1 directories still require
`--annotations`. CDS phase validates the reading frame; internal phase bases
are never trimmed from the concatenated CDS. Existing flank clipping and
strand orientation remain unchanged.

## CDS graph development workflow

The current CDS-graph/UTR-aware workflow is available through experimental scripts,
including the accelerated annotation-only restart. See [workflow and limitations](docs/cds-workflow.md).
The default construct/append commands have not been replaced by the CDS workflow.

### Run a CDS-based construct

First prepare the inputs as above, then run the experimental constructor:

```bash
python scripts/experimental_cds_construct.py \
  --prepared prepared/five_variety \
  --output results/five_variety_cds \
  --threads 32 --workers 8
```

This maps `<prefix>_cds.valid.fa` against `<prefix>_gene.fa` (gene body with
introns, no 3 kb flanks), filters coding-supported directed edges, and writes
`pre.cluster` and `last.cluster`. It then maps full `<prefix>_mrna.fa`, including
UTRs, to each cluster representative's unmasked `<prefix>_gene3000bp_gdna.fa`
for transcript models and UTR-aware deduplication. The command takes the
preprocessing directory so it can load its coordinate/model manifest. It is an
experimental CDS workflow, separate from `pantrans legacy-construct`; `--annotations`
is only needed for legacy schema-1 preprocessing directories.
