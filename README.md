# PanTrans

[中文](README.zh-CN.md)

## Installation

Run these commands on Linux or in WSL. Conda installs Python and the alignment
tools; pip installs PanTrans and its Python dependencies.

```bash
conda create -n pantrans -c conda-forge -c bioconda \
  python=3.12 pip minimap2=2.30 samtools=1.22.1
conda activate pantrans

git clone --branch codex/cds-graph-utr-aware-dedup \
  https://github.com/changQing111/PanTrans.git
cd PanTrans
python -m pip install .

pantrans --version
minimap2 --version
samtools --version
pantrans --help
```

For an editable development installation, use `python -m pip install -e .`.
TransDecoder is needed only for `--annotate-cds` in the legacy commands.

## `preprocess`: prepare genomes and annotations

Provide uncompressed genome FASTAs, GFF3 annotations containing CDS features, and
variety names in matching order. Replace the example paths with your own files.

```bash
pantrans preprocess \
  --genomes data/CS.fa data/ZM366.fa \
  --annotations data/CS.gff3 data/ZM366.gff3 \
  --varieties CS ZM366 \
  --namespace-ids --flank 3000 \
  --prefix Pan --output prepared
```

Alternatively, use three text files, one path or name per line in matching order.
Genome and annotation paths are relative to their respective list files.

```bash
pantrans preprocess \
  --genome-list genomes.txt \
  --annotation-list annotations.txt \
  --variety-list varieties.txt \
  --namespace-ids --flank 3000 \
  --prefix Pan --output prepared_from_lists
```

`--namespace-ids` prefixes IDs with variety names. Use the same ID policy for
historical and new batches, and retain the complete prepared directories.

## `construct`: build and annotate a CDS graph

```bash
pantrans construct \
  --prepared prepared \
  --output construct_results \
  --threads 32 --workers 8
```

The default annotates `last.cluster` and writes `last.gff3` and `last.cdna.fa`.
Add `--annotate-pre` to also annotate `pre.cluster`. `--threads` controls graph
alignment; `--workers` controls annotation concurrency. Use a new output directory.
The current command prefers CS as the reference and uses wheat Chr1A–Chr7D as
main chromosomes.

To rerun annotation using the existing graph and the same prepared inputs:

```bash
pantrans construct \
  --prepared prepared \
  --reuse-graph construct_results \
  --output reannotated_results \
  --workers 8
```

## `append`: add new varieties

Prepare only the new varieties, then provide the historical results and their
original prepared inputs. New variety names and gene/transcript IDs must not
overlap with the historical batch.

```bash
pantrans preprocess \
  --genomes data/JM22.fa \
  --annotations data/JM22.gff3 \
  --varieties JM22 \
  --namespace-ids --flank 3000 \
  --prefix Pan --output new_prepared

pantrans append \
  --history-results construct_results \
  --history-prepared prepared \
  --query-prepared new_prepared \
  --output append_results \
  --threads 32 --workers 8
```

For another append, prepare the next batch as `next_prepared` and use the previous
append's `merged_prepared` as the historical input:

```bash
pantrans append \
  --history-results append_results \
  --history-prepared append_results/merged_prepared \
  --query-prepared next_prepared \
  --output next_append_results \
  --threads 32 --workers 8
```

## `cnv`: screen nearby same-variety gene copies

Use pre clusters and the complete gene BED from the matching dataset.
BED column 4 contains gene IDs; column 5 contains variety names.

```bash
pantrans cnv \
  --pre-cluster construct_results/pre.cluster \
  --bed prepared/Pan_genes.bed \
  --reference CS \
  --max-gap 1000000 \
  --output cnv.genes.txt
```

## `translocation`: screen candidate chromosome changes

Create `main_chroms.txt` with one chromosome name per line, matching BED column 1.
Use clusters, edges and the complete BED from the same dataset.

```bash
pantrans translocation \
  --pre-cluster construct_results/pre.cluster \
  --bed prepared/Pan_genes.bed \
  --edges construct_results/graph.edges.tsv \
  --main-chroms main_chroms.txt \
  --reference CS \
  --min-anchor-genes 3 --max-gap 1000000 --max-skipped-genes 5 \
  --output-prefix screening
```

Outputs are `screening.translocation.genes.tsv` and
`screening.translocation.blocks.tsv`. Add `--allow-reverse` for blocks in reverse
orientation, or `--include-unplaced` to report unplaced gene pairs.

## `legacy-construct`: build with the cDNA workflow

```bash
pantrans legacy-construct \
  --cdna prepared/Pan_mrna.fa \
  --gdna prepared/Pan_gene3000bp_gdna.fa \
  --bed prepared/Pan_genes.bed \
  --gene-regions prepared/Pan_gene_regions.tsv \
  --exon-regions prepared/Pan_exon_regions.tsv \
  --reference CS --prefix Pan \
  --output legacy_results --threads 32
```

## `legacy-append`: extend a legacy cDNA graph

Use the historical unrenamed GFF3/cDNA pair and its graph package. Prepare
`combined_pre_and_new.bed` with the historical pre-cluster representatives and all
new genes; IDs and variety labels must match the sequence inputs.

```bash
pantrans legacy-append \
  --cdna legacy_results/Pan_unrenamed_cdna.refer.fasta \
  --history-gff3 legacy_results/Pan_unrenamed.gff3 \
  --history-graph legacy_results/Pan.graph.json \
  --query-cdna new_prepared/Pan_mrna.fa \
  --gdna new_prepared/Pan_gene3000bp_gdna.fa \
  --bed combined_pre_and_new.bed \
  --gene-regions prepared/Pan_gene_regions.tsv new_prepared/Pan_gene_regions.tsv \
  --exon-regions prepared/Pan_exon_regions.tsv new_prepared/Pan_exon_regions.tsv \
  --prefix Append --output legacy_append_results --threads 32
```

For optional CDS/UTR prediction, install TransDecoder separately and add
`--annotate-cds --transdecoder-dir /path/to/TransDecoder` to either legacy command.
Use `pantrans <subcommand> --help` to list all options for a command.
