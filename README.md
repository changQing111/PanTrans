# PanTrans

[中文](README.zh-CN.md)

## Installation

Run these commands on Linux or in WSL. Conda installs Python and the alignment
tools; pip installs PanTrans and its Python dependencies.

```bash
conda create -n pantrans -c conda-forge -c bioconda \
  python=3.12 pip minimap2=2.30 samtools=1.22.1
conda activate pantrans

git clone --branch codex/longest-cds-synteny \
  https://github.com/changQing111/PanTrans.git
cd PanTrans
python -m pip install .

pantrans --version
minimap2 --version
samtools --version
pantrans --help
```

For an editable development installation, use `python -m pip install -e .`.
The current CDS workflow does not require TransDecoder.

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

Preprocessing writes unflanked and flanked DNA for all genes, full cDNA, all and
valid CDS, plus `Pan_cds.longest.valid.fa` (one longest valid CDS per gene) and
`Pan_cds.longest.valid.tsv` (one audit row per gene). Equal-length valid candidates
are resolved by ascending output transcript ID. Genes without valid CDS or
transcripts retain DNA and coordinates with an explicit selection status. The
construct and append both query these longest valid CDS. See [output details](docs/preprocess-cds-inputs.md).

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

Construct now requires the longest-CDS outputs from preprocessing and queries
each variety's complete unflanked gene bodies separately, including within-variety
comparisons. Graph mapping uses `-f 20000 -N 2000` with the existing splice:hq,
forward, split-index settings. Target tasks run sequentially within the shared
`--threads` budget. Annotation still considers all eligible transcript isoforms.
Per-target BAMs, edges, commands and saturation reports are in `targets/0001/`
etc.; `graph.mapping.json` maps directories to varieties. Construct also exports
`last.gdna.fa`, `last.complete.cdna.fa` and `validation.json`, checking every
reconstructed cDNA. Representatives lacking transcript models are listed explicitly.
See [the new construct workflow](docs/construct-longest-cds.md).

Append uses the same longest-CDS, per-variety N2000 mapping profile. It reuses
historical-to-historical edges, maps new queries to each historical target, and
maps all queries to each new target. Old all-valid-CDS N100 histories require a
fresh construct before append; their mapping profiles cannot be mixed.

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

## `synteny`: pairwise gene-order synteny

Use matching complete prepared inputs, pre clusters and directed graph edges.
All unordered variety pairs are analyzed for forward and reverse blocks.

```bash
pantrans synteny \
  --prepared prepared --results construct_results \
  --max-copy 3 --distance 20 --min-anchors 5 --lift-distance 10 \
  --output synteny_results
```

For append output, use its `merged_prepared/` directory. Large runs can use staged
`--stage index/shard/pair/collect`; select pair tasks with `--pair-index`. The
command reads the existing graph without realignment or changing pre/last
clusters. Outputs include blocks, anchors, positional support, complete binary
candidates and summaries. See [synteny details](docs/synteny.md).
The `translocation` command has been removed; historical report readers remain.

Use `pantrans <subcommand> --help` for command options. The cDNA legacy commands
`legacy-construct` and `legacy-append` have been removed. Current commands are
`preprocess`, `construct`, `append`, `cnv` and `synteny`.

## Acknowledgements

OpenAI Codex provided AI-assisted support for code development, testing and
documentation. The project authors retain responsibility for reviewing the
contributions, the released code and its scientific use.
