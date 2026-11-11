# PanTrans

[English](README.md)

## 安装

在 Linux 或 WSL 中执行以下命令。Conda 安装 Python 和比对工具，pip 安装 PanTrans 及其 Python 依赖。

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

开发时可改用 `python -m pip install -e .`。只有旧流程启用 `--annotate-cds` 时才需要 TransDecoder。

## `preprocess`：准备基因组和注释

按相同顺序提供未压缩的基因组 FASTA、含 CDS 特征的 GFF3 注释和品种名。将示例路径替换为自己的文件。

```bash
pantrans preprocess \
  --genomes data/CS.fa data/ZM366.fa \
  --annotations data/CS.gff3 data/ZM366.gff3 \
  --varieties CS ZM366 \
  --namespace-ids --flank 3000 \
  --prefix Pan --output prepared
```

也可使用三个文本列表，每行一个路径或品种名，行序一一对应。基因组和注释列表中的相对路径以各自列表文件所在目录为基准。

```bash
pantrans preprocess \
  --genome-list genomes.txt \
  --annotation-list annotations.txt \
  --variety-list varieties.txt \
  --namespace-ids --flank 3000 \
  --prefix Pan --output prepared_from_lists
```

`--namespace-ids` 为 ID 添加品种名前缀。历史和新增批次采用相同的 ID 策略，并保留完整的预处理目录。

## `construct`：CDS 构图与注释

```bash
pantrans construct \
  --prepared prepared \
  --output construct_results \
  --threads 32 --workers 8
```

默认注释 `last.cluster`，生成 `last.gff3` 和 `last.cdna.fa`；添加 `--annotate-pre` 可同时注释 `pre.cluster`。
`--threads` 控制构图比对线程数，`--workers` 控制注释并发数。每次使用新的输出目录。
当前命令优先选择 CS 作为参考，主染色体采用小麦 Chr1A–Chr7D。

复用已有图，仅重新注释时，使用与原构图相同的预处理输入：

```bash
pantrans construct \
  --prepared prepared \
  --reuse-graph construct_results \
  --output reannotated_results \
  --workers 8
```

## `append`：追加新品种

只预处理新增品种，再提供历史结果及其原始预处理输入。新增品种名和 gene/transcript ID 均不能与历史批次重叠。

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

再次追加时，将下一批品种预处理到 `next_prepared`，使用上一次 append 生成的 `merged_prepared` 作为历史输入：

```bash
pantrans append \
  --history-results append_results \
  --history-prepared append_results/merged_prepared \
  --query-prepared next_prepared \
  --output next_append_results \
  --threads 32 --workers 8
```

## `cnv`：筛查同品种邻近基因拷贝

提供同一数据集的 pre 聚类和完整基因 BED。BED 第 4 列为 gene ID，第 5 列为品种名。

```bash
pantrans cnv \
  --pre-cluster construct_results/pre.cluster \
  --bed prepared/Pan_genes.bed \
  --reference CS \
  --max-gap 1000000 \
  --output cnv.genes.txt
```

## `translocation`：筛查候选染色体变化

准备 `main_chroms.txt`，每行一个染色体名，与 BED 第 1 列一致。聚类、图边和完整 BED 必须来自同一数据集。

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

生成 `screening.translocation.genes.tsv` 和 `screening.translocation.blocks.tsv`。
添加 `--allow-reverse` 可报告反向共线区块，`--include-unplaced` 可报告未定位基因对。

## `legacy-construct`：旧版 cDNA 构建流程

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

## `legacy-append`：追加旧版 cDNA 图

提供历史结果中未重命名的 GFF3/cDNA 配对及图包。另行准备 `combined_pre_and_new.bed`，包含历史 pre 聚类的代表基因和全部新增基因，gene ID 和品种标签须与序列输入一致。

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

若需预测 CDS/UTR，另行安装 TransDecoder，并在任一旧版命令中添加 `--annotate-cds --transdecoder-dir /path/to/TransDecoder`。
使用 `pantrans <子命令> --help` 查看该命令的全部参数。
