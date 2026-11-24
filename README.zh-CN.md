# PanTrans

[English](README.md)

## 安装

在 Linux 或 WSL 中执行以下命令。Conda 安装 Python 和比对工具，pip 安装 PanTrans 及其 Python 依赖。

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

开发时可改用 `python -m pip install -e .`。当前 CDS 流程不需要 TransDecoder。

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

预处理同时输出全部基因的无侧翼和带侧翼 DNA、完整 cDNA、全部及有效 CDS，以及每基因一条的 `Pan_cds.longest.valid.fa` 和逐基因选择表 `Pan_cds.longest.valid.tsv`。最长候选仅在有效 CDS 中选择，等长时按输出转录本 ID 字典升序；无有效 CDS 或无转录本的基因保留 DNA/坐标并记录状态。construct 和 append 均使用最长有效 CDS。文件及 manifest 说明见[预处理输出](docs/preprocess-cds-inputs.md)。

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

construct 读取新版 preprocess 的最长有效 CDS，依次比对到每个品种的全部无侧翼基因体，包含品种内比对。构图参数为 `-f 20000 -N 2000`，其余为 `-ax splice:hq -uf -I 2G -K 50M --secondary=yes`，并设置独立 split-prefix。`--threads` 是每个任务的共享线程预算，目标任务顺序运行，避免额外倍增线程/内存。构图仍以全部基因为节点，注释仍使用全部合格转录本。旧输入缺少最长 CDS 时，请在新目录重新运行 preprocess。

`targets/0001/` 等目录保存各目标品种的 gene.fa、raw/filtered BAM、edges.tsv、饱和统计、命令和汇总，品种对应关系见 `graph.mapping.json`。最终新增 `last.gdna.fa`、`last.complete.cdna.fa`、完整转录本列表和 `validation.json`，逐条核对 DNA/exon 重建的 cDNA；无转录本模型的代表基因列入 `last.unannotated_representatives.tsv`，不补造模型。详见[新版构建流程](docs/construct-longest-cds.md)。

append 已与 construct 统一为最长有效 CDS、分目标品种及 `-f 20000 -N 2000`。历史→历史图边直接复用；新增 CDS 分别查询每个历史目标，全部 CDS 分别查询每个新增目标。旧 N100 图包必须用新版 construct 重新构建，不能与新策略混用。详见[增量流程](docs/cds-append.md)。

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

## `synteny`：两两品种共线性

基于同一数据集的完整预处理输入、pre簇及双向同源图边，生成全部无序品种对的正向/反向共线区块。

```bash
pantrans synteny \
  --prepared prepared --results construct_results \
  --max-copy 3 --distance 20 --min-anchors 5 --lift-distance 10 \
  --output synteny_results
```

append结果使用其`merged_prepared/`作为`--prepared`。默认运行全部步骤；大数据可用
`--stage index/shard/pair/collect`分阶段调度，pair用`--pair-index`选择。只读现有图，不重新比对，
不改变pre/last。输出区块、锚点、位置支持、完整二进制候选及汇总。见[共线性说明](docs/synteny.md)。
原`translocation`命令已移除；旧分析脚本仍可读取已有易位报告。

当前仅保留 `preprocess`、`construct`、`append`、`cnv` 和 `synteny` 五个子命令。
`legacy-construct` 和 `legacy-append` 已移除。使用 `pantrans <子命令> --help` 查看参数。

## 致谢

OpenAI Codex 为代码开发、测试和文档编写提供了 AI 辅助支持。项目作者负责审查相关贡献，并对发布代码及其科学应用承担责任。
