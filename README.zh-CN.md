# PanTrans

[English](README.md)

PanTrans 是一个命令行工具，用于从多个品种构建多倍体生物的泛基因组和泛转录组。

当前实现根据 cDNA-to-gDNA 比对结果构建基因关系图，生成 `pre` 和 `last` 两级基因簇，对过滤后的 BAM 比对结果执行转录本去冗余，并输出中间结果及最终参考序列集合。

## 功能

- 使用 `construct` 从多个品种构建泛基因簇和泛转录本簇
- 使用 `construct --bam` 复用已有的 cDNA-to-gDNA BAM，跳过 `minimap2` 比对
- 使用 `construct --main-chroms` 指定主染色体
- 从 `construct` 结果保存可复用的完整有向图包
- 生成 `pre` 和 `last` 两级结果，包括 `pre.tmp.cluster`、`last.tmp.cluster`、`pre.gtf`、`gtf`、参考 FASTA 和参考 BED
- 使用 `append` 在已有参考集合中追加新品种
- 在最终 GTF/cDNA 结果中同时保留原始代表 ID 和官方 `Pan...` ID

## 安装

在项目根目录以可编辑模式安装：

```bash
git clone https://github.com/changQing111/PanTrans.git
cd PanTrans
pip install -e .
```

如果环境无法联网，但已经安装了构建所需工具，可以使用：

```bash
pip install -e . --no-build-isolation
```

## 依赖

PanTrans 依赖以下软件：

- Python 3.8+
- `pysam`
- `networkx`
- `biopython`
- `minimap2`
- `samtools`

## 输入文件

PanTrans 使用以下类型的输入：

- `cdna` FASTA：转录本序列
- `gdna` FASTA：作为比对目标的基因组基因序列
- `bed` 文件：基因坐标和链方向
- 品种名称：可以是文本文件、逗号分隔字符串或空格分隔字符串

BED 文件必须在第 4 列提供基因 ID，在第 1 列提供染色体或 contig 名称。

## Construct

使用 `construct` 从多个品种构建新的泛基因组和泛转录组。

### 命令

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

### 可选参数

- `--bam`：使用已有的 cDNA-to-gDNA BAM，跳过 `minimap2` 比对
- `--main-chroms`：包含主染色体名称的文本文件，每行一个名称
- `--threads`：`minimap2` 使用的线程数

### 示例

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

### Construct 输出

当使用 `--prefix test_pantrans` 时，PanTrans 会生成：

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
- 比对 BAM 文件

### 输出命名说明

- `pre.gtf` 保留簇中心基因 ID 和转录本 ID。
- `<prefix>_unrenamed.gtf` 与 `<prefix>_unrenamed_cdna.refer.fasta` 保留图、BED 和簇文件使用的最终 `last` 代表 ID。
- `<prefix>.gtf` 与 `<prefix>_cdna.refer.fasta` 包含相同的最终转录本模型，但会将基因和转录本重命名为官方 `Pan...` ID，例如 `Pan1A000001`。

### 从已有 Construct 结果导出图包

对于在引入图包功能之前完成的 construct 运行，不需要重新执行历史 `minimap2` 比对。可以使用已有的 PanTrans 过滤 BAM 和原始完整输入直接导出图包：

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

BAM 必须是生成历史簇时使用的过滤 BAM。其参考序列名称必须覆盖完整历史 BED 中的每个基因。三个过滤阈值和 `--filter-logic-id` 必须与生成该过滤 BAM 时使用的逻辑一致，因为 BAM 头不会保存 PanTrans 过滤参数。

导出脚本会流式读取 BAM，并在图包清单中记录源 BAM 的程序头信息；它不会重新执行比对或聚类。

## Append

使用 `append` 复用历史 construct 图包，只计算新品种到历史集合的新增比对块。合并后的图会沿用 construct 的聚类和转录本处理流程。

### 命令

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

三个转录本输入的作用不同：

- `--cdna`：上一次运行生成的 `<prefix>_unrenamed_cdna.refer.fasta`
- `--history-gtf`：与该 cDNA 配对的 `<prefix>_unrenamed.gtf`
- `--query-cdna`：新品种的完整 cDNA FASTA

历史 GTF 和 cDNA 必须包含完全一致的转录本 ID。GTF 中的基因 ID 必须是历史图包中存在的未重命名 ID；如果传入官方 `Pan...` GTF，程序会拒绝。`--bed` 是一个合并 BED，包含历史 `pre.tmp.cluster` 中每个簇的第一个基因，以及新品种的全部基因。图包提供完整的历史 gDNA、BED、过滤 BAM、图边、节点长度、品种顺序、参考品种和比对来源信息。

为了支持断点续跑，`append` 还接受 `--query-to-all-bam` 和 `--history-to-query-bam`，用于复用已经完成的两类交叉比对 BAM。

### 示例

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

### Append 流程

`append` 会执行以下步骤：

1. 校验历史图包及其源文件身份
2. 合并历史非冗余 cDNA 和查询 cDNA，并构建完整的历史加查询 gDNA 与 BED
3. 将新品种 cDNA 比对到全部历史加查询 gDNA
4. 将历史非冗余 cDNA 比对到新品种 gDNA
5. 复用历史图边，合并两类交叉比对结果及其过滤 BAM
6. 完全沿用 construct 的聚类流程重新生成 `pre` 和 `last` 簇
7. 对最终转录本去冗余使用历史 GTF 模型作为种子，但只使用其原始代表基因仍是当前 `last` 代表的模型
8. 输出官方及未重命名的最终 GTF/cDNA，以及供下一次 append 使用的图包

所有在序列中出现的历史基因和新基因都可以参与聚类。append 的 BED 是面向用户的“历史代表基因加查询基因”输入；完整历史 BED 和 gDNA 集合由历史图包提供，并在内部使用。只有当历史 GTF 模型的原始代表基因仍然是当前 `last` 簇代表时，才会复用其坐标；程序不会把转录本模型投影到另一个代表基因的坐标系上。

两类交叉比对避免了重新计算历史 cDNA 对历史 gDNA 的比对。由于 `minimap2` 只保留有限数量的次级比对，增量图与一次性 construct 得到的图预计会高度相似，但不保证逐字节完全一致。

如需继续追加下一个品种，请将本次 append 生成的图包，以及对应的 `_unrenamed.gtf` 和 `_unrenamed_cdna.refer.fasta` 作为下一次运行的历史输入。

如果 append 在任意交叉比对完成后失败，可以使用对应的 `--query-to-all-bam` 和/或 `--history-to-query-bam` 路径重新运行。BAM 必须包含 `minimap2` 的 `@PG` 记录，且版本和生成图边时使用的参数一致；PanTrans 会使用当前阈值重新过滤复用的 BAM。

## 输出说明

- `pre.tmp.cluster`：图推导出的中间簇，在最终递归细分之前生成
- `last.tmp.cluster`：递归细分和染色体感知分配后的最终簇列表
- `pre.gtf`：在 `pre` 簇级别执行转录本去冗余的结果，使用原始簇中心基因 ID
- `_unrenamed.gtf`：最终 `last` 转录本模型，使用原始代表 ID
- `gtf`：在 `last` 簇级别执行转录本去冗余的结果，使用 `Pan...` 基因命名
- `*.refer.bed`：作为参考集合的簇中心基因 BED 条目
- `*_gdna.refer.fasta`：簇中心基因的基因组序列
- `*_cdna.refer.fasta`：根据生成的 GTF 重建的 cDNA 序列
- `*_unrenamed_cdna.refer.fasta`：以原始代表转录本 ID 为序列名的最终 cDNA
- `*.graph.json`：包含源文件身份和比对来源信息的图包清单
- `*.graph.edges.tsv`：去重后的有向图边
- `*.graph.nodes.tsv`：完整的图节点长度和 BED 元数据

## 当前行为与注意事项

- 主染色体检测基于 BED 第 1 列
- 未被识别为主染色体的名称会被视为 contig
- 对 `construct` 而言，`pre` 输出保留原始 ID，最终 `last` 输出使用 `Pan...` ID
- 对 `append` 而言，历史品种和新品种会作为独立品种参与分配，与 construct 保持一致；不会额外创建只包含代表基因的 `Refer` 块
- append 图包会保留源 `minimap2` 版本、参数和过滤阈值，用于在复用历史图包前识别不兼容情况
- 图包清单会校验源文件大小和修改时间。移动清单及其旁车文件是允许的，但修改已记录的 BED、BAM、cDNA 或 gDNA 源文件会使图包失效

## 开发说明

本项目仍在持续开发中。部分流程是逐步实现的，聚类、转录本恢复和 append 逻辑仍可能继续调整。
