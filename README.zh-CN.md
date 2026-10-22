# PanTrans

[English](README.md)

PanTrans 是一个命令行工具，用于从多个品种构建多倍体生物的泛基因组和泛转录组。

当前实现根据 cDNA-to-gDNA 比对结果构建基因关系图，生成 `pre` 和 `last` 两级基因簇，对过滤后的 BAM 比对结果执行转录本去冗余，并输出中间结果及最终参考序列集合。

## 功能

**0.7.0** 将所有新生成的注释统一为 GFF3，并支持 TransDecoder CDS/UTR 预测。

- 用 `cnv` 筛选邻近基因拷贝，用 `translocation` 从 pre cluster 筛选候选染色体间位置差异

- 使用 `construct` 从多个品种构建泛基因簇和泛转录本簇
- 使用 `construct --bam` 复用已有的 cDNA-to-gDNA BAM，跳过 `minimap2` 比对
- 使用 `construct --main-chroms` 指定主染色体
- 从 `construct` 结果保存可复用的完整有向图包
- 生成 `pre` 和 `last` 两级结果，包括 `pre.tmp.cluster`、`last.tmp.cluster`、`pre.gff3`、`gff3`、参考 FASTA 和参考 BED
- 使用 `append` 在已有参考集合中追加新品种
- 在最终 GFF3/cDNA 结果中同时保留原始代表 ID 和官方 `Pan...` ID

## 易位结果可视化

独立脚本 `scripts/plot_translocations.py` 提供染色体总览和局部基因匹配图，
使用三次贝塞尔曲线连接，支持 `--flank`、可选歧义匹配，以及 PDF/PNG/SVG 输出。
绘图依赖 Matplotlib，需单独安装：

```bash
pip install matplotlib
python scripts/plot_translocations.py --help
```

参数和输入格式见[绘图说明](scripts/plot_translocations.md)。脚本不重新检测易位；
图中范围为候选锚点覆盖区间，不是断点。

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

### 可选 CDS / UTR 注释

所有新生成的注释（包括 pre、最终两套结果及 rescue 中间结果）均使用 `.gff3`。
文件包含 `##gff-version 3`，使用 `gene → mRNA → exon/CDS/UTR` 层级及
`ID`/`Parent` 属性；特殊字符按 GFF3 规则转义。正式 ID 版本重命名各层级 ID
及其 Parent，但第 1 列仍为原始代表基因序列名，坐标仍对应配套 gDNA。
append 使用 `--history-gff3`；该参数也能读取历史 GTF，旧参数 `--history-gtf`
保留为兼容别名。旧结果文件不会自动转换或删除，新结果不再生成 GTF。
`scripts/repair_official_ids.py` 现针对新的 `.gff3` 输出对进行修复。

`construct` 和 `append` 支持 `--annotate-cds`，在最终转录本去冗余和 rescue
完成后，用 TransDecoder 预测主要 ORF，并同时为 `<prefix>_unrenamed.gff3`
和 `<prefix>.gff3` 添加 `CDS`、`five_prime_UTR`、`three_prime_UTR`，以及可确定的
`start_codon` / `stop_codon`。不对 pre 结果运行预测。

在原有命令后添加：

```bash
--annotate-cds --transdecoder-dir /path/to/TransDecoder --min-orf-length 50
```

- TransDecoder 在 `PATH` 中时可省略 `--transdecoder-dir`；单独指定该目录也会启用注释。
- 支持经典 `TransDecoder.LongOrfs` + `TransDecoder.Predict` 和新版 `TransDecoder` 入口；需要可运行的 Linux/WSL 工具环境及其依赖。
- 默认关闭注释；`--min-orf-length` 默认 50，单位为氨基酸，是候选 ORF 长度阈值。
- TransDecoder 默认串行评分，不受 `-t` 影响；增加 `--transdecoder-parallel -t 8` 可使用最多 8 个评分进程，`-t 1` 保留原版串行流程。
- 并行评分目前支持评分调用结构匹配的经典 TransDecoder 5.7.1；其他版本提示警告并回退串行。各分块共享同一训练模型，调用原版评分脚本，按输入顺序合并结果。LongOrfs、模型训练和后续预测仍串行执行。
- 系统安装的 TransDecoder 不会被修改；本次运行目录保留专用 Predict 副本和各分块诊断日志。

- 对最终未重命名 cDNA 仅运行一次，使用正义链预测和 `--single_best_only`，然后沿 exon 映射回原有代表基因序列坐标。正式 GFF3 由同一结果重命名得到。
- GFF3 的 CDS 包含终止密码子，与 TransDecoder 一致；3′UTR 从终止密码子之后开始。跨外显子的 CDS 拆段、共享 ID 并计算 phase。可确定的起止密码子另加重叠特征，其 phase 为 `.`。
- 不完整 ORF 只在已确定边界的一端添加 UTR；没有预测到 ORF 的转录本保留 transcript/exon，不直接判为非编码。
- 每次运行保留 `<prefix>_transdecoder/run-*/` 中的日志及 TransDecoder 原始结果（含 GFF3、蛋白和 CDS FASTA，使用原始转录本/ORF ID），避免复用过期检查点。工具失败或结果不合法时明确报错。
- append 在启用此选项时重新预测最终转录本；历史 GFF3 的 CDS/UTR 不直接继承。继续需要这些注释时，每次 append 均启用该选项。

局部 GFF3 的转录本链方向来自所选 BAM 比对相对于输入 gDNA 的方向，而非染色体
BED 第 6 列。cDNA 按 GFF3 的 `+` 正向拼接 exon，按 `-` 在拼接后反向互补；
未知方向会报错。已按基因方向调整的 gDNA 不会因原始 BED 为负链而被重复反向互补。
相反方向的模型分别去冗余；若一个基因包含两种方向，其 gene 行使用 `.`，mRNA
及子特征保留各自方向。旧版输出的错误方向不会自动修复，需复用 BAM 重新生成
注释和 cDNA，再重跑 TransDecoder；用于 append 的历史注释也应先修复。

这些区段是基于 ORF 的预测，不是实验确认的翻译边界；当前使用标准遗传密码。

PanTrans 使用以下类型的输入：

- `cdna` FASTA：转录本序列
- `gdna` FASTA：作为比对目标的基因组基因序列
- `bed` 文件：基因坐标和链方向
- BED 第 5 列：每个基因对应的品种名
- 品种名称：从 BED 第 5 列自动推断

BED 文件必须在第 4 列提供基因 ID、第 5 列提供品种名、第 1 列提供染色体或 contig 名称；第 6 列可提供链方向，省略时按 `.` 处理。

`construct` 和 `append` 直接以 BED 第 5 列的标签定义每个基因的品种，品种列表也从该列自动生成，不再接受单独的品种名参数。

## Construct

使用 `construct` 从多个品种构建新的泛基因组和泛转录组。

### 命令

```bash
pantrans construct \
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
- `--threads`：minimap2 与 samtools 共享的线程预算
- `--rescue`：显式启用对最终未生成 GFF3 模型的簇代表基因进行恢复；默认关闭

### 示例

```bash
pantrans construct \
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
- 比对 BAM 文件

### 输出命名说明

- `pre.gff3` 保留簇中心基因 ID 和转录本 ID。
- `<prefix>_unrenamed.gff3` 与 `<prefix>_unrenamed_cdna.refer.fasta` 保留图、BED 和簇文件使用的最终 `last` 代表 ID。
- `<prefix>.gff3` 与 `<prefix>_cdna.refer.fasta` 包含相同的最终转录本模型，但会将基因和转录本重命名为官方 `Pan...` ID，例如 `PanChr1A000001`。`Pan` 之后原样保留代表基因对应的 BED 染色体名称。

### 从已有 Construct 结果导出图包

对于在引入图包功能之前完成的 construct 运行，不需要重新执行历史 `minimap2` 比对。可以使用已有的 PanTrans 过滤 BAM 和原始完整输入直接导出图包：

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

BAM 必须是生成历史簇时使用的过滤 BAM。其参考序列名称必须覆盖完整历史 BED 中的每个基因。三个过滤阈值和 `--filter-logic-id` 必须与生成该过滤 BAM 时使用的逻辑一致，因为 BAM 头不会保存 PanTrans 过滤参数。

导出脚本会流式读取 BAM，并在图包清单中记录源 BAM 的程序头信息；它不会重新执行比对或聚类。

## Append

使用 `append` 复用历史 construct 图包，只计算新品种到历史集合的新增比对块。合并后的图会沿用 construct 的聚类和转录本处理流程。

### 命令

```bash
pantrans append \
  --cdna <previous_unrenamed_final_cdna.fasta> \
  --history-gff3 <previous_unrenamed_final.gff3> \
  --query-cdna <new_variety_cdna.fasta> \
  --gdna <query_gdna.fasta> \
  --bed <combined_representatives_and_query.bed> \
  --history-graph <historical_prefix.graph.json> \
  --prefix <prefix> \
  --output <output_dir>
```

三个转录本输入的作用不同：

- `--cdna`：上一次运行生成的 `<prefix>_unrenamed_cdna.refer.fasta`
- `--history-gff3`：与该 cDNA 配对的 `<prefix>_unrenamed.gff3`
- `--query-cdna`：新品种的完整 cDNA FASTA

历史 GFF3 和 cDNA 必须包含完全一致的转录本 ID。GFF3 中的基因 ID 必须是历史图包中存在的未重命名 ID；如果传入官方 `Pan...` GFF3，程序会拒绝。`--bed` 是一个合并 BED，包含历史 `pre.tmp.cluster` 中每个簇的第一个基因，以及新品种的全部基因。图包提供完整的历史 gDNA、BED、过滤 BAM、图边、节点长度、品种顺序、参考品种和比对来源信息。

为了支持断点续跑，`append` 还接受 `--query-to-all-bam` 和 `--history-to-query-bam`，用于复用已经完成的两类交叉比对 BAM。只有显式指定 `--rescue` 才会启用缺失 GFF3 模型的恢复，默认不执行恢复。

### 示例

```bash
pantrans append \
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

### Append 流程

`append` 会执行以下步骤：

1. 校验历史图包及其源文件身份
2. 合并历史非冗余 cDNA 和查询 cDNA，并构建完整的历史加查询 gDNA 与 BED
3. 将新品种 cDNA 比对到全部历史加查询 gDNA
4. 将历史非冗余 cDNA 比对到新品种 gDNA
5. 复用历史图边，合并两类交叉比对结果及其过滤 BAM
6. 完全沿用 construct 的聚类流程重新生成 `pre` 和 `last` 簇
7. 对最终转录本去冗余使用历史 GFF3 模型作为种子，但只使用其原始代表基因仍是当前 `last` 代表的模型
8. 输出官方及未重命名的最终 GFF3/cDNA，以及供下一次 append 使用的图包

所有在序列中出现的历史基因和新基因都可以参与聚类。append 的 BED 是面向用户的“历史代表基因加查询基因”输入；完整历史 BED 和 gDNA 集合由历史图包提供，并在内部使用。只有当历史 GFF3 模型的原始代表基因仍然是当前 `last` 簇代表时，才会复用其坐标；程序不会把转录本模型投影到另一个代表基因的坐标系上。

两类交叉比对避免了重新计算历史 cDNA 对历史 gDNA 的比对。由于 `minimap2` 只保留有限数量的次级比对，增量图与一次性 construct 得到的图预计会高度相似，但不保证逐字节完全一致。

如需继续追加下一个品种，请将本次 append 生成的图包，以及对应的 `_unrenamed.gff3` 和 `_unrenamed_cdna.refer.fasta` 作为下一次运行的历史输入。

如果 append 在任意交叉比对完成后失败，可以使用对应的 `--query-to-all-bam` 和/或 `--history-to-query-bam` 路径重新运行。BAM 必须包含 `minimap2` 的 `@PG` 记录，且版本和生成图边时使用的参数一致。PanTrans 还会将 BAM 中的主查询 ID 和序列摘要与当前 cDNA FASTA 对照，然后使用当前阈值重新过滤复用的 BAM。

## CNV 筛选

使用已有 pre cluster 和完整 BED，不重新运行比对或聚类：

```bash
pantrans cnv \
  --pre-cluster five_variety_construct_out/five_variety_pre.tmp.cluster \
  --bed five_variety_genes.bed \
  --reference CS \
  --output five_variety.cnv.genes.txt
```

逐簇检查是否存在同一品种（BED 第 5 列）、同一染色体上的至少两个不同基因，
BED 区间间隔不超过 1,000,000 bp。基因按坐标排序后，间隔为
`max(0, 下一基因 start - 当前组最大 end)`；重叠间隔为 0，不限制链方向。
采用链式邻近：A–B 和 B–C 均符合条件时，三个基因归为一组，即使 A–C 超过阈值。
每个簇发现首个符合条件的组后即停止检测，不要求不同品种的拷贝数不同。

输出为**每个命中 pre cluster 一个原始 gene ID**，无表头，按输入簇顺序排列。
代表从**整个簇**选择：优先 `-r` 指定品种，其次 BED 区间长度 `end-start` 最大，
再次按 gene ID 字典序选择。代表不一定属于首个检出的邻近组。
该簇不含参考品种时选择其他成员中最长者；参考品种必须存在于完整 BED 中。
无命中时输出空文件。`-o` 指定输出文件路径，`--max-gap` 可调整间隔阈值（单位 bp）。

BED 必须覆盖所有簇成员，使用 0-based、半开区间和有效品种标签。
不同未挂载序列必须使用不同的 BED 第 1 列名称；若统一标成 `contig`，无法据此区分
不同序列上的局部坐标。此命令实现上述邻近筛选规则，不进行品种间拷贝数差异比较。

## 候选易位检测

使用同一次构建的 pre cluster、完整 BED 和有向图边文件，将每个非参考品种分别与
`-r` 指定品种比较，不重新比对、加载 BAM 或聚类：

```bash
pantrans translocation \
  --pre-cluster five_variety_construct_out/five_variety_pre.tmp.cluster \
  --bed five_variety_genes.bed \
  --edges five_variety_construct_out/five_variety.graph.edges.tsv \
  --main-chroms main_chroms \
  --reference CS \
  --output-prefix five_variety
```

BED 第 5 列定义品种，必须包含所有簇成员和图边端点。一个基因只能出现在一个
pre cluster 中。图边 TSV 为无表头的两列 `source_gene`、`target_gene`。
只接受同一 pre cluster 内直接存在两个方向图边的参考/目标品种基因对，
自连接、单向连接和间接路径不算支持，重复边不增加证据。
不含参考品种基因的簇跳过，并在日志中统计。

`--main-chroms` 必填，一行一个 BED 染色体原名；各品种染色体名称应保持一致。
此命令不把 `Chr1A` 自动转换为 `1A`。不在列表中的序列视为未定位（Ctg/contig）。
默认排除任一端为未定位基因的关系，不输出，也不参与双向匹配数和拷贝歧义统计。
显式添加 `--include-unplaced` 后，才将这些关系纳入统计并以 `unplaced` 状态报告，
仍不用于构建易位区段。

同染色体双向匹配优先视为正常匹配，不写入候选表。对于跨染色体基因对，
如果任一端点还有同染色体双向匹配，或在当前参考/目标品种比较中有多个双向匹配，
则标记为 `ambiguous_copy`，保留供检查，但不作为区段锚点。
仅启用 `--include-unplaced` 时，额外的未定位双向匹配才计入匹配数量。
首版不会仅凭邻域证据消除多拷贝歧义。

只有一对一、直接双向连接的跨染色体基因对才能作为锚点。
区段内必须是同一目标品种、相同的参考染色体与目标染色体组合，
并且两侧起始坐标顺序一致；使用 `--allow-reverse` 时也接受目标侧整体反向。
区间可以重叠，但起始坐标相同不能提供独立顺序证据。

| 参数 | 默认值 | 含义 |
| --- | --- | --- |
| `--min-anchor-genes` | `3` | 支持区段至少来自多少个不同 pre cluster；最小允许值为 2 |
| `--max-gap` | `1000000` | 相邻锚点在两侧各自的 BED 区间间隔上限，单位 bp；重叠按 0 计算 |
| `--max-skipped-genes` | `5` | 两侧相邻锚点之间各自允许跳过的基因数，按完整 BED 计数，包含非锚点基因 |
| `--allow-reverse` | 关闭 | 接受目标基因顺序反向的区段，方向由基因排列决定，不依赖 BED 链方向 |
| `--include-unplaced` | 关闭 | 纳入未定位基因对的报告及歧义统计，但不作为区段锚点 |

算法在上述限制内计算最长路径候选，再确定性地选择互不重叠的锚点链；
一个锚点最多进入一个区段，不穷举所有可能的重排。
距离限制作用于相邻锚点，不限制区段总长度。
单个大基因家族/pre cluster 即使提供多个基因对，也不能独自构成支持区段。

输出两个 TSV 文件，`-o` / `--output-prefix` 指定的是**文件前缀**：

- `five_variety.translocation.genes.tsv`：原始基因 ID、两侧 BED 坐标、品种、
  双向匹配数、是否存在同染色体替代匹配、状态、原因、区段 ID。
  `pre_000123` 表示输入 cluster 文件的第 123 行。
  状态为 `supported_block`、`single_gene_candidate`、`ambiguous_copy` 或 `unplaced`。
- `five_variety.translocation.blocks.tsv`：候选区段、排列方向、锚点数、不同簇数及支持基因。
  至少两个锚点但独立簇数未达阈值的链标为 `low_support_block`，其基因行仍标为
  `single_gene_candidate`；单个孤立基因对没有区段，`block_id` 为 `.`。

所有坐标为 BED 的 0-based 半开区间。区段范围覆盖支持基因，**不是易位断点**。
无命中时也输出表头。两份报告先写入临时文件，发布失败时恢复旧报告；
日志记录阈值和统计信息。
`supported_block` 表示存在共线性证据的候选位置差异，不等于已证实易位，
不能确定变化发生在哪个品种的进化分支，也不能排除组装、注释或拷贝歧义。

## 输出说明

- `pre.tmp.cluster`：图推导出的中间簇，在最终递归细分之前生成
- `last.tmp.cluster`：递归细分和染色体感知分配后的最终簇列表
- `pre.gff3`：在 `pre` 簇级别执行转录本去冗余的结果，使用原始簇中心基因 ID
- `_unrenamed.gff3`：最终 `last` 转录本模型，使用原始代表 ID
- `gff3`：在 `last` 簇级别执行转录本去冗余的结果，使用 `Pan...` 基因命名
- `*.refer.bed`：作为参考集合的簇中心基因 BED 条目
- `*_gdna.refer.fasta`：簇中心基因的基因组序列
- `*_cdna.refer.fasta`：根据生成的 GFF3 重建的 cDNA 序列
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

对支持的 TransDecoder 5.7.1，串行和并行均使用运行目录中的确定性 ORF 选择脚本：保持原优先级，完全并列时选择字典序最小的 ORF ID。原先随机选中的并列 ORF 因而可能改变；不支持的安装会警告无法保证并列选择确定性。
