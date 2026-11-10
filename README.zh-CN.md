# PanTrans

[English](README.md)

**0.8.0：新版CDS构图与UTR-aware注释已迁入正式 `pantrans construct/append`。**
实现位于可安装的Python包中，不依赖仓库的scripts目录；原实验脚本保留为兼容入口。

## 安装与依赖

```bash
pip install .
pantrans --version
```

需要Python 3.8+、pysam、NetworkX、Biopython，以及PATH中的minimap2和samtools。
CDS流程不需要TransDecoder；它仅用于旧流程的可选预测。

## 预处理与构建

```bash
pantrans preprocess \
  --genome-list genomes.txt --annotation-list annotations.txt \
  --variety-list varieties.txt --flank 3000 \
  --prefix sample --output prepared

pantrans construct --prepared prepared --output construct_results \
  --threads 32 --workers 8
```

三个列表按相同顺序每行填写一个基因组FASTA、GFF3注释路径和品种名。
基因组使用未压缩FASTA。preprocess生成schema-2 manifest、无侧翼gene body、
带侧翼gDNA、完整cDNA、有效CDS及坐标/模型表。gene ID须全局唯一；若使用
`--namespace-ids`，历史与新增批次应一致。详见[输入说明](docs/preprocess-cds-inputs.md)。

构图使用**所有有效CDS比对所有无侧翼gene body**，相似度、query CDS覆盖度、
target CDS覆盖度均≥90%，取消编码区剪接一致性和移码indel否决。
覆盖度仍按对齐碱基与目标转录本已注释CDS的重叠计算。
未比对、反向、supplementary及缺少NM的记录仍排除。

SCC生成pre后，按修正的图支持分配生成last。目前策略固定为参考品种CS、
小麦Chr1A–Chr7D主染色体（支持别名）；A/B/D强制划分只是独立实验，不是默认。
详见[分配规则](docs/last-assignment-graph-support.md)。

默认**仅注释last**。`--annotate-pre`可额外生成pre注释。完整cDNA比对每个hub的
带侧翼gDNA，与hub原始模型合并后进行UTR-aware去冗余。唯一兼容的partial/CDS-only
证据可并入完整结构，多重归属和无匹配模型保留明确标记。
**注释转移仍有严格编码区结构和ORF检查**；构图允许建立同源边，不代表所有异构体都能
转移成功。详见[同源图v2](docs/cds-homology-v2.md)。

## 追加新品种

```bash
pantrans append \
  --history-results construct_results --history-prepared prepared \
  --query-prepared new_prepared --output append_results \
  --threads 32 --workers 8
```

复用历史全部图边，执行新增CDS→全部gene、历史CDS→新增gene两轮比对，再重新聚类和
生成last注释。历史与新增输入均须schema 2，品种名及gene/transcript ID不重叠。
需要历史原始完整输入，不能仅提供代表序列。

连续追加时：

```bash
pantrans append \
  --history-results append_results \
  --history-prepared append_results/merged_prepared \
  --query-prepared next_prepared --output next_append_results
```

`cds.graph.json`校验输入/图边哈希和比对过滤策略。旧cDNA图包、CDS v1图包不能混用。
没有图包的历史CDS raw BAM结果可显式使用`--import-legacy-history`导入，无法追溯的
历史坐标/比对策略会记录为假定。[追加说明](docs/cds-append.md)。

由于minimap2候选输出受target集合影响，**append不保证与直接construct完全相同**。

## 输出和运行

| 文件 | 含义 |
|---|---|
| `pre.cluster`、`last.cluster` | 每行一个簇，首个gene为代表 |
| `graph.edges.tsv`、`graph.summary.json` | 有向同源图及统计 |
| `cds.graph.json` | CDS图包与输入/策略校验信息 |
| `last.gff3`、`last.cdna.fa` | 最终注释和重建cDNA |
| `last.annotation.bam` | 去冗余前接受的转移比对证据 |
| `last.model_sources.tsv`、`last.alignment_status.tsv` | 来源、合并原因和转移状态 |
| `last.evidence.sqlite`、`last.cache_stats.json` | 转移证据缓存及命中统计 |
| `last.summary.json`、`summary.json` | 数量和完成统计 |
| `progress.json`、`timings.json` | 注释进度和分阶段计时 |

GFF3为**hub对应的基因定向、带侧翼DNA的局部坐标**，不是原始染色体坐标。
转录本按`hub.N`重新编号，追溯来源请查看来源表。partial模型显式保留。
`body.fa`和`cds.fa`链接到prepared输入，需保留源文件。每次运行使用新输出目录。

`--threads`控制构图比对，`--workers`控制注释并发。append自动使用历史结果中的
`last.evidence.sqlite`；`--history-evidence-cache`指定其他缓存，`--no-evidence-cache`
禁用缓存。所有当前候选仍重新去冗余。`--local-fasta`为可选节点本地参考拷贝，默认关闭。
详见[性能与验证](docs/cds-performance.md)。

`construct --reuse-graph PATH`只重跑注释，保留历史图来源；无图包的旧结果会标记为
策略未验证、阈值未知。`--genes`子集运行及仅注释重跑不会发布可供append使用的完整新图包。

## 旧命令迁移

旧的`--cdna/--gdna/--bed/--bam/--annotate-cds`参数属于
`pantrans legacy-construct`及`pantrans legacy-append`。
不会将旧参数静默解释为新版CDS输入。
[迁移说明](docs/cli-migration.md) · [旧流程参考](LEGACY.zh-CN.md)。
`preprocess`、`cnv`、`translocation`仍可使用，绘图脚本继续独立提供。

```bash
pantrans construct --help
pantrans append --help
pantrans legacy-construct --help
```
