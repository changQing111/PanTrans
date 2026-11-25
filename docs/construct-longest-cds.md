# Construct：最长有效CDS与分品种目标构图

2026-11-23。实现于 `cds_construct.py`、`cds_target_graph.py` 和 `cds_outputs.py`。

```bash
pantrans preprocess --genome-list genomes.txt --annotation-list annotations.txt \
  --variety-list varieties.txt --flank 3000 --prefix Pan --output prepared
pantrans construct --prepared prepared --output construct_results --threads 32 --workers 8
```

## 输入与查询

新构图要求schema-2 manifest含`files.cds_longest`及`files.cds_longest_selection`。
首先按有效模型重算选择（CDS长度降序，输出转录本ID字典升序），校验逐基因表完整性，
再逐条核对最长CDS与全部有效CDS中的对应序列及模型长度。缺少新文件的旧输入明确报错，
要求在新目录重新预处理，不静默退回旧构图策略。`--genes`子集只过滤图输入，且不发布全数据图包。

`longest.cds.fa`为校验后的查询快照，`query_metadata.tsv`记录来源。原始`cds.fa`仍保留
全部有效CDS链接，以支持注释及旧接口；实际构图查询路径由每目标的command.json记录。
全部转录本模型和完整cDNA仍参与注释，最长CDS选择不删除其他异构体。

## 分目标品种比对与过滤

按输入基因表中的品种顺序，将完整gene-body DNA分流至`targets/0001/gene.fa`等目录，
数字目录避免将品种标签直接解释为路径。每个目标品种查询全部最长CDS，覆盖品种内关系。
命令为：

```text
minimap2 -ax splice:hq -uf -I 2G -K 50M --secondary=yes -f 20000 -N 2000 \
  --split-prefix <unique-temp>/mm2 -t <mapping-threads> <target>/gene.fa longest.cds.fa
```

使用原子BAM发布和独立临时split-prefix。线程按现有共享预算分给minimap2及samtools；
32线程对应minimap2 28、samtools `-@ 3`。目标任务默认依次运行，不并行倍增总资源。
`--workers`仍只控制后续转录本注释。

过滤复用construct/append的homology-v2函数：仅正向、非supplementary、有NM标签，
PID及查询/目标CDS重叠覆盖度均≥0.90；任一有效目标异构体满足即可，不施加编码结构否决。
校验每目标BAM的完整参考名称/长度、查询primary集合（含unmapped）和唯一性。
每目标保留raw.bam、filtered.bam、edges.tsv、command.json、summary.json、saturated_CDS.tsv。
饱和数统计原始mapped且非supplementary记录中的secondary，报告≥500/1000/2000；
当前上限为2000。达到上限的基因只记录，不删除。低于上限不证明穷尽比对。

`graph.mapping.json`记录品种与目录、命令策略、进度、各目标统计和唯一饱和基因数；
`query_target_variety_counts.tsv`给出所有有向品种组合，包括无有效查询的品种的零行。

## 聚类与注释

全基因节点及合并后的有向边进入现有SCC、pre递归分解和last分配，不改动生物学聚类规则。
图边分节点排序写出，避免同时复制整个边集合排序。检查pre/last完整无重复分簇、last直接
双向hub支持及主染色体同品种最多一个成员，写入`graph.validation.json`。

默认仅注释last；`--annotate-pre`仍可显式生成pre注释。保留原始hub模型种子，全部合格成员
cDNA使用单hub、N100的严格编码投影和ORF规则；UTR去冗余及证据缓存规则保持不变。

输出增加`last.gdna.fa`、`last.complete.cdna.fa`和`.fai`（非空时）、
`last.complete.transcripts.txt`、`last.unannotated_representatives.tsv`、`validation.json`。
逐条检查GFF边界/ID/Parent并以导出DNA和exon复原全部cDNA；只选择`completeness=complete`
进入完整cDNA。没有原始转录本的基因仍在图中；若最终代表无模型，则明确列入未注释表，
不补造cDNA。GFF坐标是基因方向、带侧翼DNA的局部坐标。验证成功后才写summary.json。

## 图包和复用

新全量构图发布`cds.graph.json`，分别记录过滤策略及最长CDS／分品种／N2000 mapping策略，
绑定原始完整输入、新最长CDS及选择表、合并图边、目标映射清单的SHA256。
`--reuse-graph`识别新图包并验证，复用每目标BAM/边，不要求不存在的单一cds.raw.bam。
旧N100图的注释重启仍支持；无图包的复用只检查成员划分及源文件身份，不声称策略验证。
重启目录可继续重启，并通过reused_graph.json引用源证据，不产生新的可追加图包。

append已与construct共享最长有效CDS、分品种目标及N2000参数。历史目标仅查询新增CDS，
新增目标查询全部CDS；历史→历史边直接复用。旧N100图包需重构，不能通过伪改manifest或
legacy导入绕过策略校验。详见[增量流程](cds-append.md)。
