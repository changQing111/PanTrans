# CDS append：最长有效CDS和分品种目标的增量构建

当前append与construct使用同一构图策略：每基因最长有效CDS，按目标品种分别比对，
参数为 `-ax splice:hq -uf -I 2G -K 50M --secondary=yes -f 20000 -N 2000`，
独立临时split-prefix。注释阶段仍为单hub的N100比对。

## 输入

```bash
pantrans append --history-results history_results \
  --history-prepared history_prepared --query-prepared new_prepared \
  --output append_results --threads 32 --workers 8
```

两批输入均须为新版schema-2预处理，含最长有效CDS及逐基因选择表，并保留全部原始
基因和转录本。品种名、基因ID和转录本ID不得跨批次重叠；预处理策略须一致。
历史结果要求cds.graph.json及最终summary.json。校验完整输入、新最长CDS文件/选择表、
历史图边和目标映射清单SHA256，并要求过滤策略、比对策略及minimap2版本兼容。

## 比对关系

设历史H有h个品种，新增N有n个品种：

| 查询→目标 | 处理方式 |
|---|---|
| H→H | 直接复用历史图边，包含历史品种内自比对 |
| N→每个历史品种 | 新增最长有效CDS合并集查询该品种全部gene-body DNA |
| H∪N→每个新增品种 | 全部最长有效CDS合并集查询该品种全部gene-body DNA |

新增有向品种组合数为2hn+n²，复用h²，总计(h+n)²。实际执行h+n个目标任务，
不是每个品种对单独启动任务。所有新增品种内部及新增品种之间的关系均被覆盖。
例如8品种扩到18品种，复用64个组合，新增260个组合。
任务按目标顺序执行，每任务共享threads预算；32线程对应minimap2 28、samtools -@ 3。
无有效新增CDS时，历史目标任务记录skipped_empty_cds，不生成空查询BAM；
新增目标仍由历史CDS查询。没有原始转录本的新基因也保留节点/DNA/坐标，不补造模型。

构图沿用homology-v2：PID和两侧CDS重叠覆盖度均≥0.90，任一有效目标异构体满足即可，
不施加编码结构否决。只补新关系后，在完整合并图上重新计算SCC、pre和last，不固定历史hub。

## 注释及连续追加

默认只注释last。全部合格转录本用于严格编码结构、覆盖和ORF投影。
历史last.evidence.sqlite通过查询/目标序列与模型、转移代码、minimap2二进制及依赖版本
指纹复用去冗余前证据，包括未获得合格模型的状态。原始hub种子重新生成。
无论证据命中与否，全部当前候选重新进行UTR去冗余；不复用历史最终GFF3行。
缓存来源只读，新缓存仅在注释成功后发布。可用--history-evidence-cache添加缓存、
--no-evidence-cache禁用缓存、--local-fasta选择本地参考暂存。

输出包括merged_prepared/、longest.cds.fa（全部查询）、new.longest.cds.fa（新增查询）、
graph.edges.tsv、pre.cluster、last.cluster、graph.validation.json、last.*和validation.json。
后者包括代表DNA、complete cDNA和逐条exon/cDNA重建校验，无模型的代表显式列出。

每目标targets/0001/等保存gene.fa、raw/filtered BAM、edges.tsv、command.json、summary.json、
saturated_CDS.tsv；无查询任务无BAM。历史目标edges.tsv只含本轮新来源边，复用边在父代保留。
graph.mapping.json记录query_scope（new/all）、历史边数、本轮/复用组合及任务统计。
query_target_variety_counts.tsv仅含本轮2hn+n²组合。饱和统计也仅描述本轮比对；
历史饱和及比对计数留在父代来源，不伪装成重新测量。

新图包cds.graph.json绑定完整输入、合并边、目标映射清单哈希和父代来源，构图完成后发布；
最终summary.json在注释和输出验证全部成功后发布。下一轮使用本轮结果作为history-results，
merged_prepared/作为history-prepared。所有历史输入只读；保留各代目录/BAM，不生成跨代单一raw BAM。
construct --reuse-graph可针对append结果重启注释，并识别空查询跳过的任务。

## 兼容性及验证边界

旧N100全部有效CDS图包与新版不兼容，须重新construct。
--import-legacy-history保留参数以给出明确拒绝信息，不导入无法证明策略的旧BAM。
旧cDNA图不受当前append支持，legacy-construct和legacy-append命令均已移除。
本次不自动改写既有8/18品种测试结果。

真实minimap2/SAMtools小数据覆盖连续追加、多历史/新增品种、所有方向/自比对、变更hub、
缓存复用、与直接construct比较图/簇/GFF/cDNA、空CDS、无转录本基因、损坏输入/图和失败处理。
相同目标/参数消除了旧合并目标流程中的目标竞争差异，但不承诺任意版本/环境及大数据
逐字节等价；不会从直接construct比较结果中替换增量图边。
