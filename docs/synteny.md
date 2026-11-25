# PanTrans synteny

`pantrans synteny`替代已移除的translocation命令，分析全部无序品种对的基因顺序共线性。
共线性证据生成与pre→last分配是不同步骤；本命令不修改construct/append结果或采用试验性的last优化规则。

## 输入和执行

```bash
pantrans synteny --prepared prepared --results construct_results -o synteny_results
# 显式路径，适用于没有标准results目录的已验证历史输出
pantrans synteny --prepared prepared --pre-cluster results/pre.cluster \
  --edges results/graph.edges.tsv -o synteny_results
```

--prepared要求schema-2、完整BED（品种位于第5列），基因与pre成员必须完整一致且无重复。
--results定位pre.cluster和完整graph.edges.tsv；append的prepared须使用merged_prepared/。
候选必须同pre且具有直接双向图边。只扫描一次完整图后按品种方向分片，不依赖每目标目录结构，
因此适用于construct和多代append。使用NumPy压缩候选和索引，不需要JCVI、编译器或C++运行内核。
当前实现为Python有界邻域扫描，大型高密度候选的运行时间需实际评估，不宣称原C++原型性能。

默认--max-copy 3、--distance 20、--min-anchors 5、--lift-distance 10。
--main-chroms指定主染色体名列表，默认小麦Chr1A–Chr7D。基因rank在各品种/染色体内按
坐标和ID排序，使用全部基因，而非只用锚点。双方pre拷贝数均不超过max-copy且均位于主染色体
的候选可成为种子；高拷贝和未定位候选不参与种子构建。

在每个染色体对中，先按双方rank距离均≤distance形成种子连通分量，再反复取最长严格单调链。
正向和反向都检查，相邻锚点在双方rank的步长均须>0且≤distance；每条链至少min-anchors对。
链的坐标范围为BED原始0-based half-open区间。非锚点主染色体候选可由既有区块提供位置支持：
双方距离≤lift-distance，满足区块方向和内部相邻锚点范围约束；按最近邻rank曼哈顿距离选择。
同分区块保持歧义，不强制指定。不同染色体或A/B/D亚基因组关系不会仅因名称不同而排除。

## 分阶段调度

```bash
pantrans synteny --stage index --prepared prepared --results construct_results -o synteny_results
pantrans synteny --stage shard -o synteny_results
pantrans synteny --stage pair --pair-index 0 -o synteny_results
# 对inputs.json/pairs.tsv内全部pair-index运行pair阶段后
pantrans synteny --stage collect -o synteny_results
pantrans synteny --stage export --pair-index 0 -o synteny_results
```

run为默认，顺序执行全部阶段。pair任务可独立安排到SLURM；不在命令内部自动提交。
index要求新目录；已有阶段完成标记不覆盖。阶段失败不会发布最终summary.json。
shard先核验BED和pre哈希，扫描图时记录边SHA256；不可在阶段之间更换原始输入。

## 输出

- inputs.json：品种、染色体、原始文件身份、参数和pair索引。
- pairs.tsv：全部无序品种对，8品种28对，18品种153对。
- index/、shards/：坐标/rank/pre索引与按方向分片的压缩候选。
- pairs/pair_0001/等：blocks.tsv、anchors.tsv、position_supported_pairs.tsv、summary.json。
- candidate_pairs.npy：包括未定位基因的完整候选，不静默丢弃。
- placed_candidates.tsv：默认非CS配对的主染色体候选表；涉及CS的配对默认candidates.tsv含全部候选。
  --export-candidates可对所有配对输出完整候选；export阶段也可从二进制流式导出至新文件。
- summary.tsv、summary.json：全部品种对完成且行数检查通过后发布的汇总。

数字pair目录避免品种标签被解释为路径，实际品种映射在pairs.tsv/inputs.json。
保留歧义及高拷贝信息；区块和位置支持不是一对一正交基因或确定易位事件的证明。
