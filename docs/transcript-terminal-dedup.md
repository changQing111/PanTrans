# D/E：有端部证据的转录本去冗余

2026-11-06。适用于 scripts/experimental_cds_construct.py 的CDS构图＋完整cDNA注释流程；旧 pantrans construct/append 的传统去冗余入口未迁移，不应混称为全局规则变更。

先冻结完整、非CDS-only模型的精确结构类别（完整内含子链+CDS坐标）；它们按原优先级选择代表。随后逐一匹配CDS-only或partial候选，候选不能作为桥梁将不同完整类别合并。

兼容条件：同一代表局部坐标，CDS区间完全一致；完整模型的外显子截取到较短模型实际覆盖区间后，必须与较短模型外显子完全相等；向区间外扩展的端部必须有明确缺失信息。CDS-only定义为整个exon区间列表与CDS列表相同，其UTR两端视为未注释，不解释为真实不存在UTR。跨品种partial的缺失端部从接受比对的CIGAR首尾S/H读取。历史模型若没有端部字段，除精确已观察区间兼容外不推测缺失方向。

- 唯一完整结构可解释：仅输出该完整结构代表，较短来源记为unique_partial_compatible。
- 多个完整结构兼容：较短模型保留为ambiguous，compatible_sources记录候选完整来源，不桥接合并。
- 没有兼容完整结构：较短模型保留为unresolved。
- 观察区域内有独特剪接、内含子保留、CDS变化：不兼容，保留。
- 无CDS模型不走D/E兼容折叠；沿用保守精确分组。
- 不按字符串包含合并，不增设剪接位点容差，不补接未映射序列。

输出GFF3新增utr_evidence=observed/cds_only/partial及redundancy_status=observed/ambiguous/unresolved。completeness仍表示原定义的比对完整性/原始模型种子，不声称生物学全长。model_sources.tsv保留原四列，追加utr_evidence、redundancy_status、compatible_sources和source_relations_json，明确每个来源为何保留或被折叠。

规则位于src/pantrans/transcript_redundancy.py，实验脚本调用此模块。已经完成的实验snapshot、BAM、GFF3和缓存不被更新；新的规则需要新运行才生效。旧缓存缺少missing_5prime/3prime字段时不会被擅自补推，若要重跑完整D/E规则应重新产生对应候选元数据。
