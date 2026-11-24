# preprocess：CDS和gene-body输出（schema 2）

预处理只进行输入解析、序列提取和CDS注释/序列完整性检查，不进行比对，不判断同源关系，也不执行构图。`valid`只表示满足本流程的完整CDS输入要求，不能理解为“通过构图校验”或生物学已验证。

以prefix=sample为例：

|文件|内容|
|---|---|
|sample_gene.fa|每基因一条，无侧翼gene-body gDNA，保留内含子，按基因方向输出|
|sample_cds.fa|全部具有CDS注释的转录本拼接CDS，包括未通过完整性检查的记录|
|sample_cds.valid.fa|通过完整CDS注释/序列检查的子集，供当前CDS实验流程作为query输入|
|sample_cds.longest.valid.fa|每个有有效CDS的基因仅一条最长有效CDS；FASTA ID保留所选输出转录本ID，供新版construct读取|
|sample_cds.longest.valid.tsv|逐基因选择记录，包含所选转录本、CDS长度、品种、染色体、链、基因组坐标、有效转录本数、最长等长候选数和状态|
|sample_mrna.fa|完整exon拼接cDNA，沿用原文件名|
|sample_gene3000bp_gdna.fa|gene+上下游最多3000bp，沿用原文件名|
|sample_cds_regions.tsv|原始基因组、带侧翼gDNA、无侧翼gene-body三种CDS区间及原始phase|
|sample_cds_status.tsv|逐转录本状态，不删除无CDS/无效CDS对应的gene、cDNA或BED|
|sample_transcript_models.jsonl|输出ID下的exon、CDS、cDNA-CDS区间和有效性状态，实验入口直接读取|
|sample_preprocess.json|schema_version=2；文件列表、CDS检查策略、各品种状态计数及空FASTA列表|

完整性检查：phase与累计编码长度一致；起始ATG和终止TAA/TAG/TGA；长度为3倍数；无内部终止；CDS仅含ACGT；CDS位于exon内，成熟cDNA中连续。该策略是当前实验的保守输入要求，不表示不完整CDS没有研究价值。内部CDS块的phase用于校验，不从每个块头部剪除phase个碱基。跨exon密码子保留。

结构损坏（CDS越出exon、重叠、孤立Parent）明确报错；partial_or_invalid_ORF、invalid_phase、internal_stop、ambiguous_CDS、no_CDS等写入状态表。无CDS无序列可输出，因此不在CDS FASTA中，但仍在其它输入中。没有任何合格CDS时preprocess仍可完成；空FASTA不生成.fai，实验构图入口报明确错误。

## 最长有效CDS输出

自动生成，无需额外开关。先排除未通过有效性检查的转录本，再按CDS长度降序、输出转录本ID字典升序选择。因此，无效的更长CDS不会取代有效候选，等长时`.10`先于`.2`。选择只影响新增FASTA，不缩减完整cDNA、全部CDS、全部有效CDS或转录本模型。

逐基因表列为`gene, transcript, CDS_length, variety, chrom, strand, gene_start0, gene_end0, valid_transcripts, max_length_ties, status`（制表符分隔）；坐标为原始基因组0-based half-open，序列均按基因方向输出。`selected`表示有且仅有一条序列进入最长CDS FASTA；`no_valid_CDS`表示有转录本但没有有效CDS；`no_transcripts`表示只有gene注释。后两者的transcript为空、CDS_length为0，仍保留基因体DNA、带侧翼DNA、BED、品种和染色体坐标；不补造转录本。仅含gene的注释也可预处理，空转录本/CDS FASTA不生成索引。空注释或损坏的转录本结构仍报错。

manifest中`files.cds_longest`和`files.cds_longest_selection`定位新文件，`longest_cds_selection`记录规则及有/无有效CDS、最长等长候选的基因计数，各品种也记录对应计数及`genes_without_transcripts`。原有`skipped_genes_without_transcripts`保留兼容但值为0，因为不再跳过这些基因。

本次扩展保持schema 2及`graph_query=cds_valid`（保留原始有效CDS数据集标记），并通过`files.cds_longest`提供独立的最长CDS输入。新版construct明确读取和校验该文件，按品种目标运行N2000比对；append亦使用相同最长CDS及N2000参数，按历史/新增目标选择查询集合。旧预处理目录不修改，需在新目录重新运行以获得新文件。默认侧翼为3000 bp，仍可用`--flank`调整。

append合并两个均带有新输出的预处理目录时，会同时合并最长CDS FASTA和选择表，并更新总计。若其中一个目录是没有这些文件的旧schema-2输入，合并结果仅保留共有的原有输出，不声明存在最长CDS文件或保留错误的选择统计；已有construct/append比对策略不受此兼容处理影响。

新CDS实验入口从manifest定位文件，不再固定five_variety前缀，也不要求再次提供原始GFF3列表。schema 1旧数据仍需--annotations。示例：

```bash
pantrans preprocess --genome-list genomes.txt --annotation-list annotations.txt --variety-list varieties.txt --prefix sample --output prepared
python scripts/experimental_cds_construct.py --prepared prepared --output results --threads 32 --workers 8
```

实验读取器直接读取预处理的坐标模型，并使用gene.fa和cds.valid.fa；不重新解析原始注释提取CDS。其它输入文件名、负链方向和contig侧翼截短规则不变。每个prepared目录只能有一个manifest，存在多个时明确报错。
