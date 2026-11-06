# preprocess：CDS和gene-body输出（schema 2）

预处理只进行输入解析、序列提取和CDS注释/序列完整性检查，不进行比对，不判断同源关系，也不执行构图。`valid`只表示满足本流程的完整CDS输入要求，不能理解为“通过构图校验”或生物学已验证。

以prefix=sample为例：

|文件|内容|
|---|---|
|sample_gene.fa|每基因一条，无侧翼gene-body gDNA，保留内含子，按基因方向输出|
|sample_cds.fa|全部具有CDS注释的转录本拼接CDS，包括未通过完整性检查的记录|
|sample_cds.valid.fa|通过完整CDS注释/序列检查的子集，供当前CDS实验流程作为query输入|
|sample_mrna.fa|完整exon拼接cDNA，沿用原文件名|
|sample_gene3000bp_gdna.fa|gene+上下游最多3000bp，沿用原文件名|
|sample_cds_regions.tsv|原始基因组、带侧翼gDNA、无侧翼gene-body三种CDS区间及原始phase|
|sample_cds_status.tsv|逐转录本状态，不删除无CDS/无效CDS对应的gene、cDNA或BED|
|sample_transcript_models.jsonl|输出ID下的exon、CDS、cDNA-CDS区间和有效性状态，实验入口直接读取|
|sample_preprocess.json|schema_version=2；文件列表、CDS检查策略、各品种状态计数及空FASTA列表|

完整性检查：phase与累计编码长度一致；起始ATG和终止TAA/TAG/TGA；长度为3倍数；无内部终止；CDS仅含ACGT；CDS位于exon内，成熟cDNA中连续。该策略是当前实验的保守输入要求，不表示不完整CDS没有研究价值。内部CDS块的phase用于校验，不从每个块头部剪除phase个碱基。跨exon密码子保留。

结构损坏（CDS越出exon、重叠、孤立Parent）明确报错；partial_or_invalid_ORF、invalid_phase、internal_stop、ambiguous_CDS、no_CDS等写入状态表。无CDS无序列可输出，因此不在CDS FASTA中，但仍在其它输入中。没有任何合格CDS时preprocess仍可完成；空FASTA不生成.fai，实验构图入口报明确错误。

新CDS实验入口从manifest定位文件，不再固定five_variety前缀，也不要求再次提供原始GFF3列表。schema 1旧数据仍需--annotations。示例：

```bash
pantrans preprocess --genome-list genomes.txt --annotation-list annotations.txt --variety-list varieties.txt --prefix sample --output prepared
python scripts/experimental_cds_construct.py --prepared prepared --output results --threads 32 --workers 8
```

实验读取器直接读取预处理的坐标模型，并使用gene.fa和cds.valid.fa；不重新解析原始注释提取CDS。其它输入文件名、负链方向和contig侧翼截短规则不变。每个prepared目录只能有一个manifest，存在多个时明确报错。
