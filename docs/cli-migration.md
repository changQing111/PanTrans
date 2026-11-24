# 当前CLI及历史迁移说明

当前`pantrans`仅有五个子命令：`preprocess`、`construct`、`append`、`cnv`、`synteny`。
`legacy-construct`、`legacy-append`和`translocation`已从参数注册和执行分支移除，
调用这些名称（包括--help）将报未知子命令，不重定向到新算法。

construct/append均采用最长有效CDS、按目标品种、f20000 N2000的CDS图策略。
输入和运行示例见[README](../README.zh-CN.md)、[construct](construct-longest-cds.md)、
[append](cds-append.md)。共线性分析见[synteny](synteny.md)。
旧cDNA参数不适用于当前construct/append；旧N100或无完整可验证图包的数据需要重新构建。
旧版内部Python实现和历史协议文档保留用于读取历史结果和回归验证，不构成当前命令入口。
固定的历史源码快照不因本次修改而变化。

旧版cDNA专用实验预设和SLURM验证脚本已标记退役，不能用于创建当前构建任务。

## 当前迁移方式

旧版cDNA命令没有同参数兼容入口。请从基因组和GFF3重新运行preprocess，获得
最长有效CDS及完整模型后，使用construct建立新版图。新增品种使用相同预处理策略后
通过append扩展，连续追加使用上一代merged_prepared/作为历史输入。

实验脚本experimental_cds_construct.py和experimental_cds_append.py仍是对应正式
子命令的转发入口，不是单独的算法实现。当前新构图需要新版schema-2输入；注释重启
仍通过--reuse-graph使用与历史图对应的完整预处理数据，不重新比对构图。

construct/append共用90/90/90 CDS同源过滤，默认只注释last。cDNA转移仍采用更严格的
编码结构和ORF检查，之后重新进行UTR去冗余。完整图包绑定输入/边/策略哈希，保留
父代来源；不能将旧N100图或不兼容输入伪装为新版N2000结果。

共线性命令读取完整预处理数据、pre.cluster和graph.edges.tsv，生成两两品种区块、
锚点和位置支持，不修改pre/last。当前参考/主染色体默认仍为小麦CS及Chr1A–Chr7D。
