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
下方为之前版本的迁移及验证记录，命令名仅作历史参考。

---

# Construct/append CLI migration in 0.8.0

**2026-11-23:** construct has moved to longest-valid-CDS, per-variety N2000
mapping. Append now shares this profile and accepts compatible N2000 graph packages.
Older N100 histories must be rebuilt rather than imported under the new policy.
The compatibility/equivalence checks below describe the earlier N100 workflow.
See [current construct](construct-longest-cds.md).

The official installed commands now use the CDS workflow:

| Previous entry | 0.8 command |
|---|---|
| `python scripts/experimental_cds_construct.py --prepared P --output O` | `pantrans construct --prepared P --output O` |
| `python scripts/experimental_cds_append.py ...` | `pantrans append ...` |
| `pantrans construct --cdna C --gdna G --bed B ...` | `pantrans legacy-construct --cdna C --gdna G --bed B ...` |
| `pantrans append --history-graph ... --cdna ...` | `pantrans legacy-append --history-graph ... --cdna ...` |

This is an explicit CLI contract change. Legacy arguments are not accepted or
silently routed by the new default commands. Historical function implementations
remain in `pipeline.py`; the new implementations are `pantrans.cds_construct`
and `pantrans.cds_append`. Experimental script paths are forwarding wrappers,
not separate algorithm implementations. Running immutable historical source
snapshots continues to use their original code.

`preprocess` and `cnv` retain their interfaces. The `translocation` command
has been removed; `synteny` now provides all-pair gene-order block analysis
from complete prepared inputs, pre clusters and graph edges (see synteny.md). Legacy experiment
runners/validation scripts now call the explicit legacy subcommands. Full new
constructs publish `cds.graph.json`; append checks schema, hashes and policy and
can chain through its `merged_prepared/` directory.

The CDS v2 graph uses 90% identity and both CDS-overlap coverages, without coding
splice/indel vetoes. Only last is annotated by default. UTR-aware structural
deduplication, persistent readers, bounded scheduling and transfer-evidence caches
are shared between new construct and append. Annotation transfer itself retains
strict coding-structure and ORF checks; the relaxed graph filter is not reused as
an annotation acceptance rule.

The current reference/chromosome policy remains wheat-specific (CS, Chr1A–Chr7D).
It was not generalized as part of this CLI migration. The independent ABD-only
partition experiment is not enabled by these commands.

Compatibility script options remain supported, including optional schema-1
`--annotations`, `--genes`, `--annotate-pre`, and annotation-only `--reuse-graph`.
Schema-2 full runs are the supported source of graph packages for append.
Packaged graph restarts validate current graph/input policy; unpackaged restarts
record unverified graph thresholds rather than attributing the current defaults
to an old graph. An append graph consists of multiple generations' alignments,
so the older construct raw-BAM restart interface is not a general append resume.

Validation includes installed-wheel console entrypoints outside the repository,
real minimap2 construct/append runs, repeat append, source/installed output
comparison and legacy parser regressions. Benchmark datasets and generated
experiment output are not shipped in the Python package.


Validation for this migration: 252 tests and 837 subtests passed on Linux with
real minimap2/samtools integration enabled. A built 0.8.0 wheel was installed in
an isolated virtual environment and its console commands were exercised from a
temporary directory with PYTHONPATH removed. Installed preprocess, construct,
append and a direct combined construct produced matching graph/cluster/GFF3/cDNA
and accepted alignment outputs on the smoke dataset. Console --version and
explicit legacy help were checked. This verifies package/CLI integration, not
universal full-cohort equivalence of incremental and direct alignment.
