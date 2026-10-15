# 候选易位结果可视化

`plot_translocations.py` 是独立脚本，不导入 PanTrans，不重新比对、聚类或判定易位。
依赖 Python >= 3.6、Matplotlib >= 3.3；无图形界面的服务器也可运行。

```bash
python3 -m pip install matplotlib
python3 scripts/plot_translocations.py --help
```

## 染色体总览

```bash
python3 scripts/plot_translocations.py \
  --blocks result.translocation.blocks.tsv \
  --chrom-sizes chromosome_sizes.tsv \
  --view overview \
  --formats pdf png svg \
  -o figures/translocations
```

`--chrom-sizes` 是**带表头**的三列 TSV，列名为 `variety`、`chrom`、`length`。
从各品种实际组装获取长度，不使用末端基因位置推算染色体长度。各品种的染色体名称
须与检测结果一致。文件可只包含需要展示的染色体，但必须覆盖所选 block 的两端。

```text
variety	chrom	length
CS	Chr1A	594102056
Q	Chr1A	590000000
Q	Chr3B	830000000
```

以上长度仅为格式示例。每个目标品种输出一张图：参考品种在上、目标品种在下。
两侧共用 bp 比例尺；同名染色体使用相同槽位起点，槽位间空隙仅为排版。
染色体区段连接带的边界、中心连线，以及局部图的基因匹配线均采用三次贝塞尔曲线。
曲线在轨道处垂直进入/离开；曲率仅用于排版，不表示生物学距离或断点路径。
绿色为正向排列，赭色为反向排列；实线为 supported_block，虚线为 low_support_block。
很短的区段以中心连线保持可见，连线宽度不表示区段长度。全长染色体的背景为灰色。

## 局部共线性图与 flank

```bash
python3 scripts/plot_translocations.py \
  --blocks result.translocation.blocks.tsv \
  --genes result.translocation.genes.tsv \
  --bed all_genes.bed \
  --chrom-sizes chromosome_sizes.tsv \
  --view local \
  --block-id TR000001 \
  --flank 100000 \
  --formats pdf png svg \
  -o figures/translocations
```

局部图必须提供 `--genes`，每个 block 分别输出；不指定 `--block-id` 时绘制全部选中区段。
`--flank` 默认 100000 bp，在参考侧和目标侧各自的区段范围左右扩展；设为 0 只显示区段。
左端截到 0；提供 `--chrom-sizes` 后右端截到实际染色体长度。
不提供长度时，图中明确提示右侧窗口可能超出组装末端，不会假设末端基因等于染色体末端。
两条轨道分别按自己的窗口线性缩放，标示各自实际 bp 坐标；**不可跨轨道比较视觉长度**。
目标轨道始终按坐标递增排列，反向 block 的连线保留交叉效果，不翻转坐标掩盖方向。

`--bed` 使用原检测的完整 BED，至少五列：chrom、start、end、gene、variety，无表头。
灰色显示窗口内背景基因；没有 BED 时仅能展示候选基因表中已有的基因，图内会提示。
默认只为区段锚点画匹配连线。没有连线不代表没有同源关系。

如需附加邻域匹配，可增加：

```bash
  --edges result.graph.edges.tsv --show-ambiguous
```

`--edges` 为无表头双列有向边 TSV，须同时提供完整 BED；仅叠加直接双向图边。
`--show-ambiguous` 独立控制是否叠加基因报告中标记 ambiguous_copy 的关系。
两者均仅绘制**两个当前窗口之间**的关系，灰色点线表示附加关系，不是新增支持锚点。
不会展示通向其他染色体或窗口外的匹配，不会据此重新判定正交关系或易位。

## 参数与输出

- `--query Q1 Q2`：选择目标品种；默认分别绘制全部。
- `--block-id TR000001 TR000002`：选择区段，适用于两种视图。
- `--status supported_block low_support_block`：默认显示两类区段。
- `--chrom-order FILE`：总览染色体顺序，一行一个名称；未列出者按自然顺序追加。
- `--label-blocks`：总览标注区段 ID；密集区域可能拥挤，默认关闭。
- `--label-genes`：局部图标注锚点基因 ID；长 ID 或密集基因可能重叠，默认关闭。
- `--width / --height`：图宽/高，单位英寸；不指定时自动选择。
- `--formats pdf png svg`：默认 PDF、PNG，可选 SVG，PDF/SVG 保留文本。
- `--dpi 300`：位图分辨率。

输出命名为 `<prefix>.<query>.overview.<format>` 或
`<prefix>.<query>.<block_id>.local.<format>`，并记录 `<prefix>.plot_manifest.json`，
包含参数、筛选前后数量、局部窗口、附加匹配数及输出文件。
空结果只生成 manifest，不生成虚构区段。重复运行同一输出前缀会覆盖同名文件；
前次运行但本次未选中的图不会自动清除，应使用新前缀或以本次 manifest 为准。

所有区间为 BED 的 0-based 半开区间。色块为锚点覆盖范围，**不是断点**。
此脚本只展示候选结果，不把 block 数量解释为已验证易位事件数。
