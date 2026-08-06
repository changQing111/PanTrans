# Append Construct-Flow Design

Date: 2026-08-06
Status: Proposed
Branch: `codex/append-construct-flow`

## 1. Goal

Rewrite `pantrans append` so that it does not consume an existing cluster file. Append will still use the representative cDNA and gDNA sequences produced by the existing pan-genome, combine them with the new variety sequences, and then run the same graph-to-cluster procedure as `construct` after alignment filtering.

This change addresses the observed inflation of append `last.cluster`: genes absent from the append alignment graph must not be recovered from the entire historical BED as singleton clusters.

## 2. Interface

The new append command accepts:

```text
pantrans append \
  --name <new-variety-names> \
  --cdna <new.cdna.fasta> \
  --gdna <new.gdna.fasta> \
  --bed <all-old-and-new-genes.bed> \
  --refer_cdna <existing-representatives.cdna.fasta> \
  --refer_gdna <existing-representatives.gdna.fasta> \
  [--bam <merged-cdna-to-merged-gdna.bam>] \
  --threads <N> \
  --prefix <prefix> \
  --output <directory>
```

Changes from the current development interface:

- Remove `--refer_cluster`. Append never reads or merges historical cluster membership.
- Remove the separate `--refer_bed`. `--bed` is the complete BED containing historical and new-variety genes.
- Keep `--refer_cdna` and `--refer_gdna`; they contain only historical representative sequences.
- Support optional `--bam` with the same meaning as construct: skip minimap2, but still run filtering and all downstream steps. The BAM header and target names must correspond to the merged gDNA input.

## 3. Input Semantics

The historical pan-genome is intentionally represented only by its representative cDNA and gDNA sequences. Therefore, append can cluster only:

1. historical representative genes present in `--refer_gdna`; and
2. new-variety genes present in the query gDNA.

The full BED is retained for chromosome assignment, strand lookup, ordering, and final BED output. It is metadata, not the source of the append gene universe. Historical non-representative genes that exist only in the BED cannot be reconstructed without `--refer_cluster` and must not appear as newly created singleton clusters.

Reference variety prefixes are inferred from the full BED by subtracting the names supplied through `--name`. These prefixes are passed to representative selection as the historical reference set. For chromosome assignment, all historical prefixes continue to act as the single logical `Refer` group, while each appended variety remains a separate group.

## 4. Processing Flow

### 4.1 Prepare and align

1. Normalize the appended variety names.
2. Concatenate representative cDNA with query cDNA.
3. Concatenate representative gDNA with query gDNA.
4. Read the complete BED directly; do not create a BED by combining separate reference and query BED files.
5. Run merged cDNA-to-merged gDNA alignment, or use `--bam`.
6. Run the same BAM filtering used by construct.

### 4.2 Shared construct clustering

Construct and append will call one shared post-filter clustering helper. The helper performs these steps without append-specific reprocessing:

1. Build the directed graph with `di_graph_from_pair()`.
2. Find strongly connected components with `get_conn_comp()`.
3. Call `assign_sccs()` once. This includes `unit_recursion()` on each non-singleton SCC, chromosome-based assignment, and the second `unit_recursion()` used by construct inside assigned chromosome groups.
4. Recover sequence-backed genes that are absent from the filtered graph as singleton clusters in both `pre_clusters` and `last_clusters`.
5. Convert the two returned lists to cluster dictionaries and write them directly.

There is no historical-cluster merge and no second append-only derivation of `last_clusters` from `pre_clusters`.

The helper receives an explicit eligible gene set so that recovery is correct in both modes:

- construct: genes selected from the BED by the requested construct variety names, preserving current construct behavior;
- append: gene IDs present in the merged gDNA FASTA and in the full BED.

The append eligible set is the critical invariant. The full BED may contain hundreds of thousands of historical non-representative genes, but those IDs cannot be recovered unless their genomic sequence was actually part of the append comparison.

### 4.3 Outputs

Append keeps the existing prefix-based output naming, including:

- `<prefix>_pre.tmp.cluster`
- `<prefix>_last.tmp.cluster`
- merged FASTA and filtered BAM intermediates
- pre and final GTF, representative cDNA, representative gDNA, and BED outputs

Cluster member order remains deterministic through the existing graph and cluster dictionary code. Representative selection continues to prefer a historical representative when graph degree and sequence length are tied.

## 5. Validation and Errors

- Fail early if no appended variety name is supplied.
- Fail if reference variety prefixes cannot be inferred from the full BED.
- Fail if an optional BAM does not exist.
- Report sequence IDs missing from the BED clearly; they cannot be assigned by chromosome or emitted consistently.
- Log the graph node count, eligible gene count, and recovered singleton count so unexpected cluster inflation is visible in normal runs.

## 6. Tests

Implementation will be driven by focused automated tests:

1. CLI parsing succeeds without `--refer_cluster` or `--refer_bed`, and accepts optional `--bam`.
2. Append merges representative and query FASTA files, then invokes the shared construct clustering path.
3. `assign_sccs()` output is used directly; append does not merge old clusters or independently rebuild last clusters.
4. A gene in merged gDNA but absent from the graph is recovered in both pre and last output.
5. A historical non-representative gene present only in the full BED is not recovered.
6. Historical reference prefixes are treated as `Refer`, and new varieties remain separately assigned.
7. Existing construct behavior remains unchanged through regression tests of the shared helper.
8. A small end-to-end fixture compares construct-on-all-sequences with append-on-representatives-plus-query under the expected information-loss constraint.

After unit tests pass, the implementation will be validated against the existing six-variety construct and five-plus-JM22 append datasets. Comparisons will use representative gene ID plus order-independent member sets for both pre and last cluster files, and will explicitly report identical, changed, construct-only, and append-only clusters.

## 7. Non-Goals

- Reconstructing historical non-representative membership without a cluster file.
- Guaranteeing byte-identical results between six-variety construct and representative-only append; append starts from a deliberately reduced historical sequence set.
- Changing the biological thresholds in alignment filtering, representative selection, chromosome assignment, or `unit_recursion()`.
