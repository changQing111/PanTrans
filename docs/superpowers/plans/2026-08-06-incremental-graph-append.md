# Incremental Graph Append Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Persist the complete graph from construct and make append reuse it while calculating only new-to-all and history-to-new alignments.

**Architecture:** A focused graph-package module owns manifest, edge, node, and source identity serialization. Construct builds the graph once, clusters it, and writes a reusable package. Append validates the package, creates full merged metadata, runs two cross-alignment blocks, merges historical and cross edges/BAM records, and calls the same graph-to-cluster/output path as construct.

**Tech Stack:** Python 3.8+, `argparse`, JSON, `pysam`, `networkx`, existing minimap2/samtools wrappers, `unittest`/`unittest.mock`.

---

## File Map

- Create `src/pantrans/graph_package.py`: graph manifest and sidecar I/O, source identity validation, BED merge, edge iteration.
- Modify `src/pantrans/align_filter.py`: merge BAM records into a chosen full target header while remapping reference IDs by name.
- Modify `src/pantrans/pipeline.py`: expose graph-plus-cluster derivation, write construct graph packages, replace representative append with incremental cross alignment and graph union, write chained graph packages.
- Modify `src/pantrans/main.py`: replace representative arguments with required `--history-graph` and pass the new append signature.
- Modify `tests/test_append_construct_flow.py`: new CLI, package, construct, append, and BAM-merge tests.
- Modify `README.md`: document graph outputs, new append command, incremental computation, and validation caveat.

## Task 1: Specify Graph Package I/O

- [ ] Add failing tests that write two directed edges, two node records, source paths, varieties, reference, main chromosomes, and a filtered BAM path; load the manifest and assert all resolved values.
- [ ] Run `PYTHONPATH=src python -m unittest tests.test_append_construct_flow.GraphPackageTest -v` and verify import/API failures.
- [ ] Implement `file_identity()`, `write_graph_package()`, `load_graph_package()`, and `iter_graph_edges()` in `graph_package.py`.
- [ ] Add tests proving a changed source size raises `ValueError`, duplicate edges serialize once, and missing node lengths fail before writing.
- [ ] Re-run `GraphPackageTest` and verify it passes.

## Task 2: Specify Complete BED Merge

- [ ] Add failing tests with a historical full BED and append representative-plus-new BED; assert the result contains historical rows once followed by new rows.
- [ ] Add failing tests for a conflicting historical duplicate, a missing new-gene BED row, and an append BED ID absent from both history and new gDNA.
- [ ] Run the focused tests and verify failures originate from the missing merge helper.
- [ ] Implement `merge_history_and_query_bed()` using parsed BED rows and exact gene IDs.
- [ ] Re-run the focused tests and verify all merge contracts pass.

## Task 3: Specify BAM Remapping

- [ ] Add a small real-pysam test with a full header `[Old.g1, New.g1]` and a subset BAM header `[New.g1]`; merge it and assert the resulting alignment still targets `New.g1`, not numeric reference ID zero.
- [ ] Run the focused test and verify `merge_bams_with_full_header()` is missing.
- [ ] Implement the helper in `align_filter.py`, copying reads into the selected header after mapping `reference_name` and `next_reference_name` to full-header IDs.
- [ ] Re-run the focused test and verify the remapped BAM record and header.

## Task 4: Make Construct Emit the Historical Graph

- [ ] Add a failing construct orchestration test asserting `<prefix>.graph.json`, edges, and nodes are requested with the original full input paths and filtered BAM.
- [ ] Refactor clustering into `_derive_graph_and_clusters()` returning `(graph, pre_clusters, last_clusters)` while retaining `_derive_clusters_from_alignment()` as a compatibility wrapper.
- [ ] Run existing shared-helper tests to ensure singleton recovery and `assign_sccs()` semantics remain unchanged.
- [ ] Update `unit_construct()` to use the returned graph and call `write_graph_package()` after the filtered graph is finalized.
- [ ] Run construct and graph-package tests.

## Task 5: Specify Incremental Append Orchestration

- [ ] Replace old representative CLI tests with a failing contract requiring `--history-graph` and rejecting `--refer_cdna`/`--refer_gdna`.
- [ ] Add a failing orchestration test whose loaded package contains full historical sequence parts, full BED, edge file, filtered BAM, variety order, reference, and main chromosomes.
- [ ] Assert append concatenates historical parts plus query for full output, runs `query cDNA -> merged gDNA` and `historical cDNA -> query gDNA`, and filters both BAMs.
- [ ] Assert graph input is the historical-edge iterator chained with both returned cross-edge lists, and construct metadata is passed unchanged into shared clustering.
- [ ] Assert the combined filtered BAM contains historical and both cross-filtered inputs and a new graph package records the expanded history.
- [ ] Run the focused append tests and verify they fail against representative-only `unit_append()`.

## Task 6: Implement Incremental Append

- [ ] Change `unit_append()` to load and validate `--history-graph`, normalize non-overlapping new variety names, and resolve historical source files.
- [ ] Build historical cDNA, merged cDNA/gDNA, and full merged BED working files.
- [ ] Run and filter the two cross alignments with existing thresholds.
- [ ] Merge filtered BAMs into the full merged target header and chain history/cross directed edges into `_derive_graph_and_clusters()`.
- [ ] Use historical plus new varieties, original reference, and stored main chromosomes for cluster assignment and rename generation.
- [ ] Preserve existing pre/last cluster and GTF/FASTA/BED output naming using full merged inputs.
- [ ] Write the expanded graph package for chained append.
- [ ] Run all append-focused tests and fix only implementation defects.

## Task 7: Documentation and Full Verification

- [ ] Update README construct outputs to include `.graph.json`, `.graph.edges.tsv`, and `.graph.nodes.tsv`.
- [ ] Replace representative-only append documentation with the approved `--history-graph` command and two-block cross-alignment explanation.
- [ ] Run `PYTHONPATH=src python -m unittest discover -s tests -v` and require zero failures.
- [ ] Run `python -m compileall -q src tests` and require exit code zero.
- [ ] Inspect `git diff --check`, `git status --short`, and the complete diff while preserving untracked `validation/` data.

## Task 8: Existing Six-Variety Validation

- [ ] Generate a graph package from the existing five-variety construct inputs and filtered BAM without repeating minimap2.
- [ ] Run incremental append for JM22 with the new graph package and record wall time and generated edge counts.
- [ ] Compare incremental directed edges and reciprocal pairs with the existing six-variety filtered BAM.
- [ ] Compare pre and last representative IDs and exact order-independent members with the six-variety construct outputs.
- [ ] Save the comparison summary under `validation/append_JM22/` without deleting previous validation artifacts.
