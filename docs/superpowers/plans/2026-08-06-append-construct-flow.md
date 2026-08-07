# Append Construct-Flow Implementation Plan (Superseded)

> This plan describes the earlier representative-only append implementation and
> is retained as development history. The active implementation plan is
> `2026-08-06-incremental-graph-append.md`.

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace append's historical-cluster merge with the same post-filter graph/SCC/assignment/recovery flow used by construct, using pre-cluster representative BED plus new-variety BED and no `--refer_cluster`.

**Architecture:** Extract the post-filter graph-to-cluster and independent pre/last missing-gene recovery into one pipeline helper. Construct passes its existing BED-selected input gene set; append passes the merged gDNA/append-BED intersection. Append concatenates representative and query FASTA files, copies the already-combined BED to its output workspace, filters either a new or supplied merged BAM, calls the shared helper once, and writes the returned clusters directly.

**Tech Stack:** Python 3.8+, `argparse`, `unittest`/`unittest.mock`, existing `networkx`, `pysam`, Biopython, minimap2/samtools runtime tools.

---

## File Map

- Modify `src/pantrans/main.py`: remove `--refer_cluster` and `--refer_bed`, add optional append `--bam`, pass the combined BED and BAM to the new `unit_append` signature, and allow `read_parameters(argv=None)` for parser tests.
- Modify `src/pantrans/pipeline.py`: add shared post-filter clustering, append BED/reference-prefix validation, refactor construct to use the helper, and rewrite append orchestration.
- Create `tests/test_append_construct_flow.py`: dependency-isolated unit tests for parser behavior, shared recovery, append input validation, and append orchestration.
- Modify `README.md`: document the new append command and the representative-plus-query BED contract.

## Task 1: Add failing tests for shared post-filter clustering

**Files:**
- Create: `tests/test_append_construct_flow.py`

- [ ] **Step 1: Create dependency-isolated test scaffolding.**

Create a test module that inserts minimal `pysam`, `networkx`, and `Bio` modules into `sys.modules` before importing `pantrans`. The tests patch graph/alignment/output boundaries, so they do not require minimap2, BAM parsing, or FASTA indexing.

```python
import os
import sys
import tempfile
import types
import unittest
from unittest import mock

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")

def install_dependency_stubs():
    sys.modules.setdefault("pysam", types.ModuleType("pysam"))
    networkx_stub = types.ModuleType("networkx")
    networkx_stub.DiGraph = lambda: None
    networkx_stub.strongly_connected_components = lambda graph: []
    sys.modules.setdefault("networkx", networkx_stub)
    bio_stub = types.ModuleType("Bio")
    bio_stub.__path__ = []
    sys.modules.setdefault("Bio", bio_stub)
    for name in ("Bio.SeqIO", "Bio.Seq", "Bio.SeqRecord"):
        sys.modules.setdefault(name, types.ModuleType(name))
    sys.modules["Bio.Seq"].Seq = str
    sys.modules["Bio.SeqRecord"].SeqRecord = object

install_dependency_stubs()
sys.path.insert(0, SRC)
from pantrans import pipeline
from pantrans import main as pantrans_main
```

- [ ] **Step 2: Write the failing shared-helper test.**

Add a test with mocked `di_graph_from_pair`, `get_conn_comp`, and `assign_sccs`. The mock returns a pre cluster containing `Ref.g1` and `New.g1`, but a last cluster containing only `Ref.g1`. The eligible set contains `Ref.g1`, `New.g1`, `New.g2`, and excludes BED-only `Ref.nonrep`. Assert that the helper adds `New.g2` independently to both outputs, does not duplicate `New.g1`, excludes `Ref.nonrep`, and calls `assign_sccs` exactly once.

```python
def test_shared_clustering_recovers_pre_and_last_independently():
    graph = mock.Mock()
    graph.number_of_nodes.return_value = 2
    with mock.patch.object(pipeline, "di_graph_from_pair", return_value=graph), \
         mock.patch.object(pipeline, "get_conn_comp", return_value=[("Ref.g1", "New.g1")]), \
         mock.patch.object(
             pipeline,
             "assign_sccs",
             return_value=([['Ref.g1', 'New.g1']], [['Ref.g1']]),
         ) as assign_mock:
        pre, last = pipeline._derive_clusters_from_alignment(
            aligned_gene_li=[("Ref.g1", "New.g1")],
            gene_len_dic={"Ref.g1": 100, "New.g1": 90, "New.g2": 80},
            bed_dic={"Ref.g1": [], "New.g1": [], "New.g2": [], "Ref.nonrep": []},
            variety_li=["New"],
            refer_name="Ref",
            eligible_gene_set={"Ref.g1", "New.g1", "New.g2"},
        )

    assign_mock.assert_called_once()
    assert pre == [["Ref.g1", "New.g1"], ["New.g2"]]
    assert last == [["Ref.g1"], ["New.g1"], ["New.g2"]]
```

- [ ] **Step 3: Run the test and verify the intended failure.**

Run:

```bash
PYTHONPATH=src python -m unittest tests.test_append_construct_flow.AppendFlowUnitTest.test_shared_clustering_recovers_pre_and_last_independently -v
```

Expected: `AttributeError` because `_derive_clusters_from_alignment` does not exist on the clean baseline branch.

## Task 2: Implement and integrate the shared clustering helper

**Files:**
- Modify: `src/pantrans/pipeline.py` near `_collect_input_genes()` and `unit_construct()`
- Test: `tests/test_append_construct_flow.py`

- [ ] **Step 1: Add the minimal helper required by the failing test.**

Add:

```python
def _derive_clusters_from_alignment(
    aligned_gene_li,
    gene_len_dic,
    bed_dic,
    variety_li,
    refer_name,
    eligible_gene_set,
    main_chroms=None,
    refer_prefixes=None,
):
    graph = di_graph_from_pair(aligned_gene_li)
    sccs = get_conn_comp(graph)
    logger.info(
        "Cluster graph contains %d nodes; eligible gene set contains %d genes.",
        graph.number_of_nodes(),
        len(set(eligible_gene_set)),
    )
    pre_clusters, last_clusters = assign_sccs(
        sccs,
        graph,
        gene_len_dic,
        bed_dic,
        variety_li,
        refer_name,
        main_chroms=main_chroms,
        refer_prefixes=refer_prefixes,
    )
    pre_clustered = {gene for cluster in pre_clusters for gene in cluster}
    last_clustered = {gene for cluster in last_clusters for gene in cluster}
    missing_pre = [[gene] for gene in sorted(set(eligible_gene_set) - pre_clustered)]
    missing_last = [[gene] for gene in sorted(set(eligible_gene_set) - last_clustered)]
    if missing_pre or missing_last:
        logger.warning(
            "Recovered %d pre-cluster genes and %d last-cluster genes absent from the filtered graph.",
            len(missing_pre),
            len(missing_last),
        )
        pre_clusters.extend(missing_pre)
        last_clusters.extend(missing_last)
    return pre_clusters, last_clusters
```

- [ ] **Step 2: Run the focused test and verify it passes.**

Run the command from Task 1. Expected: one passing test, with `assign_sccs` called once.

- [ ] **Step 3: Refactor construct to use the helper.**

In `unit_construct()`, replace the local graph/SCC/`assign_sccs()`/missing-gene block with:

```python
input_gene_set = _collect_input_genes(bed_dic, variety_li)
pre_clusters, last_clusters = _derive_clusters_from_alignment(
    aligned_gene_li=aligned_gene_li,
    gene_len_dic=gene_len_dic,
    bed_dic=bed_dic,
    variety_li=variety_li,
    refer_name=refer_name,
    eligible_gene_set=input_gene_set,
    main_chroms=main_chroms,
)
```

Leave output writing, rename-map generation, rescue, and construct CLI behavior unchanged.

- [ ] **Step 4: Add a construct regression test for independent recovery.**

Patch `filter_bam`, `assign_sccs`, and `_write_cluster_outputs`, run `unit_construct()` on three small BED/FASTA fixture IDs, and assert a gene present in pre but absent in last is added only to last while a graph-absent eligible gene is added to both. This protects the shared helper from reintroducing the old “derive pre missing from last” behavior.

- [ ] **Step 5: Run all tests written so far.**

Run:

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
```

Expected: all parser-independent helper and construct regression tests pass.

## Task 3: Add failing tests for append input semantics and orchestration

**Files:**
- Test: `tests/test_append_construct_flow.py`

- [ ] **Step 1: Add reference-prefix inference tests.**

Use a combined BED containing `RefA.g1`, `JM22A1.g1`, and `JM22Ctg52.g2`. Assert `_infer_reference_variety_names_from_bed(path, ["JM22"]) == ["RefA"]`; this verifies query prefixes are excluded even when chromosome suffixes are embedded in the gene prefix.

- [ ] **Step 2: Add the append gene-set validation test.**

Call `_validate_append_gene_set()` with merged gDNA IDs `{"RefA.g1", "JM22.g1"}` and BED IDs containing one extra `RefA.nonrep`. Assert the returned eligible set contains only the two sequence-backed IDs. Add a second assertion that a sequence ID absent from BED raises `ValueError` containing `"missing from append BED"`.

- [ ] **Step 3: Add the append orchestration test.**

Create temporary representative/query FASTA paths, a combined representative-plus-query BED, and an existing BAM path. Patch `concat_fasta_files`, `get_fasta_len`, `filter_bam`, `_derive_clusters_from_alignment`, and `_write_cluster_outputs`. Call `unit_append()` with `bam_path` and assert:

```python
pipeline.minimap2_map.assert_not_called()
pipeline._derive_clusters_from_alignment.assert_called_once()
call.kwargs["eligible_gene_set"] == {"RefA.g1", "JM22.g1"}
open(os.path.join(out_dir, "Append_merged.bed")).read() == open(input_bed).read()
```

Also assert no `refer_cluster_path` argument is accepted by the new function signature.

Add a separate test passing a nonexistent `bam_path` and assert `unit_append()` raises `FileNotFoundError` before `filter_bam()` is called.

- [ ] **Step 4: Run the new append tests and verify they fail for missing behavior.**

Run:

```bash
PYTHONPATH=src python -m unittest tests.test_append_construct_flow.AppendFlowUnitTest -v
```

Expected failures: missing reference-prefix/validation helpers and the baseline `unit_append()` requiring the old cluster argument.

## Task 4: Implement new append pipeline behavior

**Files:**
- Modify: `src/pantrans/pipeline.py`
- Test: `tests/test_append_construct_flow.py`

- [ ] **Step 1: Add reference-prefix and BED/gDNA validation helpers.**

Implement:

```python
def _infer_reference_variety_names_from_bed(bed_path, query_variety_li):
    query_prefixes = tuple(query_variety_li)
    reference_names = [
        name for name in _infer_variety_names_from_bed(bed_path)
        if not name.startswith(query_prefixes)
    ]
    if not reference_names:
        raise ValueError(
            f"No reference variety names could be inferred from {bed_path}; "
            "append --bed must contain pre-cluster representatives and query genes."
        )
    return reference_names

def _validate_append_gene_set(sequence_gene_set, bed_dic):
    sequence_gene_set = set(sequence_gene_set)
    bed_gene_set = set(bed_dic)
    missing_bed = sorted(sequence_gene_set - bed_gene_set)
    if missing_bed:
        raise ValueError(
            "The following merged gDNA genes are missing from append BED: "
            + ", ".join(missing_bed[:10])
        )
    bed_only = sorted(bed_gene_set - sequence_gene_set)
    if bed_only:
        logger.warning(
            "Ignoring %d BED genes without merged gDNA sequence; they will not form clusters.",
            len(bed_only),
        )
    return sequence_gene_set
```

- [ ] **Step 2: Run the validation tests and verify they pass.**

Run the focused append test command from Task 3. Expected: prefix and validation tests pass; orchestration remains red until the next step.

- [ ] **Step 3: Rewrite `unit_append()` with the new signature.**

Use this signature:

```python
def unit_append(
    query_cdna_path,
    query_gdna_path,
    all_bed_path,
    refer_cdna_path,
    refer_gdna_path,
    bam_path,
    variety_name,
    threads,
    out_dir,
    prefix="Append",
):
```

The body must:

1. Normalize names and create output paths for merged cDNA/gDNA/BED/BAM and filtered BAM.
2. Concatenate only representative and query FASTA files.
3. Copy `all_bed_path` to `<prefix>_merged.bed`; do not concatenate reference/query BEDs.
4. If `bam_path` is supplied but does not exist, raise `FileNotFoundError`. Otherwise use the supplied BAM or call `minimap2_map(merged_cdna_path, merged_gdna_path, threads, merged_bam_path)`.
5. Load the copied BED, cDNA lengths, and filtered BAM.
6. Infer historical reference prefixes from the combined BED.
7. Set `eligible_gene_set = _validate_append_gene_set(get_fasta_len(merged_gdna_path), bed_dic)`.
8. Call `_derive_clusters_from_alignment()` exactly once with `variety_li=new_variety_li`, `refer_prefixes=refer_variety_li`, and `main_chroms=None`.
9. Write the returned pre and last cluster dictionaries directly; do not read a cluster file, call `derive_last_clusters_from_pre()`, or merge historical clusters.
10. Preserve existing transcript/FASTA/BED output calls and return the five intermediate paths.

- [ ] **Step 4: Run append tests until green.**

Run:

```bash
PYTHONPATH=src python -m unittest tests.test_append_construct_flow.AppendFlowUnitTest -v
```

Expected: all append helper/orchestration tests pass and `minimap2_map` is skipped when `--bam` is supplied.

## Task 5: Update CLI and documentation

**Files:**
- Modify: `src/pantrans/main.py`
- Modify: `README.md`
- Test: `tests/test_append_construct_flow.py`

- [ ] **Step 1: Add the failing parser test.**

Call `read_parameters([...])` with `append`, `--bed`, `--refer_cdna`, `--refer_gdna`, and optional `--bam`, but no cluster argument. Assert the parsed namespace has `bam`, has no `refer_cluster` or `refer_bed`, and retains the supplied BED path.

- [ ] **Step 2: Implement parser changes.**

Change `read_parameters()` to `read_parameters(argv=None)` and call `parser.parse_args(argv)`. For append, update the BED help to “pre-cluster representative plus new-variety BED”, remove `--refer_cluster` and `--refer_bed`, add optional `--bam`, and pass `all_bed_path=args.bed` and `bam_path=args.bam` to `unit_append()`.

- [ ] **Step 3: Run the parser test and command help.**

Run:

```bash
PYTHONPATH=src python -m unittest tests.test_append_construct_flow.AppendFlowUnitTest.test_append_parser_uses_combined_bed_without_cluster -v
PYTHONPATH=src python -m pantrans.main append --help
```

Expected: parser test passes; help contains `--refer_cdna`, `--refer_gdna`, and `--bam`, and does not contain `--refer_cluster` or `--refer_bed`.

- [ ] **Step 4: Update README append examples and caveats.**

Replace the old command that passes `--refer_cluster` with a command using `--bed <pre-representative-plus-new-variety.bed>`, `--refer_cdna`, `--refer_gdna`, and optional `--bam`. Explain that the historical side of BED is the first gene ID of each pre cluster, that old non-representatives are intentionally absent, and that append reruns graph/SCC/`assign_sccs()` directly.

## Task 6: Full verification and versioned handoff

**Files:**
- Modify: `tests/test_append_construct_flow.py` if a regression test is needed.
- Modify: `README.md` only for verified behavior corrections.

- [ ] **Step 1: Run syntax and unit verification.**

Run:

```bash
python -m compileall -q src tests
PYTHONPATH=src python -m unittest discover -s tests -v
git diff --check
```

Expected: compile succeeds, every test passes, and `git diff --check` emits no diagnostics. Do not use the currently crashing pytest command as the release gate.

- [ ] **Step 2: Run append against the existing BAM-backed fixture.**

Use the existing representative FASTA, JM22 FASTA, combined representative-plus-JM22 BED, and merged BAM under `/data/changq/PanTrans/test/`. Write all outputs under the isolated worktree, not the main test directories:

```bash
PYTHONPATH=src python -m pantrans.main append \
  --name JM22 \
  --cdna /data/changq/PanTrans/test/JM22_cdna.fasta \
  --gdna /data/changq/PanTrans/test/JM22.gdna.fasta \
  --bed /data/changq/PanTrans/test/test_pantrans_pre.refer_append_JM22.bed \
  --refer_cdna /data/changq/PanTrans/test/pantrans_construct_laste/test_pantrans_pre_cdna.refer.fasta \
  --refer_gdna /data/changq/PanTrans/test/pantrans_construct_laste/test_pantrans_pre_gdna.refer.fasta \
  --bam /data/changq/PanTrans/test/pantrans_append_laste_pre_refer_append_JM22/test_append_JM22_merged_cdna_align_gdna.bam \
  --threads 1 \
  --prefix append_construct_flow_JM22 \
  --output /data/changq/PanTrans/.worktrees/append-construct-flow/validation/append_JM22
```

Check that both cluster files contain only IDs present in the combined BED and that no historical BED-only IDs are introduced. Compare cluster files by representative first column and order-independent member sets; record counts and any changed representatives without overwriting existing test outputs.

- [ ] **Step 3: Commit in reviewable changes.**

Use separate commits for the test/helper refactor, append/CLI implementation, and documentation/verification adjustments:

```bash
git add tests/test_append_construct_flow.py src/pantrans/pipeline.py src/pantrans/main.py
git commit -m "feat: route append through construct clustering flow"
git add README.md
git commit -m "docs: document representative BED append inputs"
```

Before reporting completion, verify the worktree is clean and the main worktree still contains its pre-existing unrelated modifications.
