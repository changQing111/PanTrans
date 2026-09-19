# Batched BAM Merge Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Prevent rescue BAM merges from exceeding the process file-descriptor limit by hierarchically merging at most 128 input BAMs per command.

**Architecture:** Keep `_merge_bam_files(bam_paths, merged_bam_path)` as the sole public-internal boundary. For multiple inputs, repeatedly partition paths into groups of 128, materialize each group as a temporary BAM in a `TemporaryDirectory` next to the destination, and use the last temporary BAM as the final output. A one-item group is copied to a temporary BAM; this avoids moving or mutating a caller's source input.

**Tech Stack:** Python 3.8+, `subprocess`, `tempfile`, `shutil`, `unittest.mock`, `samtools merge`.

---

### Task 1: Add a Regression Test for Command Partitioning

**Files:**
- Modify: `tests/test_append_construct_flow.py`
- Test: `tests/test_append_construct_flow.py::BatchedBamMergeTest::test_hierarchically_merges_more_than_one_batch`

- [ ] **Step 1: Write the failing test**

Add `BatchedBamMergeTest` using the existing imported `pipeline` module. Create 129 real one-byte input files in a `TemporaryDirectory`, then patch `pipeline.subprocess.run` with a side effect that creates the requested command output. Assert that `_merge_bam_files()` returns `True`, creates the requested destination, executes two merge commands across two rounds, and that every command has at most 128 input paths:

```python
input_count = len(command) - 4  # samtools merge -f OUTPUT INPUT...
self.assertLessEqual(input_count, 128)
```

- [ ] **Step 2: Run the test to verify it fails**

Run:

```bash
python -m unittest tests.test_append_construct_flow.BatchedBamMergeTest.test_hierarchically_merges_more_than_one_batch -v
```

Expected: FAIL because the current implementation makes one `samtools merge` call with 129 inputs.

### Task 2: Implement Bounded Hierarchical Merge

**Files:**
- Modify: `src/pantrans/pipeline.py:426-433`
- Test: `tests/test_append_construct_flow.py::BatchedBamMergeTest::test_hierarchically_merges_more_than_one_batch`

- [ ] **Step 1: Add the bounded batch constant and path grouping loop**

Import `tempfile` and define a private batch size constant of `128`. Preserve these existing early returns:

```python
if not bam_paths:
    return False
if len(bam_paths) == 1:
    shutil.copyfile(bam_paths[0], merged_bam_path)
    return True
```

For multiple paths, create a `tempfile.TemporaryDirectory` adjacent to the absolute destination directory with a stable `<output-stem>.merge-` prefix. While more than one path remains, partition the current paths into groups of at most 128. Copy one-item groups to a uniquely named temporary BAM; merge larger groups with:

```python
subprocess.run(["samtools", "merge", "-f", intermediate_path] + group, check=True)
```

Replace `current_paths` with the generated intermediate paths after every round. Finally use `os.replace()` to place the last intermediate at `merged_bam_path`. Let `TemporaryDirectory` clean intermediates if a merge command fails.

- [ ] **Step 2: Run the regression test to verify it passes**

Run:

```bash
python -m unittest tests.test_append_construct_flow.BatchedBamMergeTest.test_hierarchically_merges_more_than_one_batch -v
```

Expected: PASS; the two `samtools merge` calls use 128 and 2 inputs respectively.

### Task 3: Lock Compatibility and Failure Cleanup

**Files:**
- Modify: `tests/test_append_construct_flow.py`
- Test: `tests/test_append_construct_flow.py::BatchedBamMergeTest`

- [ ] **Step 1: Add compatibility and cleanup tests**

Add focused tests that verify:

- empty input returns `False` without calling `subprocess.run`;
- one input is copied verbatim without calling `subprocess.run`;
- the successful 129-input merge leaves no `<output-stem>.merge-*` temporary directory;
- a mocked failing `samtools merge` raises `subprocess.CalledProcessError` and leaves no temporary merge directory next to the destination.

- [ ] **Step 2: Run the focused test class**

Run:

```bash
python -m unittest tests.test_append_construct_flow.BatchedBamMergeTest -v
```

Expected: PASS.

### Task 4: Verify the Repository

**Files:**
- Verify: `src/pantrans/pipeline.py`
- Verify: `tests/test_append_construct_flow.py`

- [ ] **Step 1: Run all tests**

Run:

```bash
python -m unittest discover -s tests -v
```

Expected: PASS.

- [ ] **Step 2: Check the changed files**

Run:

```bash
git diff -- src/pantrans/pipeline.py tests/test_append_construct_flow.py
```

Expected: only the bounded merge implementation and its tests. The repository is currently not a Git worktree, so record that limitation rather than attempting a commit.
