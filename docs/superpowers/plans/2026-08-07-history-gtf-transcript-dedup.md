# Historical GTF-Aware Append Transcript Deduplication Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make construct retain an unrenamed final transcriptome and make append reuse the previous unrenamed GTF/cDNA during final transcript deduplication.

**Architecture:** Extend transcript processing with validated GTF transcript models and deterministic GTF ID renaming. Generate the unrenamed final output first in both workflows, derive the renamed output from it, and pass historical GTF models as authoritative seed candidates only to append's last-cluster deduplication.

**Tech Stack:** Python 3, argparse, pysam, Biopython, unittest/pytest, existing PanTrans graph and transcript-processing modules.

---

### Task 1: Historical GTF model parsing and seeded deduplication

**Files:**
- Modify: `src/pantrans/transcript_processor.py`
- Create: `tests/test_historical_gtf_transcripts.py`

- [ ] **Step 1: Write failing tests for GTF model loading**

Create a small unrenamed GTF containing transcript and exon rows and assert that a new `load_gtf_transcript_models()` function returns transcript-to-gene, splice-site, and exon-coordinate mappings. Include a malformed transcript whose `gene_id` differs from column 1 and assert a clear `ValueError`.

```python
models = load_gtf_transcript_models(history_gtf)
self.assertEqual(models["transcript_gene"]["Ref.g1.1"], "Ref.g1")
self.assertEqual(models["splice_sites"]["Ref.g1.1"], [(11, 20)])
self.assertEqual(models["exon_coords"]["Ref.g1.1"], [(1, 10), (21, 30)])
```

- [ ] **Step 2: Run the focused tests and verify RED**

Run:

```bash
/data/changq/miniforge3/bin/python -m pytest tests/test_historical_gtf_transcripts.py -q
```

Expected: import or attribute failure because `load_gtf_transcript_models` does not exist.

- [ ] **Step 3: Implement GTF model loading**

Add a parser that validates transcript/exon attributes, groups exons by transcript, sorts coordinates, derives introns as `(previous_exon_end + 1, next_exon_start - 1)`, and rejects inconsistent gene IDs or mixed strands/sequence names.

- [ ] **Step 4: Write and verify a failing seeded-dedup test**

Create a BAM containing only a new-variety transcript aligned to `Ref.g1`, provide a historical GTF model for `Ref.g1.1`, and call:

```python
transcript_dedup(
    bam_path,
    cluster_dic={"Ref.g1": ["Ref.g1", "JM22.g1"]},
    trans_len_dic={"Ref.g1.1": 20, "JM22.g1.1": 25},
    gene_len_dic={"Ref.g1": 100},
    gene_strand_dic={"Ref.g1": "+"},
    rename_map=None,
    gtf_path=output_gtf,
    seed_gtf_path=history_gtf,
)
```

Assert that historical and new splice patterns both appear. Expected RED: `transcript_dedup` does not accept `seed_gtf_path`.

- [ ] **Step 5: Implement seeded candidate merging**

Add optional `seed_gtf_path=None`. Merge seed splice/exon dictionaries after BAM extraction so seed records win transcript-ID collisions. Include only seed genes that are current cluster keys, log the skipped count, then use the existing `get_last_trans()` and `generate_gtf()` flow.

- [ ] **Step 6: Run focused and transcript tests**

Run:

```bash
/data/changq/miniforge3/bin/python -m pytest tests/test_historical_gtf_transcripts.py -q
```

Expected: all tests pass.

### Task 2: Unrenamed and renamed final GTF/cDNA outputs

**Files:**
- Modify: `src/pantrans/transcript_processor.py`
- Modify: `src/pantrans/pipeline.py`
- Modify: `tests/test_historical_gtf_transcripts.py`
- Modify: `tests/test_append_construct_flow.py`

- [ ] **Step 1: Write a failing deterministic-renaming test**

Test a new `rename_gtf_ids(source_path, output_path, rename_map)` function. The output must keep column 1 and coordinates unchanged while converting:

```text
gene_id "Ref.g1"; transcript_id "Ref.g1.2";
```

to:

```text
gene_id "Pan1A000001"; transcript_id "Pan1A000001.2";
```

Reject a source transcript whose numeric suffix cannot be preserved.

- [ ] **Step 2: Run the renaming test and verify RED**

Run the focused test and confirm failure because the helper is missing.

- [ ] **Step 3: Implement deterministic renaming**

Parse attributes structurally, require the unrenamed `gene_id` to match the sequence name, require `<gene_id>.<numeric_suffix>`, replace only `gene_id` and `transcript_id`, and write the destination atomically enough for the existing local workflow.

- [ ] **Step 4: Write a failing output-orchestration test**

Call `_write_cluster_outputs()` with a non-empty rename map while patching `transcript_dedup`, rescue, FASTA, and BED helpers. Assert:

- dedup is called once with `rename_map=None` and `<prefix>_unrenamed.gtf`;
- `rename_gtf_ids` creates `<prefix>.gtf`;
- cDNA extraction runs for both `<prefix>_unrenamed_cdna.refer.fasta` and the existing official cDNA path;
- the function still returns the existing four official output paths.

- [ ] **Step 5: Implement dual final output orchestration**

When `rename_map` is present, write/rescue the unrenamed GTF first, rename it to the official GTF, sort the official GTF, and derive both cDNA FASTAs. Keep pre output behavior and the four-value return contract unchanged. Add optional `seed_gtf_path` forwarding to `transcript_dedup`.

- [ ] **Step 6: Run focused and existing construct-flow tests**

Run:

```bash
/data/changq/miniforge3/bin/python -m pytest tests/test_historical_gtf_transcripts.py tests/test_append_construct_flow.py -q
```

Expected: all tests pass.

### Task 3: Append historical transcriptome input contract

**Files:**
- Modify: `src/pantrans/main.py`
- Modify: `src/pantrans/pipeline.py`
- Modify: `tests/test_append_construct_flow.py`
- Modify: `tests/test_incremental_graph_append.py`

- [ ] **Step 1: Write failing CLI tests**

Require this append contract:

```bash
pantrans append \
  --name JM22 \
  --cdna history_unrenamed_cdna.refer.fasta \
  --history-gtf history_unrenamed.gtf \
  --query-cdna JM22_cdna.fasta \
  --gdna JM22.gdna.fasta \
  --bed combined.bed \
  --history-graph history.graph.json \
  --output append_out
```

Assert all three transcript inputs are exposed under distinct argument names and omission of either new option is rejected.

- [ ] **Step 2: Run CLI tests and verify RED**

Expected: argparse rejects `--history-gtf` and `--query-cdna` as unknown.

- [ ] **Step 3: Write failing historical transcriptome validation tests**

Add `_validate_history_transcriptome(history_gtf_path, history_cdna_path, history_gene_ids)` tests for exact transcript-set agreement, graph membership, `gene_id`/sequence-name agreement, and rejection of `Pan...` renamed IDs.

- [ ] **Step 4: Implement CLI and validation**

Change `unit_append` parameters to begin with `history_cdna_path`, `history_gtf_path`, `query_cdna_path`, and `query_gdna_path`. Parse `--cdna` as history cDNA, add required `--history-gtf` and `--query-cdna`, and call the validator before creating BAMs.

- [ ] **Step 5: Write a failing append data-flow test**

Update the incremental append mock test to assert:

- query-to-all minimap2 receives `query_cdna_path`;
- history-to-query minimap2 receives `history_cdna_path`, not graph-package `cdna_paths`;
- merged working cDNA concatenates `[history_cdna_path, query_cdna_path]`;
- final `_write_cluster_outputs` receives `seed_gtf_path=history_gtf_path`;
- pre output receives no seed GTF;
- append graph metadata records the historical nonredundant and query cDNA paths.

- [ ] **Step 6: Implement append data flow**

Remove concatenation of graph-package cDNA as the history alignment query. Retain graph-package gDNA/BED/BAM/edges. Use the explicit history cDNA for history-to-query, combine history plus query cDNA for transcript lengths/rescue, and forward the historical GTF only to final last-cluster output.

- [ ] **Step 7: Run all append unit tests**

Run:

```bash
/data/changq/miniforge3/bin/python -m pytest tests/test_append_construct_flow.py tests/test_incremental_graph_append.py tests/test_historical_gtf_transcripts.py -q
```

Expected: all tests pass.

### Task 4: Documentation and real-validation command update

**Files:**
- Modify: `README.md`
- Modify: `validation/run_incremental_graph_append_JM22.slurm`
- Modify: `tests/test_validation_scripts.py`
- Modify: `pyproject.toml`
- Modify: `src/pantrans/version.py`

- [ ] **Step 1: Write failing validation-script assertions**

Assert the SLURM command contains `--query-cdna`, `--history-gtf`, and an unrenamed historical `--cdna` input. Keep assertions for graph package, combined BED, comparison outputs, resources, and explicit thresholds.

- [ ] **Step 2: Run validation-script tests and verify RED**

Expected: missing new options.

- [ ] **Step 3: Update documentation and scripts**

Document both construct final output variants, the new append inputs, chaining behavior, validation failures, and the rule that historical GTF coordinates are reused only for unchanged representative IDs. Update the SLURM script with variables for the historical unrenamed GTF/cDNA and new JM22 cDNA; do not submit it. Bump the feature version consistently.

- [ ] **Step 4: Run documentation/CLI/script tests**

Run:

```bash
/data/changq/miniforge3/bin/python -m pytest tests/test_validation_scripts.py tests/test_append_construct_flow.py -q
```

Expected: all tests pass.

### Task 5: Full verification and review

**Files:**
- Review all modified source, tests, docs, and validation files.

- [ ] **Step 1: Run the complete test suite**

```bash
/data/changq/miniforge3/bin/python -m pytest -q
```

Expected: zero failures.

- [ ] **Step 2: Run syntax and whitespace checks**

```bash
/data/changq/miniforge3/bin/python -m compileall -q src scripts tests
git diff --check
bash -n validation/run_incremental_graph_append_JM22.slurm
```

Expected: all commands exit 0.

- [ ] **Step 3: Run a small real-file end-to-end test**

Use temporary miniature FASTA, BED, GTF, and BAM fixtures to execute construct final-output generation followed by append final-output generation without mocking transcript parsing. Verify official and unrenamed GTF/cDNA pairs exist, transcript IDs match their GTFs, and historical seed models survive when absent from the append BAM.

- [ ] **Step 4: Review the final diff against the design**

Confirm no graph clustering behavior changed, no full historical cDNA realignment was reintroduced, CLI help matches implementation, and all generated-path names are documented.
