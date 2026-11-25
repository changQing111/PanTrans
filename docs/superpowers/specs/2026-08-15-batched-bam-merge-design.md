# Batched BAM Merge Design

## Problem

The rescue workflow can produce thousands of per-cluster BAM files. The current
`_merge_bam_files()` passes every input to one `samtools merge` process, which
requires too many simultaneously open files and fails under common `nofile`
limits.

## Design

Keep the existing `_merge_bam_files(bam_paths, merged_bam_path)` interface and
its current empty-list and single-input behavior. For multiple inputs, merge in
rounds with at most 128 input BAMs per `samtools merge` invocation. Each group
produces an intermediate BAM; subsequent rounds merge those intermediates until
one file remains, which becomes `merged_bam_path`.

Intermediate files live in a temporary directory next to the requested output
so the final move remains on the same filesystem. The temporary directory is
removed on success and on any exception. A failed `samtools` invocation remains
visible to callers through `subprocess.CalledProcessError`.

## Compatibility

- Empty input returns `False` without creating an output.
- One input is copied to the output and returns `True`.
- Multiple inputs return `True` after a successful hierarchical merge.
- Existing callers and output names remain unchanged.
- Every `samtools merge` command receives no more than 128 input BAMs.

## Tests

Add focused unit tests using a controlled `subprocess.run` replacement because
the behavior under test is command partitioning, while actual BAM merge behavior
belongs to samtools. Verify that 129 inputs require two rounds, no invocation
exceeds 128 inputs, the requested output is produced, and temporary files are
removed. Also verify cleanup when a merge invocation fails.
