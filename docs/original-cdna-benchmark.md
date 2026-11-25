# Full original-cDNA append benchmark (2026-11-02)

`pantrans append --history-graph-cdna original.fa` generates old-to-new graph
edges from the supplied original historical cDNA. The FASTA must cover exactly
the historical graph gene set. `--cdna` and `--history-gff3` retain their meaning:
final unrenamed historical annotation inputs. The default behavior is unchanged.

The new alignment is graph-only, saved as `*_history_graph_to_query[.filtered].bam`.
The reconstructed-cDNA alignment is still computed and used for annotation.
This prevents same-named original and reconstructed transcripts with different
sequences from being confused. The merged annotation BAM therefore does not
necessarily encode the graph's complete edge set; use the saved graph edges to
reproduce graph clustering. Both query sources are recorded in graph provenance.

This first full benchmark measures that conservative separation, including the
additional alignment cost. It is not a benchmark of simply replacing the input
of one alignment call, and it does not measure a future optimized implementation.

On PKU-IAAS, `/data/changq/Project/pantrans_test/full_benchmark_20261102/` contains
a frozen source checkout and orchestration scripts. All runs use 32 CPUs, 240G,
rescue, CDS annotation, and full alignments without reused BAMs. Construct5 and
construct10 run independently; after construct5 succeeds, a dependent launcher
prepares/submits append using its new final annotation and graph plus the original
five-variety cDNA for graph edges. Every analysis has a resource collector.
The summary waits for collectors and produces `comparison.json` with separate
total construct5+append and append-only time ratios against construct10.

Jobs initially submitted: construct5 18695301; construct10 18695304;
append preparation/submission launcher 18695306. The append job ID will be stored
in the run manifest at the directory recorded in `append.path` after dispatch.
Runs are recorded in the shared experiment_runs/runs.tsv.

The pipelines run on allocated cluster nodes; node/hardware and scheduling
differences should be considered when interpreting a single timing replicate.
RSS values retain the runner's documented scope, not simultaneous aggregate
process-tree peak. Initial validation: 52 tests and 804 subtests passed.
