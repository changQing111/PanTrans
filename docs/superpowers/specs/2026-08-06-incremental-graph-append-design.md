# Incremental Graph Append Design

Date: 2026-08-06
Status: Approved
Branch: `codex/append-construct-flow`

## Goal

Make `append` reuse the complete directed gene graph produced by an earlier
`construct` run and calculate only relationships involving newly appended
genes. The merged graph must then pass through the same SCC,
`unit_recursion()`, chromosome assignment, and output flow as a direct
`construct` run.

## Why Representative-Only Append Is Insufficient

The five-variety graph contains 695,746 historical genes, while the current
representative-only append alignment contains 154,776 historical pre-cluster
representatives plus 136,734 JM22 genes. Merging those two graphs preserves
historical edges but cannot create JM22 relationships to historical
non-representative genes.

Comparison against the six-variety construct graph showed that 2,598,274 of
2,598,592 missing directed edges and 994,051 of 994,087 missing reciprocal
pairs involve JM22. An oracle union of the five-variety graph and all JM22
edges from the six-variety graph reached 98.26% directed-edge Jaccard and
98.45% reciprocal-pair Jaccard. This establishes that complete cross-set
edges, rather than representative-only edges, are the required increment.

## Graph Package

Every successful `construct` writes `<prefix>.graph.json` plus sidecar files:

- `<prefix>.graph.edges.tsv`: unique directed `query_gene -> target_gene`
  edges after the normal alignment filters;
- `<prefix>.graph.nodes.tsv`: every historical gene with target length and BED
  metadata;
- a JSON manifest containing the schema version, historical variety order,
  reference name, main chromosomes, sidecar paths, filtered BAM path, and
  validated source cDNA/gDNA/BED paths.

External data paths are recorded with file size and nanosecond modification
time. Loading a package fails if a required file is missing or no longer
matches its recorded identity. Sidecar files are resolved relative to the
manifest so the graph metadata can be moved as a unit.

## Append Interface

```text
pantrans append \
  --name <new-variety-names> \
  --cdna <new.cdna.fasta> \
  --gdna <new.gdna.fasta> \
  --bed <historical-pre-representatives-plus-new.bed> \
  --history-graph <historical-prefix.graph.json> \
  --threads <N> \
  --prefix <prefix> \
  --output <directory>
```

`--refer_cdna` and `--refer_gdna` are removed. The graph manifest is the
single source of truth for historical full-sequence inputs. `--bed` retains
the established contract: it contains historical pre-cluster representatives
and all new genes. The package supplies the complete historical BED needed by
the merged full graph.

## Incremental Alignment Flow

For historical set `H` and new set `N`, append performs:

1. `N cDNA -> (H + N) gDNA`, producing `N -> H` and `N -> N` edges while
   retaining full-target competition for new queries.
2. `H cDNA -> N gDNA`, producing `H -> N` edges.
3. Reuse historical `H -> H` edges from the graph package.
4. Merge and deduplicate all directed edges.

The historical `H -> H` minimap2 comparison is not repeated. With the current
695,746 historical and 136,734 new genes, this reduces the approximate
query-target pair space to 30.15% of direct six-variety construction.

The existing minimap2 `-N 100` limit means separately evaluated target sets
can retain a small number of edges that a monolithic run would displace. The
observed oracle comparison bounds this effect at roughly 1.5-1.8% of the
merged graph. The project therefore targets high similarity rather than
byte-identical graphs.

## Merged Metadata and BAM

Append builds full merged cDNA, gDNA, and BED working files from the historical
source list and new inputs. Historical BED rows take precedence for duplicate
historical representative IDs; conflicting duplicate rows are rejected. Every
new gDNA ID must be present in append `--bed`, and every non-new BED ID must
already exist in the historical node set.

The historical filtered BAM and both filtered cross-alignment BAMs are copied
into one BAM using the full merged target header. Reference IDs are remapped by
reference name, not assumed to share numeric order. This combined filtered BAM
is used by transcript processing and recorded in the next graph package.

## Shared Cluster Flow

The merged graph contains every historical and new gene. Append invokes the
same graph clustering helper as construct with:

- historical varieties followed by new varieties;
- the original construct reference name;
- the original main-chromosome configuration;
- the union of historical and new sequence-backed gene IDs as eligible nodes.

No append-specific representative preference or logical `Refer` grouping is
used. The returned pre and last clusters are written directly, and the normal
construct transcript/rescue/output functions operate on the full merged
sequence and BED files.

## Chained Append

Append writes a new graph package whose edge set and node metadata describe
the merged history plus the new varieties. Its manifest retains the historical
source sequence parts and adds the new cDNA/gDNA inputs, allowing a later
append to repeat the same incremental process without recomputing any older
pairwise block.

## Validation

Automated tests cover manifest validation, edge serialization, BED merging,
BAM reference remapping, construct package creation, append cross-alignment
orchestration, graph union, and unchanged construct clustering behavior.

Real-data validation compares the incremental merged graph with the existing
six-variety construct graph using:

- unique directed-edge precision, recall, and Jaccard;
- reciprocal-pair precision, recall, and Jaccard;
- pre/last first-column representative overlap;
- exact order-independent cluster membership matches.
