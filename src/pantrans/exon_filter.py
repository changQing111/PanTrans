"""Validated preprocessing coordinates for masked mapping and exon-aware edges."""
import csv
import hashlib
import json
from collections import Counter

from Bio import SeqIO

LOGIC_ID = "pantrans-exon-filter-v1"


class ExonPolicy:
    def __init__(self, gene_tables, exon_tables, overlap=0.90, query_overlap=0.80):
        if not 0 < overlap <= 1 or not 0 < query_overlap <= 1:
            raise ValueError("Exon thresholds must be in (0, 1]")
        self.overlap, self.query_overlap = overlap, query_overlap
        self.genes, self.transcripts = {}, {}
        self.files = []
        for paths, kind in ((gene_tables, "gene"), (exon_tables, "exon")):
            for path in paths:
                with open(path, "rb") as handle:
                    digest = hashlib.sha256(handle.read()).hexdigest()
                self.files.append(dict(kind=kind, path=str(path), sha256=digest))
                with open(path, encoding="utf-8-sig") as handle:
                    for row in csv.DictReader(handle, delimiter="\t"):
                        gene = row["gene"]
                        if kind == "gene":
                            length, a, b = [int(row[k]) for k in ("gdna_length", "body_start0", "body_end0")]
                            if gene in self.genes or not 0 <= a < b <= length:
                                raise ValueError(f"Duplicate/invalid gene region: {gene}")
                            self.genes[gene] = (length, a, b)
                        else:
                            tid = row["transcript"]
                            a, b = int(row["gdna_start0"]), int(row["gdna_end0"])
                            strand = row["gdna_strand"]
                            if gene not in self.genes or strand not in ("+", "-"):
                                raise ValueError(f"Unknown gene/strand in exon table: {gene}")
                            _, left, right = self.genes[gene]
                            if not left <= a < b <= right:
                                raise ValueError(f"Exon outside gene body: {tid}")
                            self.transcripts.setdefault(gene, {}).setdefault((tid, strand), set()).add((a, b))
        if not self.genes or set(self.genes) != set(self.transcripts):
            raise ValueError("Every gene region must have target transcript exons")
        for gene, transcripts in self.transcripts.items():
            for key, intervals in list(transcripts.items()):
                intervals = sorted(intervals)
                if any(x[1] > y[0] for x, y in zip(intervals, intervals[1:])):
                    raise ValueError(f"Overlapping exons: {key[0]}")
                transcripts[key] = intervals
        self.digests = {g: hashlib.sha256(json.dumps([self.genes[g],
                       sorted((t, strand, e) for (t, strand), e in self.transcripts[g].items())],
                       sort_keys=True).encode()).hexdigest() for g in self.genes}

    def provenance(self):
        return dict(filter_logic_id=LOGIC_ID, target_mask="gene-body-only-N-flanks-v1",
                    exon_thresholds=dict(target_fraction=self.overlap, query_fraction=self.query_overlap),
                    annotation_sources=self.files, gene_region_digests=self.digests)

    def validate_targets(self, lengths):
        for gene, length in lengths.items():
            if gene not in self.genes or self.genes[gene][0] != length:
                raise ValueError(f"Missing/mismatched target coordinate table: {gene}")

    def mask(self, source, destination):
        seen = set()
        with open(source) as src, open(destination, "w") as dst:
            for record in SeqIO.parse(src, "fasta"):
                if record.id in seen:
                    raise ValueError(f"Duplicate target: {record.id}")
                seen.add(record.id)
                self.validate_targets({record.id: len(record)})
                _, a, b = self.genes[record.id]
                record.seq = record.seq.__class__("N" * a + str(record.seq[a:b]) + "N" * (len(record) - b))
                SeqIO.write(record, dst, "fasta")
        if not seen:
            raise ValueError("Empty mapping target FASTA")

    def accepts(self, read, query_length):
        blocks, pos, matched = [], read.reference_start, 0
        for op, length in read.cigartuples or []:
            if op in (0, 7, 8):
                blocks.append((pos, pos + length))
                matched += length
                pos += length
            elif op in (2, 3):
                pos += length
        if not matched or not query_length:
            return False
        for (_, strand), exons in self.transcripts[read.reference_name].items():
            if read.is_reverse != (strand == "-"):
                continue
            overlap = sum(max(0, min(b, y) - max(a, x))
                          for a, b in blocks for x, y in exons)
            if overlap / matched >= self.overlap and overlap / query_length >= self.query_overlap:
                return True
        return False
