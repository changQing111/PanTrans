"""Ordered, bounded parallel execution of the original TransDecoder scorer.

This module is also invoked directly as a script by a private Predict adapter.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from time import perf_counter

HEADER = b"#acc\tMarkov_order\tseq_length\tscore_1\tscore_2\tscore_3\tscore_4\tscore_5\tscore_6\n"


def records(path):
    """Stream FASTA records, retaining their original bytes and order."""
    record = []
    with open(path, "rb") as handle:
        for line in handle:
            if line.startswith(b">"):
                if record:
                    yield record
                record = [line]
            elif record:
                record.append(line)
            elif line.strip():
                raise ValueError("Sequence before FASTA header")
    if record:
        yield record


def score_parallel(cds, model, scorer, workers, output):
    if workers < 1:
        raise ValueError("Scoring workers must be positive")
    started = perf_counter()
    total = sum(sum(len(line.strip()) for line in record[1:]) for record in records(cds))
    # Consecutive chunks permit a simple ordered merge without sorting IDs.
    target = max(1, (total + workers - 1) // workers)
    root = Path(tempfile.mkdtemp(prefix="parallel-scores-", dir=str(Path(cds).resolve().parent)))
    chunks = []
    handle = None
    bases = 0
    try:
        for record in records(cds):
            if handle is None or (bases >= target and len(chunks) < workers):
                if handle:
                    handle.close()
                chunk = root / ("chunk-%04d.fa" % len(chunks))
                chunks.append(chunk)
                handle = chunk.open("wb")
                bases = 0
            handle.writelines(record)
            # Ensure a terminal newline even when the original last record lacks one.
            if not record[-1].endswith(b"\n"):
                handle.write(b"\n")
            bases += sum(len(line.strip()) for line in record[1:])
    finally:
        if handle:
            handle.close()

    def run(chunk):
        scores = chunk.with_suffix(".scores")
        with scores.open("wb") as dest, chunk.with_suffix(".stderr").open("wb") as err:
            subprocess.run([str(scorer), str(chunk), str(model)], stdout=dest, stderr=err, check=True)
        # Validate coverage and order before publishing any score data.
        with scores.open("rb") as source:
            if source.readline() != HEADER:
                raise ValueError("Unexpected TransDecoder score header: " + str(scores))
            for record in records(chunk):
                fields = source.readline().rstrip(b"\r\n").split(b"\t")
                if len(fields) != 9 or fields[0] != record[0][1:].split()[0]:
                    raise ValueError("Missing or reordered TransDecoder score: " + str(scores))
            if source.read():
                raise ValueError("Extra TransDecoder scores: " + str(scores))
        return scores

    print("Parallel CDS scoring: %d workers, %d chunks; diagnostics: %s" %
          (min(workers, len(chunks)), len(chunks), root), file=sys.stderr, flush=True)
    with ThreadPoolExecutor(max_workers=workers) as pool:
        scores = list(pool.map(run, chunks))
    output.write(HEADER)
    for path in scores:
        with path.open("rb") as source:
            source.readline()
            shutil.copyfileobj(source, output)
    output.flush()
    print("Parallel CDS scoring finished in %.3f s" % (perf_counter() - started), file=sys.stderr, flush=True)
    # Retain small diagnostics, but discard duplicate sequences/scores on success.
    for chunk, score in zip(chunks, scores):
        chunk.unlink()
        score.unlink()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", required=True, type=int)
    parser.add_argument("--scorer", required=True)
    parser.add_argument("cds")
    parser.add_argument("model")
    args = parser.parse_args()
    score_parallel(args.cds, args.model, args.scorer, args.workers, sys.stdout.buffer)


if __name__ == "__main__":
    main()
