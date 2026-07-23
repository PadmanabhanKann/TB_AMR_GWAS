#!/usr/bin/env python3
import csv
import sys
import numpy as np
from pathlib import Path

IN  = Path("/data/biol-micro-genomics/pkannan/TB/treewas/kmer_matrix")
OUT = Path("/data/biol-micro-genomics/pkannan/TB/treewas/kmer_matrix_T.csv")
MM  = Path("/data/biol-micro-genomics/pkannan/TB/treewas/kmer_matrix_T.uint8.memmap")
KM  = Path("/data/biol-micro-genomics/pkannan/TB/treewas/kmer_ids.txt")

print("[1/4] Reading header (sample IDs) + counting kmers...", file=sys.stderr)
with open(IN, "r") as f:
    reader = csv.reader(f)
    header = next(reader)
    samples = header[1:]                 # REAL sample IDs (3080)
    n_samples = len(samples)

    n_kmers = 0
    for _ in reader:
        n_kmers += 1

print(f"  samples={n_samples:,}", file=sys.stderr)
print(f"  kmers  ={n_kmers:,}", file=sys.stderr)

print("[2/4] Creating disk-backed matrix X (samples × kmers)...", file=sys.stderr)
X = np.memmap(MM, dtype=np.uint8, mode="w+", shape=(n_samples, n_kmers))
X[:] = 0
X.flush()

print("[3/4] Streaming input: fill X and write kmer IDs...", file=sys.stderr)
with open(IN, "r") as f, open(KM, "w") as km_out:
    reader = csv.reader(f)
    header = next(reader)

    for kmer_j, row in enumerate(reader):
        if (kmer_j + 1) % 10000 == 0:
            print(f"  processed kmers: {kmer_j+1:,}/{n_kmers:,}", file=sys.stderr)

        kmer_id = row[0]
        km_out.write(kmer_id + "\n")

        vals = row[1:]
        # fill column kmer_j across samples
        for sample_i, v in enumerate(vals):
            if v == "1":
                X[sample_i, kmer_j] = 1

X.flush()

print("[4/4] Writing transposed CSV (rows=samples, cols=kmers)...", file=sys.stderr)
# Read kmer IDs for header
with open(KM, "r") as km_in:
    kmer_ids = [line.rstrip("\n") for line in km_in]

with open(OUT, "w", newline="") as out:
    writer = csv.writer(out)
    writer.writerow(["sample"] + kmer_ids)

    for i, s in enumerate(samples):
        writer.writerow([s] + X[i, :].tolist())

print("[DONE] Wrote:", OUT, file=sys.stderr)
print("[DONE] Sanity check expected: header NF =", 1 + n_kmers, file=sys.stderr)
