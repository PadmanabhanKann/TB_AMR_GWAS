#!/usr/bin/env python3
import sys, math
import pandas as pd

if len(sys.argv) < 4:
    print("Usage: python kmer_filter.py in.csv out.csv min_prev [max_prev]", file=sys.stderr)
    sys.exit(1)

IN, OUT = sys.argv[1], sys.argv[2]
minp = float(sys.argv[3])
maxp = float(sys.argv[4]) if len(sys.argv) > 4 else 0.99

COL_CHUNK = 20000
ROW_CHUNK = 50000

# header
cols = pd.read_csv(IN, nrows=0).columns.tolist()
idcol, kmers = cols[0], cols[1:]
K = len(kmers)
if K == 0:
    raise SystemExit("No k-mer columns")

# N (rows excluding header)
with open(IN, "rb") as f:
    N = sum(1 for _ in f) - 1
if N <= 0:
    raise SystemExit("No rows")

minc = math.ceil(minp * N)
maxc = math.floor(maxp * N)
print(f"N={N} K={K} keep count in [{minc},{maxc}]")

kept = []
for i in range(0, K, COL_CHUNK):
    chunk = kmers[i:i+COL_CHUNK]
    print(f"sum cols {i+1}-{i+len(chunk)}"); sys.stdout.flush()
    df = pd.read_csv(IN, usecols=[idcol] + chunk, dtype={c: "uint8" for c in chunk})
    s = df[chunk].sum(axis=0)
    kept += [c for c in chunk if (minc <= int(s[c]) <= maxc)]

print(f"kept {len(kept)} / {K}")
if not kept:
    raise SystemExit("No columns kept")

written, first = 0, True
for df in pd.read_csv(IN, usecols=[idcol] + kept, chunksize=ROW_CHUNK):
    df.to_csv(OUT, mode="w" if first else "a", index=False, header=first)
    first = False
    written += len(df)
    if written % (ROW_CHUNK * 2) == 0:
        print(f"wrote {written}/{N}"); sys.stdout.flush()

print(f"done: wrote {written}/{N}")
if written != N:
    raise SystemExit("ERROR: output rows != input rows")
