#!/usr/bin/env python3
import os
import glob
import pandas as pd
import matplotlib.pyplot as plt

BASE = "."                   # parent directory
OUTDIR = "genehits_maxp_dist"      # output folder
BINS = 40

os.makedirs(OUTDIR, exist_ok=True)

files = sorted(glob.glob(os.path.join(BASE, "*", "gene_hits.sorted.txt")))

if not files:
    raise SystemExit("No gene_hits.sorted.txt files found")

for path in files:
    antibiotic = os.path.basename(os.path.dirname(path))

    df = pd.read_csv(path, sep="\t")

    if "maxp" not in df.columns:
        print(f"[SKIP] {antibiotic}: no maxp column")
        continue

    maxp = pd.to_numeric(df["maxp"], errors="coerce").dropna()

    plt.figure()
    plt.hist(maxp, bins=BINS)
    plt.xlabel("−log10(p)  [maxp]")
    plt.ylabel("Number of genes")
    plt.title(f"{antibiotic}: gene-level GWAS signal distribution")

    out = os.path.join(OUTDIR, f"{antibiotic}_maxp_distribution.png")
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()

    print(f"[OK] {antibiotic}: {len(maxp)} genes → {out}")

print("Done.")
