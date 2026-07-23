#!/usr/bin/env python3
import os
import pandas as pd
from collections import defaultdict

# Run from inside pyseer2/

ANTIBIOTICS = [
    "amikacin",
    "bedaquiline",
    "clofazimine",
    "delamanid",
    "ethambutol",
    "ethionamide",
    "isoniazid",
    "kanamycin",
    "levofloxacin",
    "linezolid",
    "moxifloxacin",
    "rifabutin",
    "rifampicin",
]

OUTFILE = "gene_occurrence_across_antibiotics.csv"

gene_to_antibiotics = defaultdict(set)

for ab in ANTIBIOTICS:
    infile = os.path.join(ab, "gene_hits.sorted.txt")

    if not os.path.isfile(infile):
        print(f"[SKIP] Missing: {infile}")
        continue

    df = pd.read_csv(infile, sep="\t")

    if "gene" not in df.columns:
        print(f"[SKIP] {ab}: no 'gene' column")
        continue

    # clean gene names
    df = df[df["gene"].notna()].copy()
    df["gene"] = df["gene"].astype(str).str.strip()
    df = df[df["gene"] != ""]

    # unique genes per antibiotic
    unique_genes = set(df["gene"])

    for g in unique_genes:
        gene_to_antibiotics[g].add(ab)

    print(f"[OK] {ab}: {len(unique_genes)} genes")

# build output
rows = []
for gene, abs_set in gene_to_antibiotics.items():
    abs_sorted = [ab for ab in ANTIBIOTICS if ab in abs_set]
    rows.append({
        "gene": gene,
        "n_antibiotics": len(abs_sorted),
        "antibiotics": ",".join(abs_sorted)
    })

# write file
if rows:
    df_out = pd.DataFrame(rows).sort_values(
        by=["n_antibiotics", "gene"],
        ascending=[False, True]
    )
    df_out.to_csv(OUTFILE, index=False)
    print(f"[OK] Wrote {OUTFILE}")
else:
    print("No genes found.")
