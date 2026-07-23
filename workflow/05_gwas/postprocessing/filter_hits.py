#!/usr/bin/env python3
import os
import pandas as pd
from collections import defaultdict

# Run from inside pyseer2/

MAXP_THRESH = {
    "amikacin": 20,
    "bedaquiline": 15,
    "clofazimine": 10,
    "delamanid": 15,
    "ethambutol": 20,
    "ethionamide": 15,
    "isoniazid": 20,
    "kanamycin": 15,
    "levofloxacin": 20,
    "linezolid": 10,
    "moxifloxacin": 20,
    "rifabutin": 15,
    "rifampicin": 20,
}

MAF_MAX = 0.2  # keep avg_maf <= 0.2

OUT_SUMMARY_TSV = "gene_hits_filter_summary.tsv"
OUT_COMMON_CSV  = "filtered_gene_hits_most_common.csv"

summary_rows = []

# Cross-antibiotic tracking
gene_to_antibiotics = defaultdict(set)

# Per-gene, per-antibiotic values (only for genes that PASS filters in that antibiotic)
# store: gene -> { ab -> {"maxp": val, "avg_beta": val, "hits": val, "avg_maf": val} }
gene_ab_values = defaultdict(dict)

for ab, thr in MAXP_THRESH.items():
    infile = os.path.join(ab, "gene_hits.sorted.txt")
    if not os.path.isfile(infile):
        print(f"[SKIP] Missing: {infile}")
        continue

    df = pd.read_csv(infile, sep="\t")

    required = {"gene", "maxp", "avg_maf"}
    missing = required - set(df.columns)
    if missing:
        print(f"[SKIP] {ab}: missing columns {missing} in {infile}")
        continue

    # numeric cleanup
    df["maxp"] = pd.to_numeric(df["maxp"], errors="coerce")
    df["avg_maf"] = pd.to_numeric(df["avg_maf"], errors="coerce")
    if "avg_beta" in df.columns:
        df["avg_beta"] = pd.to_numeric(df["avg_beta"], errors="coerce")
    if "hits" in df.columns:
        df["hits"] = pd.to_numeric(df["hits"], errors="coerce")

    before = len(df)

    # Apply filters
    filt = df[(df["maxp"] >= thr) & (df["avg_maf"] >= MAF_MAX)].copy()
    after = len(filt)

    outfile = os.path.join(ab, f"gene_hits.filtered_maxp{thr}_maf{MAF_MAX}.txt")
    filt.to_csv(outfile, sep="\t", index=False)

    summary_rows.append({
        "antibiotic": ab,
        "maxp_threshold": thr,
        "maf_max": MAF_MAX,
        "n_total_genes": before,
        "n_pass": after,
        "outfile": outfile
    })

    print(f"[OK] {ab}: {after}/{before} genes passed → {outfile}")

    # Collect for "most common" and per-antibiotic columns
    for _, row in filt.iterrows():
        g = str(row["gene"])
        gene_to_antibiotics[g].add(ab)

        gene_ab_values[g][ab] = {
            "maxp": row.get("maxp", None),
            "avg_beta": row.get("avg_beta", None),
            "hits": row.get("hits", None),
            "avg_maf": row.get("avg_maf", None),
        }

# Write filtering summary
if summary_rows:
    sdf = pd.DataFrame(summary_rows).sort_values(
        ["n_pass", "antibiotic"], ascending=[False, True]
    )
    sdf.to_csv(OUT_SUMMARY_TSV, sep="\t", index=False)
    print(f"[OK] Wrote {OUT_SUMMARY_TSV}")
else:
    print("No filtering outputs generated (no inputs found or all skipped).")

# Build "most common genes" table
antibiotics_order = list(MAXP_THRESH.keys())  # keep your defined order

common_rows = []
for gene, abs_set in gene_to_antibiotics.items():
    abs_sorted = [ab for ab in antibiotics_order if ab in abs_set]
    n_ab = len(abs_sorted)

    # aggregates across antibiotics where gene passed
    maxp_vals = []
    beta_vals = []
    hits_vals = []
    maf_vals  = []

    for ab in abs_sorted:
        v = gene_ab_values[gene].get(ab, {})
        if v.get("maxp") is not None and pd.notna(v.get("maxp")):
            maxp_vals.append(float(v["maxp"]))
        if v.get("avg_beta") is not None and pd.notna(v.get("avg_beta")):
            beta_vals.append(float(v["avg_beta"]))
        if v.get("hits") is not None and pd.notna(v.get("hits")):
            hits_vals.append(float(v["hits"]))
        if v.get("avg_maf") is not None and pd.notna(v.get("avg_maf")):
            maf_vals.append(float(v["avg_maf"]))

    row_out = {
        "gene": gene,
        "n_antibiotics": n_ab,
        "antibiotics": ",".join(abs_sorted),
        "hits_sum_across_antibiotics": (sum(hits_vals) if hits_vals else None),
        "maxp_mean_across_antibiotics": (sum(maxp_vals)/len(maxp_vals) if maxp_vals else None),
        "avg_beta_mean_across_antibiotics": (sum(beta_vals)/len(beta_vals) if beta_vals else None),
        "avg_maf_mean_across_antibiotics": (sum(maf_vals)/len(maf_vals) if maf_vals else None),
    }

    # add per-antibiotic columns: <ab>_maxp and <ab>_avg_beta
    for ab in antibiotics_order:
        v = gene_ab_values[gene].get(ab)
        row_out[f"{ab}_maxp"] = (float(v["maxp"]) if v and pd.notna(v.get("maxp")) else "")
        row_out[f"{ab}_avg_beta"] = (float(v["avg_beta"]) if v and pd.notna(v.get("avg_beta")) else "")

    common_rows.append(row_out)

# Write most-common CSV sorted by n_antibiotics then mean maxp
if common_rows:
    cdf = pd.DataFrame(common_rows)

    # ensure numeric sort works even with blanks
    cdf["n_antibiotics"] = pd.to_numeric(cdf["n_antibiotics"], errors="coerce")
    cdf["maxp_mean_across_antibiotics"] = pd.to_numeric(cdf["maxp_mean_across_antibiotics"], errors="coerce")

    cdf = cdf.sort_values(
        by=["n_antibiotics", "maxp_mean_across_antibiotics", "gene"],
        ascending=[False, False, True],
        na_position="last"
    )

    cdf.to_csv(OUT_COMMON_CSV, index=False)
    print(f"[OK] Wrote {OUT_COMMON_CSV} (rows={len(cdf)})")
else:
    print("No common-gene table generated (no genes passed filters).")
