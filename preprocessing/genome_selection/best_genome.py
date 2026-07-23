#!/usr/bin/env python3
import pandas as pd
import sys
import os
import argparse

parser = argparse.ArgumentParser(
    description="Filter assemblies by QC thresholds and report pass/fail counts."
)
parser.add_argument("-i", "--input", required=True, help="Path to input TSV file")
parser.add_argument("-o", "--output", help="Path to output file with passing genome IDs (optional)")

# QC thresholds
parser.add_argument("--max-contigs", type=int, default=1000,
                    help="Exclude assemblies with more than this many contigs (default: 1000)")
parser.add_argument("--min-n50", type=int, default=10_000,
                    help="Exclude assemblies with N50 below this (bp). Default: 10000")
parser.add_argument("--expected-size", type=float, default=None,
                    help="Expected genome size in bp (e.g., 2.2e6). If set, apply tolerance filter.")
parser.add_argument("--tolerance-pct", type=float, default=15.0,
                    help="Allowed ±%% around expected size (default: 15%%)")

args = parser.parse_args()

if not os.path.exists(args.input):
    print(f"Error: Input file '{args.input}' does not exist.", file=sys.stderr)
    sys.exit(1)

# Load
df = pd.read_csv(args.input, sep="\t", dtype={"Assembly": str})

# Required columns
req_cols = ['Assembly', 'N50', 'L50', '# contigs', "# N's per 100 kbp", 'Total length']
missing = [c for c in req_cols if c not in df.columns]
if missing:
    print(f"Error: Missing required columns: {missing}", file=sys.stderr)
    sys.exit(1)

# Cast numeric
num_cols = [c for c in req_cols if c != 'Assembly']
df[num_cols] = df[num_cols].apply(pd.to_numeric, errors='coerce')

# Add status column
df['Status'] = 'PASS'

# Apply filters step by step
df.loc[df['# contigs'] > args.max_contigs, 'Status'] = 'FAIL_contigs'
df.loc[df['N50'] < args.min_n50, 'Status'] = 'FAIL_N50'

if args.expected_size is not None:
    tol_frac = args.tolerance_pct / 100.0
    lower = args.expected_size * (1 - tol_frac)
    upper = args.expected_size * (1 + tol_frac)
    df.loc[(df['Total length'] < lower) | (df['Total length'] > upper),
           'Status'] = 'FAIL_size'

# Passed assemblies
passed = df[df['Status'] == 'PASS']

# Write only if output is specified
if args.output:
    passed['Assembly'].to_csv(args.output, index=False, header=False)
    print(f"Passing assembly IDs written to: {args.output}")

# Print summary
total = len(df)
n_pass = len(passed)
n_fail = total - n_pass
reason_counts = df[df['Status'] != 'PASS']['Status'].value_counts()

print(f"✅ QC filtering complete.")
print(f"  Total assemblies: {total}")
print(f"  Passed: {n_pass}")
print(f"  Failed: {n_fail}")
if n_fail > 0:
    print("  Fail reasons:")
    for reason, count in reason_counts.items():
        print(f"    {reason}: {count}")
