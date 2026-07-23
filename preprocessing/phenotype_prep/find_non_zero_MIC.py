#!/usr/bin/env python3
import argparse
import pandas as pd
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(
        description="Extract Genome IDs with non-zero Measurement Value."
    )
    parser.add_argument("-i", "--input", required=True, help="Input MIC .csv file")
    parser.add_argument("-o", "--output", help="Optional output CSV (Genome ID + MIC)")
    parser.add_argument("--id-column", default="Genome ID", help="Column name for Genome IDs (default: Genome ID)")
    parser.add_argument("--mic-column", default="Measurement Value", help="Column name for MIC values (default: Measurement Value)")
    args = parser.parse_args()

    # Load CSV as strings
    df = pd.read_csv(args.input, dtype=str)

    if args.id_column not in df.columns or args.mic_column not in df.columns:
        raise ValueError(f"Columns not found. Available: {list(df.columns)}")

    total_ids = df[args.id_column].notna().sum()

    # Convert MIC to numeric
    df[args.mic_column] = pd.to_numeric(df[args.mic_column], errors="coerce")

    # Filter MIC > 0
    filtered = df[df[args.mic_column] > 0][[args.id_column, args.mic_column]].dropna()

    # Replace '.' with '_' in IDs
    filtered[args.id_column] = filtered[args.id_column].astype(str).str.replace(".", "_", regex=False)

    # Report counts
    print(f"Total Genome IDs in file: {total_ids}")
    print(f"Number of Genome IDs with non-zero {args.mic_column}: {len(filtered)}")

    # Only write file if requested
    if args.output:
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        filtered.to_csv(args.output, index=False)
        print(f"Saved IDs + MICs to {args.output}")

if __name__ == "__main__":
    main()
