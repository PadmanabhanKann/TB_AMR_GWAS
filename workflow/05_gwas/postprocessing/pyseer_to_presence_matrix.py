#!/usr/bin/env python3
import argparse
import gzip
import sys


def convert_pyseer_to_dense_csv(pyseer_file, output_csv, progress_every=100000):
    """
    Convert pyseer sparse presence format to a dense presence/absence CSV matrix
    WITHOUT holding the full matrix in memory.

    Rows: unitigs (sequences)
    Cols: samples
    Values: 0/1
    """

    # -------------------------
    # PASS 1: collect all samples
    # -------------------------
    print("[pass1] Collecting sample IDs...", file=sys.stderr)
    all_samples = set()
    n_unitigs = 0

    with gzip.open(pyseer_file, "rt") as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue

            parts = line.split(" | ")
            if len(parts) != 2:
                continue

            _, samples_str = parts
            for token in samples_str.split():
                all_samples.add(token.split(":")[0])

            n_unitigs += 1
            if progress_every and (line_num % progress_every == 0):
                print(f"[pass1] lines={line_num:,} unitigs={n_unitigs:,} samples={len(all_samples):,}", file=sys.stderr)

    if n_unitigs == 0:
        raise RuntimeError("No valid pyseer lines found. Check the input format.")

    samples = sorted(all_samples)
    ncols = len(samples)
    col_index = {s: i for i, s in enumerate(samples)}  # 0-based

    print(f"[pass1] Done: unitigs={n_unitigs:,}, samples={ncols:,}", file=sys.stderr)

    # -------------------------
    # PASS 2: stream rows -> CSV
    # -------------------------
    print("[pass2] Writing dense CSV (streaming)...", file=sys.stderr)

    # Output can be .csv or .csv.gz (auto based on extension)
    if output_csv.endswith(".gz"):
        out_handle = gzip.open(output_csv, "wt")
    else:
        out_handle = open(output_csv, "wt")

    with out_handle as out, gzip.open(pyseer_file, "rt") as f:
        # Write header
        out.write("unitig," + ",".join(samples) + "\n")

        written = 0
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue

            parts = line.split(" | ")
            if len(parts) != 2:
                continue

            unitig_seq, samples_str = parts

            # Initialize row as all 0s (strings for fast join)
            row = ["0"] * ncols

            # Set 1s where present
            for token in samples_str.split():
                s = token.split(":")[0]
                j = col_index.get(s)
                if j is not None:
                    row[j] = "1"

            # Write one CSV row (unitig + dense vector)
            out.write(unitig_seq + "," + ",".join(row) + "\n")

            written += 1
            if progress_every and (written % max(1, progress_every // 10) == 0):
                # More meaningful progress: based on rows written
                print(f"[pass2] wrote {written:,} / {n_unitigs:,} rows", file=sys.stderr)

    print(f"[done] Wrote dense matrix CSV: {output_csv}", file=sys.stderr)
    print(f"[done] Shape: {n_unitigs:,} unitigs × {ncols:,} samples", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(
        description="Convert pyseer presence (.pyseer.gz) to a dense 0/1 CSV matrix (streaming, memory-efficient)."
    )
    ap.add_argument(
        "-i", "--input",
        required=True,
        help="Input pyseer gzipped file (e.g. sigkmer_presence.ordered.pyseer.gz)"
    )
    ap.add_argument(
        "-o", "--output",
        required=True,
        help="Output CSV file path (.csv or .csv.gz)"
    )
    ap.add_argument(
        "--progress-every",
        type=int,
        default=100000,
        help="Progress reporting frequency (lines/rows). Default: 100000. Use 0 to disable."
    )
    args = ap.parse_args()

    convert_pyseer_to_dense_csv(
        pyseer_file=args.input,
        output_csv=args.output,
        progress_every=args.progress_every
    )


if __name__ == "__main__":
    main()
