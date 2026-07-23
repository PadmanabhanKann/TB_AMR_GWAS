#!/usr/bin/env python3
import csv, sys, os

def die(msg):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)

if len(sys.argv) != 5:
    die("Usage: subset_kmer_matrix.py <kmer_matrix_T.csv> <kmers.list> <out.csv> <id_col_name>\n"
        "Example: subset_kmer_matrix.py kmer_matrix_T.csv ethambutol.kmers.list ethambutol_subset.csv sample")

matrix_path, kmers_list_path, out_path, id_col = sys.argv[1:]

if not os.path.exists(matrix_path): die(f"Matrix not found: {matrix_path}")
if not os.path.exists(kmers_list_path): die(f"Kmer list not found: {kmers_list_path}")

# Load desired k-mers (one per line)
wanted = []
wanted_set = set()
with open(kmers_list_path, "r") as f:
    for line in f:
        k = line.strip()
        if not k: 
            continue
        if k not in wanted_set:
            wanted.append(k)
            wanted_set.add(k)

if not wanted:
    die("Kmer list is empty.")

with open(matrix_path, newline="") as fin:
    reader = csv.reader(fin)
    try:
        header = next(reader)
    except StopIteration:
        die("Matrix file is empty.")

    # Map header -> index
    idx = {name: i for i, name in enumerate(header)}

    if id_col not in idx:
        die(f"ID column '{id_col}' not found in matrix header. First few headers: {header[:5]}")

    keep_indices = [idx[id_col]]
    missing = []
    for k in wanted:
        if k in idx:
            keep_indices.append(idx[k])
        else:
            missing.append(k)

    # Write output
    with open(out_path, "w", newline="") as fout:
        writer = csv.writer(fout)
        writer.writerow([header[i] for i in keep_indices])

        for row in reader:
            # Some rows could be shorter; guard
            outrow = []
            for i in keep_indices:
                outrow.append(row[i] if i < len(row) else "")
            writer.writerow(outrow)

print(f"Wrote: {out_path}", file=sys.stderr)
print(f"Requested kmers: {len(wanted)}", file=sys.stderr)
print(f"Found in matrix: {len(keep_indices)-1}", file=sys.stderr)
print(f"Missing kmers: {len(missing)}", file=sys.stderr)
if missing:
    # Print first few missing to stderr (avoid huge spam)
    print("First missing examples:", ", ".join(missing[:10]), file=sys.stderr)
