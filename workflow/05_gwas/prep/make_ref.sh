#!/usr/bin/env bash
set -euo pipefail

# --- Paths ---
BASE="/data/biol-micro-genomics/pkannan/TB/pyseer"

REF_DIR="$BASE/ref"
REF_FA="$REF_DIR/ref.fna"
REF_GFF="$REF_DIR/ref.gff"

FA_DIR="/data/biol-micro-genomics/pkannan/TB/assemblies"
GFF_DIR="/data/biol-micro-genomics/pkannan/TB/panaroo/combined_gff"

OUT="$BASE/references.txt"
# -------------

# sanity checks
[[ -s "$REF_FA" ]]  || { echo "Missing $REF_FA"; exit 1; }
[[ -s "$REF_GFF" ]] || { echo "Missing $REF_GFF"; exit 1; }

# start fresh
: > "$OUT"

# reference genome (first line)
printf "%s\t%s\tref\n" "$REF_FA" "$REF_GFF" >> "$OUT"

# pair draft assemblies with GFFs
shopt -s nullglob
for fa in "$FA_DIR"/*.fas; do
    stem="$(basename "$fa" .fas)"
    gff="$GFF_DIR/$stem.gff"

    if [[ -s "$gff" ]]; then
        printf "%s\t%s\tdraft\n" "$fa" "$gff" >> "$OUT"
    else
        echo "WARNING: No matching GFF for $stem" >&2
    fi
done

echo "✔ Wrote $(wc -l < "$OUT") lines to $OUT"
