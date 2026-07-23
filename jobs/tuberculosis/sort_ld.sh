#!/bin/bash
#SBATCH --job-name=sort_ld_r2
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --partition=short
#SBATCH --time=10:00:00
#SBATCH --output=/data/biol-micro-genomics/pkannan/logs/sort_ld_r2_%j.out
#SBATCH --error=/data/biol-micro-genomics/pkannan/logs/sort_ld_r2_%j.err

set -euo pipefail

# ---- paths ----
WORKDIR=/data/biol-micro-genomics/pkannan/TB/plink
INFILE=mtb_ld.ld
OUTFILE=mtb_ld_sorted_maxR2.ld

cd "$WORKDIR"

echo "Sorting LD file by max R2..."
echo "Input : $INFILE"
echo "Output: $OUTFILE"

# Preserve header, sort numerically by R2 (column 7, descending)
awk 'NR==1{print; next} {print | "sort -k7,7gr"}' "$INFILE" > "$OUTFILE"

echo "Done."
echo "Top 5 rows:"
head -n 6 "$OUTFILE"
