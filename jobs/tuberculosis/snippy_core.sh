#!/bin/bash
#SBATCH --job-name=snippy_core
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --partition=short,medium,long
#SBATCH --output=/data/biol-micro-genomics/pkannan/logs/snippy_core_%j.out
#SBATCH --error=/data/biol-micro-genomics/pkannan/logs/snippy_core_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=padmanabhan21@iisertvm.ac.in

set -euo pipefail

module purge
module load snippy/4.6.0-foss-2021b-R-4.1.2

REF=/data/biol-micro-genomics/pkannan/TB/epiML/snippy/ref/GCF_000195955.2_ASM19595v2_genomic.fna
WORK=/data/biol-micro-genomics/pkannan/TB/epiML/snippy
OUT=$WORK/results

cd "$WORK"

# sanity check: ensure snippy outputs exist
n=$(find "$OUT" -type f -name "snps.tab" | wc -l)
echo "Found $n per-sample snps.tab files"
if [[ "$n" -lt 1 ]]; then
  echo "ERROR: No snps.tab found in $OUT (per-sample snippy outputs missing?)"
  exit 1
fi

# Build core outputs (written into $WORK)
snippy-core --ref "$REF" "$OUT"/*

echo "snippy-core finished successfully."
echo "Final outputs in: $WORK"
ls -lh core.*
