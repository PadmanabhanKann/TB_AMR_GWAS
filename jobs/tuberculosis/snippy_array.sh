#!/bin/bash
#SBATCH --job-name=snippy_array
#SBATCH --array=1-3080%80          # adjust 3080 to total assemblies; %80 = max parallel jobs
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --partition=medium,long
#SBATCH --output=/data/biol-micro-genomics/pkannan/logs/snippy_%A_%a.out
#SBATCH --error=/data/biol-micro-genomics/pkannan/logs/snippy_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=padmanabhan21@iisertvm.ac.in

set -euo pipefail

module purge
module load snippy/4.6.0-foss-2021b-R-4.1.2

# --------------------
# Paths
# --------------------
REF=/data/biol-micro-genomics/pkannan/TB/epiML/snippy/ref/GCF_000195955.2_ASM19595v2_genomic.fna
WORK=/data/biol-micro-genomics/pkannan/TB/epiML/snippy
ASSEMBLIES=/data/biol-micro-genomics/pkannan/TB/assemblies
TAB=$WORK/samples.tab
OUT=$WORK/results

mkdir -p "$OUT"

# --------------------
# Create samples.tab ONCE (race-safe)
# --------------------
if [[ ! -f "$TAB" ]]; then
  echo "Creating samples.tab from .fas assemblies"
  for f in "$ASSEMBLIES"/*.fas; do
    sample=$(basename "$f" .fas)
    echo -e "${sample}\t${f}"
  done > "$TAB"
fi

# --------------------
# Safeguard: array bounds
# --------------------
TOTAL=$(wc -l < "$TAB")
if (( SLURM_ARRAY_TASK_ID > TOTAL )); then
    echo "Array index $SLURM_ARRAY_TASK_ID exceeds number of samples ($TOTAL). Exiting."
    exit 0
fi

# --------------------
# Read sample line
# samples.tab format:
# sample_id<TAB>/path/to/assembly.fasta
# --------------------
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$TAB")
sample=$(echo "$line" | cut -f1)
assembly=$(echo "$line" | cut -f2)

echo "Running Snippy for sample: $sample"
echo "Assembly: $assembly"

# --------------------
# Run Snippy (ASSEMBLIES)
# --------------------
snippy \
  --outdir "$OUT/$sample" \
  --ref "$REF" \
  --ctgs "$assembly" \
  --cpus "$SLURM_CPUS_PER_TASK" \
  --force

echo "Finished sample: $sample"
