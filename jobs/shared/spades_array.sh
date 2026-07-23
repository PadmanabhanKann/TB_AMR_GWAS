#!/bin/bash
#SBATCH --job-name=spades_refguided
#SBATCH --array=1-730%50           # adjust 730 to wc -l of your list; %50 = max concurrent tasks
#SBATCH --cpus-per-task=12
#SBATCH --mem=16G                  # no spaces around '='
#SBATCH --time=04:00:00
#SBATCH --partition=short,medium,long         # pick ONE partition; remove/adjust as needed
#SBATCH --output=/data/biol-micro-genomics/pkannan/logs/spades_%A_%a.out
#SBATCH --error=/data/biol-micro-genomics/pkannan/logs/spades_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=padmanabhan21@iisertvm.ac.in

set -euo pipefail

module purge
module load SPAdes/3.15.2-GCC-10.2.0

REF=/data/biol-micro-genomics/pkannan/sra/spades/reference.fa
LIST=/data/biol-micro-genomics/pkannan/sra/spades/fixed_paired.txt

THREADS="${SLURM_CPUS_PER_TASK}"
MEM_GB=64

mkdir -p /data/biol-micro-genomics/pkannan/logs
mkdir -p /data/biol-micro-genomics/pkannan/sra/spades/assemblies
mkdir -p /data/biol-micro-genomics/pkannan/sra/spades/scaffolds_renamed

# Count number of lines in LIST
TOTAL=$(wc -l < "$LIST")

# Guard: exit if this task ID is larger than the number of lines
if (( SLURM_ARRAY_TASK_ID > TOTAL )); then
  echo "[$(date)] Task ${SLURM_ARRAY_TASK_ID} > ${TOTAL} lines in $LIST; exiting."
  exit 0
fi

# Pull the R1,R2 pair for this array index
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST")
R1=$(echo "$LINE" | cut -d, -f1)
R2=$(echo "$LINE" | cut -d, -f2)

# Derive sample name from R1; handles .fastq or .fastq.gz
fname=$(basename "$R1")
sample="${fname%_1.fastq}"
sample="${sample%_1.fastq.gz}"

outdir="/data/biol-micro-genomics/pkannan/sra/spades/assemblies/spades_refguided_${sample}"

echo "[$(date)] Assembling ${sample}"
echo "R1: $R1"
echo "R2: $R2"
echo "Out: $outdir"

# Skip if already finished
if [[ -s "${outdir}/scaffolds.fasta" ]]; then
  echo "scaffolds.fasta exists for ${sample}; renaming/copying and skipping SPAdes."
else
  spades.py \
    -1 "$R1" \
    -2 "$R2" \
    --trusted-contigs "$REF" \
    --cov-cutoff 10 \
    -o "$outdir" \
    --threads "$THREADS" \
    --memory "$MEM_GB" \
    --phred-offset 33 
fi

# Rename/copy scaffolds to include the sample id
if [[ -s "${outdir}/scaffolds.fasta" ]]; then
  dest="/data/biol-micro-genomics/pkannan/sra/spades/scaffolds_renamed/scaffolds_${sample}.fasta"
  cp -f "${outdir}/scaffolds.fasta" "$dest"
  echo "[$(date)] Finished ${sample}: $(realpath "$dest")"
else
  echo "[$(date)] ERROR: ${outdir}/scaffolds.fasta not found or empty." >&2
  exit 2
fi
