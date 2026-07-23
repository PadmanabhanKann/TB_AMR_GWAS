#!/bin/bash
#SBATCH --job-name=quast_ng_batches
#SBATCH --array=1-15%8             # 7369/500 ≈ 15 batches; run up to 4 at once
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
#SBATCH --partition=medium,long
#SBATCH --output=/data/biol-micro-genomics/pkannan/logs/quast_batch_%A_%a.out
#SBATCH --error=/data/biol-micro-genomics/pkannan/logs/quast_batch_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=padmanabhan21@iisertvm.ac.in

set -euo pipefail

module purge
module load QUAST/5.0.2-foss-2020a-Python-3.8.2

# ABSOLUTE paths
GENOMES_LIST=/data/biol-micro-genomics/pkannan/Azithro_N.gon/genome_paths.txt
BATCH_SIZE=500
BASE_OUT=/data/biol-micro-genomics/pkannan/Azithro_N.gon/quast_results

# ensure dirs exist
mkdir -p /data/biol-micro-genomics/pkannan/logs "$BASE_OUT"

# preflight
echo "Host=$(hostname)"
echo "Task=$SLURM_ARRAY_TASK_ID  CPUs=$SLURM_CPUS_PER_TASK"
echo "GENOMES_LIST=$GENOMES_LIST  BATCH_SIZE=$BATCH_SIZE"

# validate input list
if [[ ! -r "$GENOMES_LIST" ]]; then
  echo "ERROR: Cannot read GENOMES_LIST at $GENOMES_LIST" >&2
  exit 2
fi

N=$(wc -l < "$GENOMES_LIST")
echo "Total genomes (N)=$N"

START=$(( (SLURM_ARRAY_TASK_ID - 1) * BATCH_SIZE + 1 ))
if (( START > N )); then
  echo "No files for batch ${SLURM_ARRAY_TASK_ID} (start=$START > N=$N)"; exit 0
fi
END=$(( SLURM_ARRAY_TASK_ID * BATCH_SIZE ))
if (( END > N )); then END=$N; fi
echo "Slice: START=$START END=$END"

mapfile -t FILES < <(sed -n "${START},${END}p" "$GENOMES_LIST")
echo "FILES in this batch=${#FILES[@]}"

OUTDIR="${BASE_OUT}/batch_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$OUTDIR"

# run QUAST
quast "${FILES[@]}" -o "$OUTDIR" -t "$SLURM_CPUS_PER_TASK" --silent --fast --min-contig 500
