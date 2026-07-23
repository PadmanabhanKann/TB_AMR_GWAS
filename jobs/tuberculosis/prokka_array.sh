#!/bin/bash
#SBATCH --job-name=prokka_ng
#SBATCH --array=1-350%53
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G
#SBATCH --time=10:00:00
#SBATCH --partition=short,medium,long
#SBATCH --output=/data/biol-micro-genomics/pkannan/logs/prokka_%A_%a.out
#SBATCH --error=/data/biol-micro-genomics/pkannan/logs/prokka_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=padmanabhan21@iisertvm.ac.in
set -euo pipefail

# --- HARD-CODED PATHS ---
INPUT_DIR=/data/biol-micro-genomics/pkannan/TB/assemblies
LIST_FILE=/data/biol-micro-genomics/pkannan/TB/genomepaths.txt
OUT_ROOT=/data/biol-micro-genomics/pkannan/TB/prokka
PANAROO_GFF_DIR=/data/biol-micro-genomics/pkannan/TB/panaroo/combined_gff
LOG_DIR=/data/biol-micro-genomics/pkannan/logs


module purge
module load Mamba/4.14.0-0
source activate $DATA/env/prokka

mkdir -p "$OUT_ROOT" "$PANAROO_GFF_DIR" "$LOG_DIR"

# Build list once (absolute paths of .fna)
if [[ ! -s "$LIST_FILE" ]]; then
  find "$INPUT_DIR" -maxdepth 1 -type f -regextype posix-extended \
    -iregex '.*\.(fas|fa|fna|fasta)(\.gz)?$' | sort > "$LIST_FILE"
fi
# ---- Batch K genomes per task ----
BATCH_SIZE=10                                           # << change if you want 3, 5, 10, ...
N=$(wc -l < "$LIST_FILE")
START=$(( (SLURM_ARRAY_TASK_ID - 1) * BATCH_SIZE + 1 ))
END=$(( START + BATCH_SIZE - 1 ))
(( END > N )) && END=$N

# Load this slice
mapfile -t FILES < <(sed -n "${START},${END}p" "$LIST_FILE")

# Nothing to do? exit cleanly
((${#FILES[@]}==0)) && { echo "No files for task ${SLURM_ARRAY_TASK_ID} (start=$START > N=$N)"; exit 0; }

echo "Node=$(hostname)"
echo "Task=${SLURM_ARRAY_TASK_ID}  Slice: $START..$END  Files=${#FILES[@]}"

# ---- Process each file in this batch sequentially ----
for FILE in "${FILES[@]}"; do
  [[ -z "${FILE:-}" ]] && continue

  base=$(basename "$FILE")
  sample=${base%.*}

  SAMPLE_OUT="${OUT_ROOT}/${sample}"
  mkdir -p "$SAMPLE_OUT"

  echo ">>> Prokka: $sample"
prokka --outdir "$SAMPLE_OUT" \
       --prefix "$sample" \
       --cpus "$SLURM_CPUS_PER_TASK" \
       --genus Mycobacterium --species tuberculosis \
       --force \
       "$FILE"

  gff="${SAMPLE_OUT}/${sample}.gff"
  if [[ -f "$gff" ]]; then
    cp -f "$gff" "$PANAROO_GFF_DIR"/
  else
    echo "WARNING: No GFF produced for $sample" >&2
  fi
  echo "<<< Done: $sample"
done
