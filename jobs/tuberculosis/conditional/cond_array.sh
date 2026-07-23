#!/bin/bash 
#SBATCH --job-name=pyseer_cond
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --partition=short,medium,long
#SBATCH --time=0-10:00:00
#SBATCH --error=/data/biol-micro-genomics/pkannan/logs/pyseer_cond_%A_%a.err
#SBATCH --output=/data/biol-micro-genomics/pkannan/logs/pyseer_cond_%A_%a.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=padmanabhan21@iisertvm.ac.in
#SBATCH --array=0-12

set -euo pipefail

BASE="/data/biol-micro-genomics/pkannan/TB/cond_gwas2"
KMER="${BASE}/unitig.pyseer.gz"
SIM="${BASE}/phylogeny_K.tsv"
CPU="${SLURM_CPUS_PER_TASK}"

module purge || true
module load Anaconda3/2023.09-0 || true
source activate /data/biol-micro-genomics/pkannan/env/pyseer

cd "$BASE"

# Antibiotic folders (phenotype-free discovery)
mapfile -t ABDIRS < <(find "$BASE" -mindepth 1 -maxdepth 1 -type d | sort)
AB_DIR="${ABDIRS[$SLURM_ARRAY_TASK_ID]}"
AB="$(basename "$AB_DIR")"

PHENO="${AB_DIR}/${AB}.tsv"
COV="${AB_DIR}/${AB}_subset.covariates.tsv"

# Skip non-antibiotic dirs if any exist at BASE level
[[ "$AB" =~ ^(plot|data|ref|logs|bubbleplots|esxQ_rank)$ ]] && exit 0

[[ -s "$PHENO" ]] || { echo "Missing phenotype: $PHENO"; exit 1; }
[[ -s "$COV"   ]] || { echo "Missing covariates: $COV"; exit 1; }

cd "$AB_DIR"

# build "2 3 ... N" for --use-covariates (skip sample column)
N=$(head -n1 "$COV" | awk -F'\t' '{print NF}')
[[ "$N" -ge 2 ]] || { echo "Covariates file has <2 columns: $COV"; exit 1; }
USE=$(seq 2 "$N" | paste -sd' ' -)

OUT_PREFIX="${AB}_kmers"

pyseer \
  --lmm \
  --phenotypes "$PHENO" \
  --similarity "$SIM" \
  --kmers "$KMER" \
  --covariates "$COV" \
  --use-covariates $USE \
  --output-patterns "${OUT_PREFIX}_patterns.txt" \
  --cpu "$CPU" \
  > "${OUT_PREFIX}_results.txt"
