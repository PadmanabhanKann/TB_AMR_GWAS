#!/bin/bash
#SBATCH --job-name=pyseer_tb
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --partition=short,medium,long
#SBATCH --time=0-10:00:00
#SBATCH --error=/data/biol-micro-genomics/pkannan/logs/pyseer_tb_%A_%a.err
#SBATCH --output=/data/biol-micro-genomics/pkannan/logs/pyseer_tb_%A_%a.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=padmanabhan21@iisertvm.ac.in
#SBATCH --array=0-13

set -euo pipefail

BASE="/data/biol-micro-genomics/pkannan/TB/pyseer2"
KMER="${BASE}/unitig.pyseer.gz"
SIM="${BASE}/phylogeny_K.tsv"
CPU="${SLURM_CPUS_PER_TASK}"

module purge
module load Anaconda3/2023.09-0
source activate /data/biol-micro-genomics/pkannan/env/pyseer

cd "$BASE"

# one CSV per antibiotic folder
mapfile -t PHENOS < <(find "$BASE" -mindepth 2 -maxdepth 2 -type f -name "*.tsv" | sort)

PHENO="${PHENOS[$SLURM_ARRAY_TASK_ID]}"
AB_DIR="$(dirname "$PHENO")"
AB="$(basename "$AB_DIR")"

cd "$AB_DIR"

OUT_PREFIX="${AB}_kmers"
pyseer \
  --lmm \
  --phenotypes "$PHENO" \
  --similarity "$SIM" \
  --kmers "$KMER" \
  --output-patterns "${OUT_PREFIX}_patterns.txt" \
  --cpu "$CPU" \
  > "${OUT_PREFIX}_results.txt"
