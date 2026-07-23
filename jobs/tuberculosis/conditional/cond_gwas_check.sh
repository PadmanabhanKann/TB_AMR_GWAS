#!/bin/bash
#SBATCH --job-name=pyseer_cond_clofmerge
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --partition=short,medium,long
#SBATCH --time=0-10:00:00
#SBATCH --error=/data/biol-micro-genomics/pkannan/logs/pyseer_cond_clofmerge_%A_%a.err
#SBATCH --output=/data/biol-micro-genomics/pkannan/logs/pyseer_cond_clofmerge_%A_%a.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=padmanabhan21@iisertvm.ac.in
#SBATCH --array=0-2

set -euo pipefail




#tocheck effect of metK

BASE="/data/biol-micro-genomics/pkannan/TB/cond_gwas2"
KMER="${BASE}/unitig.pyseer.gz"
SIM="${BASE}/phylogeny_K.tsv"
CPU="${SLURM_CPUS_PER_TASK}"

module purge || true
module load Anaconda3/2023.09-0 || true
source activate /data/biol-micro-genomics/pkannan/env/pyseer

ABLIST=(ethambutol linezolid moxifloxacin)
AB="${ABLIST[$SLURM_ARRAY_TASK_ID]}"
AB_DIR="${BASE}/${AB}"

PHENO="${AB_DIR}/${AB}.tsv"
COV="${AB_DIR}/check/${AB}_subset.covariates_plus_clofazimine.tsv"

[[ -s "$PHENO" ]] || { echo "Missing phenotype: $PHENO"; exit 1; }
[[ -s "$COV"   ]] || { echo "Missing covariates: $COV"; exit 1; }

N=$(head -n1 "$COV" | awk -F'\t' '{print NF}')
[[ "$N" -ge 2 ]] || { echo "Covariates file has <2 columns: $COV"; exit 1; }
USE=$(seq 2 "$N" | paste -sd' ' -)

OUT_PREFIX="${AB_DIR}/check/${AB}_kmers_cond_plus_clof"

cd "$AB_DIR"

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
