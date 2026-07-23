#!/bin/bash
#SBATCH --job-name=pyseer_cond_rifamycins
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --partition=short,medium,long
#SBATCH --time=0-10:00:00
#SBATCH --error=/data/biol-micro-genomics/pkannan/logs/pyseer_cond_rifamycins_%j.err
#SBATCH --output=/data/biol-micro-genomics/pkannan/logs/pyseer_cond_rifamycins_%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=padmanabhan21@iisertvm.ac.in
#SBATCH --array=0-1

set -euo pipefail

BASE="/data/biol-micro-genomics/pkannan/TB/cond_gwas2"
KMER="${BASE}/unitig.pyseer.gz"
SIM="${BASE}/phylogeny_K.tsv"
CPU="${SLURM_CPUS_PER_TASK}"

ANTIBIOTICS=(rifampicin rifabutin)
AB="${ANTIBIOTICS[$SLURM_ARRAY_TASK_ID]}"
AB_DIR="${BASE}/${AB}"

PHENO="${AB_DIR}/${AB}.tsv"
COV="${AB_DIR}/${AB}_subset.covariates.tsv"

module purge || true
module load Anaconda3/2023.09-0 || true
source activate /data/biol-micro-genomics/pkannan/env/pyseer

# Sanity checks
[[ -s "$PHENO" ]] || { echo "ERROR: Missing phenotype: $PHENO"; exit 1; }
[[ -s "$COV"   ]] || { echo "ERROR: Missing covariates: $COV";  exit 1; }
[[ -s "$KMER"  ]] || { echo "ERROR: Missing kmers: $KMER";       exit 1; }
[[ -s "$SIM"   ]] || { echo "ERROR: Missing similarity: $SIM";   exit 1; }

echo "Antibiotic: $AB"
echo "Covariate file: $COV"
echo "Columns: $(head -n1 "$COV" | awk -F'\t' '{print NF}')"
echo "First covariate name: $(head -n1 "$COV" | awk -F'\t' '{print $2}')"
echo "Sample count: $(tail -n+2 "$COV" | wc -l)"

cd "$AB_DIR"

N=$(head -n1 "$COV" | awk -F'\t' '{print NF}')
[[ "$N" -ge 2 ]] || { echo "ERROR: Covariates file has <2 columns"; exit 1; }
USE=$(seq 2 "$N" | paste -sd' ' -)

echo "Using covariate columns: $USE"
echo "Running pyseer conditional GWAS for ${AB}..."

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

echo "Done. Output: ${AB_DIR}/${OUT_PREFIX}_results.txt"

KMER_SEQ="ATTAGGGTGGGGTCGTTCCGGGGTCGGTGGG"
echo ""
echo "--- esxQ k-mer check in new ${AB} results ---"
grep "$KMER_SEQ" "${OUT_PREFIX}_results.txt" \
  || echo "esxQ k-mer not found (did not pass filter)"

echo ""
echo "--- Comparison with standard GWAS ---"
STD_FILE="/data/biol-micro-genomics/pkannan/TB/pyseer2/${AB}/annotated_kmers.txt"
if [[ -f "$STD_FILE" ]]; then
    STD=$(grep -m1 "$KMER_SEQ" "$STD_FILE" | awk '{print $3,$4,$5,$6}')
    COND=$(grep -m1 "$KMER_SEQ" "${OUT_PREFIX}_results.txt" | awk '{print $3,$4,$5,$6}')
    echo "Standard    (filter_p lrt_p beta se): $STD"
    echo "Conditional (filter_p lrt_p beta se): $COND"
    if [[ "$STD" == "$COND" ]]; then
        echo "WARNING: values identical - covariate may not have taken effect"
    else
        echo "OK: values differ - conditional GWAS ran correctly"
    fi
else
    echo "Standard GWAS file not found at $STD_FILE - skipping comparison"
fi
