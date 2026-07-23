#!/bin/bash
#SBATCH --job-name=gemma
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --partition=short,medium,long
#SBATCH --time=0-10:00:00
#SBATCH --output=/data/biol-micro-genomics/pkannan/logs/gemma_tb_%j.out
#SBATCH --error=/data/biol-micro-genomics/pkannan/logs/gemma_tb_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=padmanabhan21@iisertvm.ac.in

set -euo pipefail

BASE="/data/biol-micro-genomics/pkannan/TB/gemma"
BFILE="${BASE}/plink/mtb_nomaf"
KIN="${BASE}/output/mtb_kinship_nomaf.cXX.txt"
PHENODIR="${BASE}/phenotypes"
RESULTS="${BASE}/results"

mkdir -p "${RESULTS}" /data/biol-micro-genomics/pkannan/logs

[[ -f "${BFILE}.bed" ]] || { echo "Missing ${BFILE}.bed"; exit 1; }
[[ -f "${KIN}" ]] || { echo "Missing kinship ${KIN}"; exit 1; }

for PHENO in "${PHENODIR}"/plink_*_log.pheno; do
  F="$(basename "${PHENO}")"
  AB="${F#plink_}"
  AB="${AB%_log.pheno}"

  AB_DIR="${RESULTS}/${AB}"
  mkdir -p "${AB_DIR}/work" "${AB_DIR}/output"

  echo "===== Running ${AB} ====="

  plink --bfile "${BFILE}" --allow-extra-chr \
        --pheno "${PHENO}" --make-bed \
        --out "${AB_DIR}/work/mtb_${AB}"

  ( cd "${AB_DIR}" && \
    gemma -bfile "${AB_DIR}/work/mtb_${AB}" \
          -k "${KIN}" -lmm 4 -o "${AB}_lmm" )
done
