#!/bin/bash
#SBATCH --job-name=quast_test
#SBATCH --partition=devel
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=quast_test_%j.out
#SBATCH --error=quast_test_%j.err

set -euo pipefail

module purge
module load QUAST/5.0.2-foss-2020a-Python-3.8.2

GENOMES_LIST=/data/biol-micro-genomics/pkannan/Azithro_N.gon_10k/genomes.txt
TEST_OUT=/data/biol-micro-genomics/pkannan/quast_results/test_run

mkdir -p "$TEST_OUT"

# Take just the first 5 genomes to test
mapfile -t FILES < <(head -n 5 "$GENOMES_LIST")

echo "Running QUAST on ${#FILES[@]} files:"
printf '%s\n' "${FILES[@]}"

which quast
quast --version

quast "${FILES[@]}" -o "$TEST_OUT" -t $SLURM_CPUS_PER_TASK --silent --fast
