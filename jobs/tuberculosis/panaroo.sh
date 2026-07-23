#!/bin/bash
#SBATCH --job-name=panaroo
#SBATCH --time=5-00:00:00            # adjust as needed (5 days here)
#SBATCH --partition=long
#SBATCH --cpus-per-task=16           # threads   
#SBATCH --output=/data/biol-micro-genomics/pkannan/logs/panaroo_%j.out
#SBATCH --error=/data/biol-micro-genomics/pkannan/logs/panaroo_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=padmanabhan21@iisertvm.ac.in

set -euo pipefail

# Load environment
module purge
module load Mamba/4.14.0-0
source activate $DATA/env/panaroo_py

# Define paths
OUT_DIR=/data/biol-micro-genomics/pkannan/TB/panaroo 
INPUT_DIR=/data/biol-micro-genomics/pkannan/TB/panaroo/combined_gff/

mkdir -p "$OUT_DIR"

# Run Panaroo
panaroo -i ${INPUT_DIR}/*.gff \
        -o $OUT_DIR \
        -a core \
        --clean-mode strict \
        --aligner mafft \
        -t $SLURM_CPUS_PER_TASK
