#!/bin/bash

#SBATCH --job-name=unitig-caller
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --partition=medium,long
#SBATCH --time=48:00:00
#SBATCH --error=/data/biol-micro-genomics/pkannan/logs/unitig_error_%j.txt
#SBATCH --output=/data/biol-micro-genomics/pkannan//logs/unitig_output_%j.txt
#SBATCH --mail-type=ALL
#SBATCH --mail-user=padmanabhan21@iisertvm.ac.in

# --- Load environment ---
module purge
module load Anaconda3/2023.09-0
source activate $DATA/env/unitig-caller

# --- Define paths ---
INPUT_FILE="/data/biol-micro-genomics/pkannan/TB/unitigs/paths.txt"
RESULT_DIR="/data/biol-micro-genomics/pkannan/TB/unitigs"

# --- Run unitig-caller ---
unitig-caller --call --refs "$INPUT_FILE" --out "$RESULT_DIR/unitig" --pyseer --kmer 31 --threads 16
