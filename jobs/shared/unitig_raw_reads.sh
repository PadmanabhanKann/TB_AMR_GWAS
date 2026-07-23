#!/bin/bash
#SBATCH --job-name=unitig_raw_reads
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=16
#SBATCH --partition=long
#SBATCH  --mem=100G 
#SBATCH --time=7-00:00:00 
#SBATCH --error=/data/biol-micro-genomics/pkannan/logs/unitig_raw_error_%j.txt 
#SBATCH --output=/data/biol-micro-genomics/pkannan/logs/unitig_raw_output_%j.txt 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=padmanabhan21@iisertvm.ac.in
# --- Load environment ---
module purge 
module load Anaconda3/2023.09-0 
source activate $DATA/env/unitig-caller
# --- Define paths ---
INPUT_FILE="/data/biol-micro-genomics/pkannan/sra/dashing2/paths.txt" 
RESULT_DIR="/data/biol-micro-genomics/pkannan/sra/unitig"
# --- Run unitig-caller ---
unitig-caller --call --reads "$INPUT_FILE" --out "$RESULT_DIR/unitig" --kmer 31 --threads 16  --pyseer --rtab 
