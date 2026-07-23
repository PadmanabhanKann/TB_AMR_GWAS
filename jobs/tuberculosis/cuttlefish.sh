#!/bin/bash
#SBATCH --job-name=cuttlefish
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --partition=medium,long
#SBATCH --time=48:00:00
#SBATCH --error=/data/biol-micro-genomics/pkannan/logs/cuttlefish_%j.txt
#SBATCH --output=/data/biol-micro-genomics/pkannan/logs/cuttlefish_%j.txt
#SBATCH --mail-type=ALL
#SBATCH --mail-user=padmanabhan21@iisertvm.ac.in
# --- Load environment ---
module purge
module load Anaconda3/2023.09-0
source activate pangwes

cd /data/biol-micro-genomics/pkannan/TB/pangwes

cuttlefish build   --list genome_paths.txt    --format 1  --kmer-len 31   --output cdbg_k31_dbg
