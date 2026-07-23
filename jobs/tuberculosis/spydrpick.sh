#!/bin/bash
#SBATCH --job-name=spydrpick
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --partition=short,medium,long
#SBATCH --time=10:00:00
#SBATCH --error=/data/biol-micro-genomics/pkannan/logs/spydrpick_%j.txt
#SBATCH --output=/data/biol-micro-genomics/pkannan/logs/spydrpick_%j.txt
#SBATCH --mail-type=ALL
#SBATCH --mail-user=padmanabhan21@iisertvm.ac.in
# --- Load environment ---
module purge
module load Anaconda3/2023.09-0
source activate pangwes

cd /data/biol-micro-genomics/pkannan/TB/pangwes2

SpydrPick --alignmentfile  cdbg_31.fasta \
            --ld-threshold    20000      \
            --threads        8              \
            --verbose  
