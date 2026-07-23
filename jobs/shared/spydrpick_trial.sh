#!/bin/bash
#SBATCH --job-name=spyderpick_trial
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --partition=short,medium,long
#SBATCH --time=03:00:00
#SBATCH --error=/data/biol-micro-genomics/pkannan/logs/spyderpick_trial.err
#SBATCH --output=/data/biol-micro-genomics/pkannan/logs/spyderpick_trial.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=padmanabhan21@iisertvm.ac.in
# --- Load environment ---
module purge
module load Anaconda3/2023.09-0
source activate pangwes

cd /data/biol-micro-genomics/pkannan/pangwes/trial

SpydrPick --alignmentfile  trial.fa  --verbose        --threads 12
  
