#!/bin/bash

#SBATCH --job-name=pyseer_test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --partition=devel
#SBATCH --time=10:00
#SBATCH --error=/data/biol-micro-genomics/pkannan/logs/test_pyseer_%j.err
#SBATCH --output=/data/biol-micro-genomics/pkannan/logs/test_pyseer_%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=padmanabhan21@iisertvm.ac.in

# --- Load environment ---
module purge
module load Anaconda3/2023.09-0
source activate $DATA/env/pyseer

pyseer
