#!/bin/bash

#SBATCH --job-name=ld_decay
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --partition=short,medium,long
#SBATCH --time=0-10:00:00
#SBATCH --error=/data/biol-micro-genomics/pkannan/logs/ld_decay_%j.err
#SBATCH --output=/data/biol-micro-genomics/pkannan/logs/ld_decay_%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=padmanabhan21@iisertvm.ac.in

# --- Load environment ---
module purge
module load Anaconda3/2023.09-0
source activate $DATA/env/python2

cd /data/biol-micro-genomics/pkannan/TB/plink
/data/biol-micro-genomics/pkannan/env/python2/bin/python ld_decay_calc.py -i mtb_ld.ld.gz -o decay
