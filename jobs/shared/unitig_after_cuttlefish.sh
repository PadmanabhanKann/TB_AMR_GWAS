#!/bin/bash
#SBATCH --job-name=unitig_after_cuttlefish
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --partition=medium,long
#SBATCH --time=48:00:00
#SBATCH --error=/data/biol-micro-genomics/pkannan/logs/unitig_%j.txt
#SBATCH --output=/data/biol-micro-genomics/pkannan/logs/unitig_%j.txt
#SBATCH --mail-type=ALL
#SBATCH --mail-user=padmanabhan21@iisertvm.ac.in
# --- Load environment ---
module purge
module load Anaconda3/2023.09-0

cd /data/biol-micro-genomics/pkannan/pangwes2/pyseer

source activate unitig-caller
unitig-caller --call --refs trial.txt --unitigs unitigs.all.fa  --pyseer --rtab
