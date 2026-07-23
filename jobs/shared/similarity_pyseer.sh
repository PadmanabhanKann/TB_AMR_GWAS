#!/bin/bash

#SBATCH --job-name=sim_pyseer
#SBATCH --nodes=1
#SBATCH  --mem=100G
#SBATCH --ntasks-per-node=12
#SBATCH --partition=short,medium,long
#SBATCH --time=48:00:00
#SBATCH --error=/data/biol-micro-genomics/pkannan/logs/sim_pyseer_%j.err
#SBATCH --output=/data/biol-micro-genomics/pkannan/logs/sim_pyseer_%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=padmanabhan21@iisertvm.ac.in

# --- Load environment ---
module purge
module load Anaconda3/2023.09-0
source activate $DATA/env/pyseer

cd /data/biol-micro-genomics/pkannan/sra/pyseer2
similarity_pyseer --kmers sra_final.pyseer.gz  genome_ids.txt > kinship.tsv
