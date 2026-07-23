#!/bin/bash

#SBATCH --job-name=superdca
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --partition=short,medium,long
#SBATCH --time=0-10:00:00
#SBATCH --error=/data/biol-micro-genomics/pkannan/logs/superdca_%j.err
#SBATCH --output=/data/biol-micro-genomics/pkannan/logs/superdca_%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=padmanabhan21@iisertvm.ac.in

# --- Load environment ---
module purge
module load Anaconda3/2023.09-0

cd /data/biol-micro-genomics/pkannan/TB/superdca
/data/biol-micro-genomics/pkannan/SuperDCA/binaries/SuperDCA_v0.1.1_master-84e7342/SuperDCA  -v /data/biol-micro-genomics/pkannan/TB/panaroo/core_gene_alignment.aln
