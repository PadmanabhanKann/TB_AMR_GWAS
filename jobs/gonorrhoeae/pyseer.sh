#!/bin/bash

#SBATCH --job-name=pyseer
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --partition=short,medium,long
#SBATCH --time=0-10:00:00
#SBATCH --error=/data/biol-micro-genomics/pkannan/logs/pyseer_%j.err
#SBATCH --output=/data/biol-micro-genomics/pkannan/logs/pyseer_%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=padmanabhan21@iisertvm.ac.in

# --- Load environment ---
module purge
module load Anaconda3/2023.09-0
source activate $DATA/env/pyseer

cd /data/biol-micro-genomics/pkannan/Azi_ngon2/clonal_pyseer/pyseer
pyseer --lmm --phenotypes phenotype.txt --similarity phylogeny_k.tsv --kmers unitig.pyseer.gz --output-patterns azithro_kmer_patterns.txt --cpu 16 > azithro_kmers_lmm_results.txt
