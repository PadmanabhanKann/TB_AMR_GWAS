#!/bin/bash

#SBATCH --job-name=cond_pyseer
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --partition=short,medium,long
#SBATCH --time=0-10:00:00
#SBATCH --error=/data/biol-micro-genomics/pkannan/logs/cond_pyseer_%j.err
#SBATCH --output=/data/biol-micro-genomics/pkannan/logs/cond_pyseer_%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=padmanabhan21@iisertvm.ac.in

# --- Load environment ---
module purge
module load Anaconda3/2023.09-0
source activate $DATA/env/pyseer

cd /data/biol-micro-genomics/pkannan/Azi_ngon2/clonal_pyseer/cond_pyseer

pyseer --lmm --use-covariates 2 3   --phenotypes phenotype.txt   --kmers unitig.pyseer.gz   --similarity  phylogeny_k.tsv   --covariates covariate.tsv   --output-patterns patterns.txt   > conditional_gwas_results.txt
