#!/bin/bash
#SBATCH --job-name=genepi
#SBATCH --time=10:00:00            # adjust as needed (5 days here)
#SBATCH --nodes=1
#SBATCH --mem=100G
#SBATCH --partition=short,medium,long         
#SBATCH --output=/data/biol-micro-genomics/pkannan/logs/genepi_%j.out
#SBATCH --error=/data/biol-micro-genomics/pkannan/logs/genepi_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=padmanabhan21@iisertvm.ac.in


module purge
source activate genepi

cd /data/biol-micro-genomics/pkannan/TB/epiML

GenEpi -g mtb.gen \
       -p test.csv \
       -o ./genepi_RIF \
       -m r 
