#!/bin/bash
#SBATCH --job-name=blast
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --partition=short,medium,long
#SBATCH --time=08:00:00
#SBATCH --error=/data/biol-micro-genomics/pkannan/logs/blast_%j.txt
#SBATCH --output=/data/biol-micro-genomics/pkannan/logs/blast_%j.txt
#SBATCH --mail-type=ALL
#SBATCH --mail-user=padmanabhan21@iisertvm.ac.in

module purge
module load BLAST+/2.14.0-gompi-2022a

cd /data/biol-micro-genomics/pkannan/pangwes2/map_uni
mkdir -p results

blastn -task blastn-short \
  -query unitigs.fa \
  -db rrna_db \
-perc_identity 97 \
  -evalue 1e-10 \
  -qcov_hsp_perc 80 \
  -num_threads 8 \
  -outfmt 6 \
  >results/blast_rrna.tsv
