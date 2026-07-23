#!/bin/bash
#SBATCH --job-name=iqtree
#SBATCH --time=48:00:00            # adjust as needed (5 days here)
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --partition=medium,long         
#SBATCH --output=/data/biol-micro-genomics/pkannan/logs/iqtree_%j.out
#SBATCH --error=/data/biol-micro-genomics/pkannan/logs/iqtree_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=padmanabhan21@iisertvm.ac.in


module purge
module load Anaconda3/2023.09-0
source activate iqtree



cd /data/biol-micro-genomics/pkannan/TB/tree
INPUT="/data/biol-micro-genomics/pkannan/TB/tree/core_gene_alignment.aln"
OUTPUT="/data/biol-micro-genomics/pkannan/TB/tree"

iqtree -s "$INPUT" -m GTR+G -T 24 -redo  -bb 1000  -pre "$(basename "$INPUT" .phy)_iqtree"

