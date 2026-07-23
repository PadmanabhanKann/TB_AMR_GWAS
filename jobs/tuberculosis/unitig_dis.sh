#!/bin/bash
#SBATCH --job-name=unitig_dis
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --partition=short,medium,long
#SBATCH --time=10:00:00
#SBATCH --error=/data/biol-micro-genomics/pkannan/logs/unitig_dis_%j.txt
#SBATCH --output=/data/biol-micro-genomics/pkannan/logs/unitig_dis_%j.txt
#SBATCH --mail-type=ALL
#SBATCH --mail-user=padmanabhan21@iisertvm.ac.in
# --- Load environment ---
module purge
module load Anaconda3/2023.09-0
source activate pangwes

cd /data/biol-micro-genomics/pkannan/TB/pangwes

 unitig_distance --unitigs-file    cdbg_31.unitigs                \
                --edges-file      cdbg_31.edges                  \
                --k-mer-length    31                     \
                --sgg-paths-file  cdbg_31.paths                  \
                --queries-file cdbg_31.filtered_ge001maf_le015gf_gt01-lt05states.L121282n3080.spydrpick_couplings.1-based.12766521edges \
                --threads         8                               \
                --queries-one-based                               \
                --run-sggs-only                                   \
                --output-stem cdgb                              \
                --verbose
          
