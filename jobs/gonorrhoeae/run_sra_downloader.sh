#!/bin/bash
#SBATCH --job-name=sra_genomesaver
#SBATCH --nodes=1
#SBATCH --partition=short,medium,long
#SBATCH --time=0-10:00:00
#SBATCH --error=/data/biol-micro-genomics/pkannan/logs/sra_genomesaver_%j.error
#SBATCH --output=/data/biol-micro-genomics/pkannan/logs/sra_genomesaver_%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=padmanabhan21@iisertvm.ac.in

module purge
cd /data/biol-micro-genomics/pkannan/scripts
source activate sra

# down1.sh was folded into the generic download_sra.sh (see preprocessing/sra_download/) — pass its azithromycin-specific input list explicitly
bash ./download_sra.sh -i /data/biol-micro-genomics/pkannan/sra/azi_rawreads.txt
