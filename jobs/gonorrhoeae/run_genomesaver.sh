#!/bin/bash

#SBATCH --job-name=genomesaver
#SBATCH --nodes=1
#SBATCH --partition=short,medium,long
#SBATCH --time=0-10:00:00
#SBATCH --error=/data/biol-micro-genomics/pkannan/logs/genomesaver.error
#SBATCH --output=/data/biol-micro-genomics/pkannan/logs/genomesaver.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=padmanabhan21@iisertvm.ac.in

module purge
cd /data/biol-micro-genomics/pkannan/scripts

bash ./genomesaver.sh /data/biol-micro-genomics/pkannan/Azithro_N.gon/genome_ids.txt  /data/biol-micro-genomics/pkannan/Azithro_N.gon/genomes 
