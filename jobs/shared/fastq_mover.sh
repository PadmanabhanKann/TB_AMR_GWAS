#!/bin/bash
#SBATCH --job-name=fastq_mover
#SBATCH --nodes=1
#SBATCH --partition=short,medium,long
#SBATCH --time=0-12:00:00
#SBATCH --error=/data/biol-micro-genomics/pkannan/logs/fastq_mover_%j.error
#SBATCH --output=/data/biol-micro-genomics/pkannan/logs/fastq_mover_%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=padmanabhan21@iisertvm.ac.in

module  purge
cd /data/biol-micro-genomics/pkannan/sra/downloads



find . -type f -name "*fastq" -exec cp -n -t ../fastq/ {} +
