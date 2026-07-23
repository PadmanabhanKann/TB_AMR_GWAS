#!/bin/bash
#SBATCH --job-name=dashing2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --time=03:00:00
#SBATCH --partition=short,medium,long
#SBATCH --error=/data/biol-micro-genomics/pkannan/logs/dashing2_error_%j.txt
#SBATCH --output=/data/biol-micro-genomics/pkannan/dashing2_output_%j.txt
#SBATCH --mail-type=ALL
#SBATCH --mail-user=padmanabhan21@iisertvm.ac.in


module purge
module load Anaconda3/2023.09-0
source activate dashing 

cd  /data/biol-micro-genomics/pkannan/sra/dashing2


realpath ../interleaved_reads/*.fastq >paths.txt

dashing2 sketch -F paths.txt -k 31 -p 12  -o all_samples.h5

dashing2 cmp --presketched all_samples.h5 --distance -p 32 -o distances_pyseer.tsv

