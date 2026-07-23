#!/bin/bash
#SBATCH --job-name=interweave_fasq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --partition=short,medium,long
#SBATCH --time=0-10:00:00
#SBATCH --error=/data/biol-micro-genomics/pkannan/logs/interweave_%j.txt
#SBATCH --output=/data/biol-micro-genomics/pkannan/logs/interweave_%j.txt
#SBATCH --mail-type=ALL
#SBATCH --mail-user=padmanabhan21@iisertvm.ac.in




 module purge 
module load Anaconda3/2023.09-0 
source activate dashing 

cd /data/biol-micro-genomics/pkannan/sra/fastq
for fq1 in /data/biol-micro-genomics/pkannan/sra/fastq/*_1.fastq; do
  base="$(basename "${fq1%_1.fastq}")" # get sample name without path 
  fq2="/data/biol-micro-genomics/pkannan/sra/fastq/${base}_2.fastq" 
  out="../interleaved_reads/${base}.interleaved.fastq" 
  echo "Repairing/interleaving 
  $fq1 + $fq2 → $out"
  repair.sh in="$fq1" in2="$fq2" out="$out" overwrite=true
done
