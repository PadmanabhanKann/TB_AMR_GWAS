#!/bin/bash

#SBATCH --job-name=raw_read_unitig-caller
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH 
#SBATCH --partition=short,medium,long
#SBATCH --time=0-01:00:00
#SBATCH --error=/data/biol-micro-genomics/pkannan/unitig_error.txt
#SBATCH --output=/data/biol-micro-genomics/pkannan/unitig_output.txt
#SBATCH --mail-type=ALL
#SBATCH --mail-user=padmanabhan21@iisertvm.ac.in
