#!/bin/bash

#SBATCH --job-name=ClonalframeML
#SBATCH --time=48:00:00            # adjust as needed (5 days here)
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=medium,long         
#SBATCH --output=/data/biol-micro-genomics/pkannan/logs/ClonalFrameML_%j.out
#SBATCH --error=/data/biol-micro-genomics/pkannan/logs/ClonalFrameML_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=padmanabhan21@iisertvm.ac.in


module purge

source activate clonalframeml


cd /data/biol-micro-genomics/pkannan/Azi_ngon2/clonal_pyseer

ClonalFrameML core_gene_alignment.tree.nwk core_gene_alignment.aln corrected.nwk
