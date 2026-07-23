#!/bin/bash
#SBATCH --job-name=Fasttree_1k_multi_thread
#SBATCH --time=48:00:00            # adjust as needed (5 days here)
#SBATCH --mem=200G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=medium,long         
#SBATCH --output=/data/biol-micro-genomics/pkannan/logs/fastree_1k_%j.out
#SBATCH --error=/data/biol-micro-genomics/pkannan/logs/fastree_1k_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=padmanabhan21@iisertvm.ac.in


module purge

module load FastTree/2.1.11-GCCcore-11.3.0

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK


ALIGNMENT="/data/biol-micro-genomics/pkannan/TB/fasttree/core_gene_alignment.aln"

OUTPUT="/data/biol-micro-genomics/pkannan/TB/fasttree/core_gene_alignment.tree.nwk"

FastTree -nt -gtr -gamma "$ALIGNMENT" > "$OUTPUT"
