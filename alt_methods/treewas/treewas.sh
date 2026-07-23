#!/bin/bash
#SBATCH --job-name=treewas
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --partition=short,medium,long
#SBATCH --time=0-10:00:00
#SBATCH --mem=250G
#SBATCH --array=1-13
#SBATCH --error=/data/biol-micro-genomics/pkannan/logs/treewas_%A_%a.err
#SBATCH --output=/data/biol-micro-genomics/pkannan/logs/treewas_%A_%a.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=padmanabhan21@iisertvm.ac.in

set -eo pipefail
# --- Load environment ---
module purge
module load Anaconda3/2023.09-0
source activate $DATA/env/treewas

ABS=(amikacin bedaquiline clofazimine delamanid ethambutol ethionamide isoniazid kanamycin levofloxacin linezolid moxifloxacin rifabutin rifampicin)
AB="${ABS[$((SLURM_ARRAY_TASK_ID-1))]}"


SCRIPT="/data/biol-micro-genomics/pkannan/scripts/treewas_no_chunks.r"
Rscript "${SCRIPT}" "${AB}"
