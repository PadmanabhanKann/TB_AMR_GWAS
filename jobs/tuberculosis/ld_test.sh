#!/bin/bash
#SBATCH --job-name=ld_decay_debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=devel
#SBATCH --time=00:10:00
#SBATCH --error=/data/biol-micro-genomics/pkannan/logs/ld_decay_debug_%j.err
#SBATCH --output=/data/biol-micro-genomics/pkannan/logs/ld_decay_debug_%j.out

echo "========== BASIC JOB INFO =========="
hostname
date
whoami
pwd
echo

echo "========== ENV VARS =========="
echo "DATA=$DATA"
echo "PATH=$PATH"
echo

echo "========== MODULES =========="
module purge
module load Anaconda3/2023.09-0
module list
echo

echo "========== CONDA CHECK =========="
which conda || echo "conda NOT FOUND"
conda --version || echo "conda version failed"
echo

echo "========== TRY SOURCE ACTIVATE =========="
source activate /data/biol-micro-genomics/pkannan/env/python2 && echo "source activate worked" || echo "source activate FAILED"
echo

echo "========== PYTHON AFTER SOURCE =========="
which python || echo "python not found"
python --version || echo "python version failed"
echo

echo "========== TRY PROPER CONDA INIT =========="
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /data/biol-micro-genomics/pkannan/env/python2 && echo "conda activate worked" || echo "conda activate FAILED"
echo

echo "========== PYTHON AFTER CONDA ACTIVATE =========="
which python || echo "python not found"
python --version || echo "python version failed"
echo

echo "========== TEST PYTHON 2 PRINT =========="
python - <<EOF
print "THIS IS PYTHON 2 PRINT"
EOF
