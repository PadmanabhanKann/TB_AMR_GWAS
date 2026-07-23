#!/bin/bash

# Usage: ./download_sra.sh -i <biosample_list.txt> [-o <output_dir>]
# Downloads FASTQ files for every BioSample accession in the input list.

# Disable exit on error so the script continues even if a command fails
set +e

# Defaults match the original trial run; override per-organism with -i/-o
input_file="/data/biol-micro-genomics/pkannan/sra/trial.txt"
output_base="/data/biol-micro-genomics/pkannan/sra/fastq"

while getopts "i:o:" opt; do
  case $opt in
    i) input_file="$OPTARG" ;;
    o) output_base="$OPTARG" ;;
    *) echo "Usage: $0 -i <biosample_list.txt> [-o <output_dir>]" >&2; exit 1 ;;
  esac
done

# Ensure enaBrowserTools is in your PATH
export PATH="$PATH:/data/biol-micro-genomics/pkannan/sra/enaBrowserTools/python3"

# Create base output directory if it doesn't exist
mkdir -p "$output_base"

# Loop through each BioSample accession
while IFS= read -r biosample; do
    # Skip empty lines
    [[ -z "$biosample" ]] && continue

    echo "=========================================="
    echo "Processing BioSample: $biosample"

    # Create subdirectory for this BioSample
    output_dir="$output_base/$biosample"
    mkdir -p "$output_dir"

    # Fetch associated Run accessions
    run_accessions=$(esearch -db biosample -query "$biosample" \
        | elink -target sra \
        | efetch -format runinfo \
        | cut -d',' -f1 \
        | grep -E "SRR|ERR|DRR")

    # Check if any runs were found
    if [ -z "$run_accessions" ]; then
        echo "Warning: No run accessions found for $biosample"
        continue
    fi

    echo "Found runs: $run_accessions"

    # Download FASTQ files for each Run accession
    for run in $run_accessions; do
        echo "Downloading FASTQ for Run: $run"
        enaDataGet -f fastq -d "$output_dir" "$run"
        if [ $? -ne 0 ]; then
            echo "Warning: Download failed for $run, skipping..."
        else
            echo "Completed download for $run"
        fi
    done

done < "$input_file"

echo "=========================================="
echo "All BioSamples processed."
