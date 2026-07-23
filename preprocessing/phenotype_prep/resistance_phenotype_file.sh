#!/bin/bash

# Usage:
# ./resistance_phenotype_file.sh -r /path/to/res_dir -s /path/to/sus_dir -o phenotype.txt

while getopts "r:s:o:" opt; do
  case $opt in
    r) RESISTANT_DIR="$OPTARG" ;;
    s) SUSCEPTIBLE_DIR="$OPTARG" ;;
    o) OUTPUT_FILE="$OPTARG" ;;
    \?) echo "Usage: $0 -r resistant_dir -s susceptible_dir -o output_file" >&2; exit 1 ;;
  esac
done

# Check input
if [ ! -d "$RESISTANT_DIR" ] || [ ! -d "$SUSCEPTIBLE_DIR" ]; then
  echo "Error: Both -r and -s must be valid directories." >&2
  exit 1
fi

> "$OUTPUT_FILE"  # Clear output

# Resistant: label 1
for file in "$RESISTANT_DIR"/*.fna; do
  fname=$(basename "$file")
  genome_id="${fname%%.fna}"
  echo -e "${genome_id}\t1" >> "$OUTPUT_FILE"
done

# Susceptible: label 0
for file in "$SUSCEPTIBLE_DIR"/*.fna; do
  fname=$(basename "$file")
  genome_id="${fname%%.fna}"
  echo -e "${genome_id}\t0" >> "$OUTPUT_FILE"
done

echo "Phenotype file created at $OUTPUT_FILE"
