#!/bin/bash

# Usage: ./unitig-caller.sh -i strain_list.txt -o output_dir

set -euo pipefail

# Parse args
while getopts ":i:o:" opt; do
  case $opt in
    i) INPUT_FILE="$OPTARG" ;;
    o) OUTPUT_DIR="$OPTARG" ;;
    *) echo "Usage: $0 -i strain_list.txt -o output_dir" >&2; exit 1 ;;
  esac
done

# Validate input
if [[ -z "${INPUT_FILE:-}" || -z "${OUTPUT_DIR:-}" ]]; then
  echo "Usage: $0 -i strain_list.txt -o output_dir"
  exit 1
fi

# Create output folder
mkdir -p "$OUTPUT_DIR"

# Run unitig-caller – outputs will use default names in that folder
unitig-caller --call --pyseer  --refs "$INPUT_FILE" --out "$OUTPUT_DIR"
