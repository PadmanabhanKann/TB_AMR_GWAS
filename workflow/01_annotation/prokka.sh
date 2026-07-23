#!/bin/bash

# --- Parse arguments ---
while getopts ":i:o:g:" opt; do
  case ${opt} in
    i ) INPUT_DIR=$OPTARG ;;
    o ) OUTPUT_DIR=$OPTARG ;;
    g ) GFF_DIR=$OPTARG ;;
    \? ) echo "Usage: $0 -i input_folder -o prokka_output_folder -g gff_output_folder"; exit 1 ;;
  esac
done

# --- Validate input ---
if [[ -z "$INPUT_DIR" || -z "$OUTPUT_DIR" || -z "$GFF_DIR" ]]; then
  echo "❌ All of -i (input), -o (prokka output), and -g (gff output) must be provided."
  echo "Usage: $0 -i input_folder -o prokka_output_folder -g gff_output_folder"
  exit 1
fi

if [[ ! -d "$INPUT_DIR" ]]; then
  echo "❌ Input folder not found: $INPUT_DIR"
  exit 1
fi

mkdir -p "$OUTPUT_DIR"
mkdir -p "$GFF_DIR"

# --- Run Prokka on all .fna files ---
for file in "$INPUT_DIR"/*.fna; do
  base=$(basename "$file" .fna)
  echo "▶️ Running Prokka on: $base"

  prokka --outdir "${OUTPUT_DIR}/${base}" --prefix "$base" "$file"
done

# --- Collect all .gff files ---
find "$OUTPUT_DIR" -name "*.gff" -exec cp {} "$GFF_DIR/" \;

echo "✅ Done. All GFF files copied to: $GFF_DIR"
