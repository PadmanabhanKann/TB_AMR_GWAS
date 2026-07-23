#!/bin/bash

# Usage: ./run_panaroo.sh -i /path/to/gff_folder -o /path/to/output_folder

# Default values
INPUT_FOLDER=""
OUTPUT_FOLDER=""

# Parse command-line options
while getopts "i:o:" opt; do
  case $opt in
    i)
      INPUT_FOLDER="$OPTARG"
      ;;
    o)
      OUTPUT_FOLDER="$OPTARG"
      ;;
    *)
      echo "Usage: $0 -i <input_folder> -o <output_folder>"
      exit 1
      ;;
  esac
done

# Validate inputs
if [[ -z "$INPUT_FOLDER" || -z "$OUTPUT_FOLDER" ]]; then
  echo "Error: Both -i and -o parameters are required."
  echo "Usage: $0 -i <input_folder_with_gff_files> -o <output_folder>"
  exit 1
fi

# Check if Panaroo is installed
if ! command -v panaroo &> /dev/null; then
  echo "Error: Panaroo is not installed or not in PATH."
  exit 1
fi

# Create output folder if it doesn't exist
mkdir -p "$OUTPUT_FOLDER"

# Collect GFF files
GFF_FILES=("$INPUT_FOLDER"/*.gff)
if [ ${#GFF_FILES[@]} -eq 0 ]; then
  echo "Error: No .gff files found in $INPUT_FOLDER"
  exit 1
fi

# Run Panaroo
echo "Running Panaroo with the following parameters:"
echo "  Input folder : $INPUT_FOLDER"
echo "  Output folder: $OUTPUT_FOLDER"

panaroo -i "${GFF_FILES[@]}" \
        -o "$OUTPUT_FOLDER" \
        -a core \
        --clean-mode strict \
        --aligner mafft \

echo "✅ Panaroo finished. Results in: $OUTPUT_FOLDER"
