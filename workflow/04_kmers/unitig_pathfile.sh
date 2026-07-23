#!/bin/bash
#to create strain.txt containing all paths to the fast files
# Usage: ./unitig_pathfile.sh output_file path1 path2 ...

OUTPUT_FILE=$1
shift

# Clear the output file if it exists
> "$OUTPUT_FILE"

# Loop through all input directories
for DIR in "$@"; do
    find "$DIR" -type f -name "*.fna" ! -path "*/prokka_output/*" | while read -r file; do
        realpath "$file"
    done
done >> "$OUTPUT_FILE"
