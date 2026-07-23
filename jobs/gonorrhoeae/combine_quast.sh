#!/bin/bash

# Go to the directory containing batch_* folders
cd /data/biol-micro-genomics/pkannan/Azithro_N.gon/quast_results

# Output file
OUTPUT="transposed_all.tsv"

# Flag for header control
first_file=true

# Loop through all batch folders
for dir in batch_*; do
    file="$dir/transposed_report.tsv"
    if [ -f "$file" ]; then
        if $first_file; then
            cat "$file" > "$OUTPUT"
            first_file=false
        else
            # Skip header (first line)
            tail -n +2 "$file" >> "$OUTPUT"
        fi
    else
        echo "Warning: $file not found!"
    fi
done

echo "Combined file created: $OUTPUT"
