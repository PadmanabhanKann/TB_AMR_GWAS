#!/bin/bash

# --- Usage help ---
usage() {
  echo "Usage: $0 -i <alignment_file>"
  echo "  -i    Path to input alignment file (e.g., .phy, .aln, .fasta)"
  exit 1
}

# --- Parse arguments ---
while getopts ":i:" opt; do
  case $opt in
    i) INPUT="$OPTARG" ;;
    *) usage ;;
  esac
done

# --- Check input ---
if [ -z "$INPUT" ]; then
  echo "❌ Error: Input file is required."
  usage
fi

# --- Run IQ-TREE with ModelFinder Plus and auto threads ---
iqtree -s "$INPUT" -m MFP -T AUTO -bb 1000 -pre "$(basename "$INPUT" .phy)_iqtree"

# --- Done ---
echo "✅ IQ-TREE run complete. Output files prefixed with: $(basename "$INPUT" .phy)_iqtree"
