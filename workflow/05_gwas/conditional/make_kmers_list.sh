#!/bin/bash
# ------------------------------------------------------------
# Script: extract_covariate_kmers.sh
#
# Purpose:
#   For each antibiotic, extract k-mers corresponding to
#   specified covariate genes / CDS IDs from annotated_kmers.txt
#
# Directory structure expected:
#   cond_gwas2/
#     ├── amikacin/annotated_kmers.txt
#     ├── bedaquiline/annotated_kmers.txt
#     ├── ethambutol/annotated_kmers.txt
#     └── ...
#
# Output:
#   <antibiotic>/<antibiotic>_covariate_kmers.list
#
# ------------------------------------------------------------

set -euo pipefail

BASE="/data/biol-micro-genomics/pkannan/TB/cond_gwas2"
MAP="$BASE/covariates_map.tsv"

# -----------------------------
# HOW TO CREATE covariates_map.tsv
# -----------------------------
# The file MUST be TAB-separated with exactly 2 columns:
#
#   column 1: antibiotic folder name
#   column 2: comma-separated list of covariate genes / CDS IDs
#
# Example:
#
# amikacin      cds-NP_218247.1
# ethambutol   embB,rpoB,katG
# delamanid    cds-NP_215829.1,cds-NP_215909
# rifampicin   rpoB,katG
#
# Notes:
# - Use COMMAS between multiple covariates
# - NO spaces around commas
# - Gene names are matched as ;gene; in annotated_kmers.txt
# - NP_XXXX.X IDs are matched anywhere in the annotation field
# -----------------------------

[ -s "$MAP" ] || { echo "ERROR: missing $MAP"; exit 1; }

echo "Using covariate map: $MAP"
echo

while IFS=$'\t' read -r AB COVS; do
  # skip empty or malformed lines
  [ -z "${AB:-}" ] && continue
  [ -z "${COVS:-}" ] && continue

  ANN="$BASE/$AB/annotated_kmers.txt"
  OUT="$BASE/$AB/${AB}_covariate_kmers.list"

  if [ ! -s "$ANN" ]; then
    echo "WARN: $ANN not found, skipping $AB"
    continue
  fi

  echo "Processing $AB → $COVS"

  gene_re=""
  np_re=""

  # split comma-separated covariates
  IFS=',' read -ra ITEMS <<< "$COVS"
  for c in "${ITEMS[@]}"; do
    c="${c//[[:space:]]/}"

    # NP_XXXX.X identifiers
    if [[ "$c" =~ NP_[0-9]+\.[0-9]+$ ]]; then
      np_re="${np_re:+$np_re|}$c"
    else
      # gene names (matched as ;gene;)
      c_lc="$(echo "$c" | tr '[:upper:]' '[:lower:]')"
      gene_re="${gene_re:+$gene_re|}$c_lc"
    fi
  done

  awk -F'\t' -v GRE="$gene_re" -v NPRE="$np_re" '
    BEGIN{IGNORECASE=1}
    {
      ann=$NF
      hit=0
      if (GRE  != "" && ann ~ ";(" GRE ");") hit=1
      if (!hit && NPRE != "" && ann ~ "(" NPRE ")") hit=1
      if (hit) print $1
    }
  ' "$ANN" | sort -u > "$OUT"

  echo "  → $(wc -l < "$OUT") kmers written"
  echo

done < "$MAP"

echo "Done extracting covariate k-mers"
