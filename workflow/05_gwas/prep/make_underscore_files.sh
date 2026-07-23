#!/bin/bash

# Input files
PHENO="phenotype.txt"
KMERS="gent.pyseer"
SIM="phylogeny_K.tsv"

# Output files
PHENO_OUT="phenotype_underscore.txt"
KMERS_OUT="gent_underscore.pyseer"
SIM_OUT="phylogeny_K_underscore.tsv"

echo "🔧 Converting phenotype file..."
awk -F'\t' '{gsub(/\./, "_", $1); print $1 "\t" $2}' "$PHENO" > "$PHENO_OUT"

echo "🔧 Converting k-mer file header..."
awk 'NR==1 {for(i=1;i<=NF;i++) gsub(/\./,"_",$i)} {print}' "$KMERS" > "$KMERS_OUT"

echo "🔧 Converting similarity matrix..."
{
  head -n 1 "$SIM" | tr '\t' '\n' | sed 's/\./_/g' | paste -sd '\t'
  tail -n +2 "$SIM" | awk -F'\t' '{
    gsub(/\./,"_",$1);
    for (i=2; i<=NF; i++) gsub(/\./,"_",$i);
    print
  }'
} > "$SIM_OUT"

echo "✅ Done. Output files:"
echo "  - $PHENO_OUT"
echo "  - $KMERS_OUT"
echo "  - $SIM_OUT"
