#!/bin/bash
set -euo pipefail

BASE="/data/biol-micro-genomics/pkannan/TB/cond_gwas2"
MATRIX="/data/biol-micro-genomics/pkannan/TB/treewas/kmer_matrix_T.csv"

command -v gawk >/dev/null 2>&1 || { echo "ERROR: gawk not found"; exit 1; }
[ -s "$MATRIX" ] || { echo "ERROR: missing matrix: $MATRIX"; exit 1; }

for AB_DIR in "$BASE"/*/; do
  AB="$(basename "${AB_DIR%/}")"
  LIST="$AB_DIR/${AB}_covariate_kmers.list"
  OUT="$AB_DIR/${AB}_subset.covariates.tsv"

  [ -s "$LIST" ] || { echo "skip $AB (missing list $LIST)"; continue; }

  echo "==== $AB ===="
  echo "list: $LIST"
  echo "out : $OUT"

  gawk -v LISTFILE="$LIST" '
    BEGIN{
      FS=","; OFS="\t";
      while((getline line < LISTFILE) > 0){
        gsub(/\r/,"",line);
        if(line!="") want[line]=1;
      }
      close(LISTFILE);
    }

    NR==1{
      # header row: decide which columns to keep
      keepN=0;

      # always keep sample column (col 1)
      keep[++keepN]=1;
      name[keepN]=$1;

      for(i=2;i<=NF;i++){
        if($i in want){
          keep[++keepN]=i;
          name[keepN]=$i;
          found[$i]=1;
        }
      }

      # print header (tab-separated)
      for(j=1;j<=keepN;j++){
        printf "%s%s", name[j], (j<keepN?OFS:ORS);
      }

      # report how many kmers were found in header
      # (to stderr so it doesnt pollute output file)
      nWant=0; for(x in want) nWant++;
      nFound=0; for(x in found) nFound++;
      printf "Selected %d/%d kmer columns (+sample)\n", nFound, nWant > "/dev/stderr";

      next
    }

    {
      # print selected columns for each sample row
      for(j=1;j<=keepN;j++){
        printf "%s%s", $(keep[j]), (j<keepN?OFS:ORS);
      }
    }
  ' "$MATRIX" > "$OUT"

  echo
done
