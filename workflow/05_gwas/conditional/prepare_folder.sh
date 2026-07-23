SRC="/data/biol-micro-genomics/pkannan/TB/cond_gwas"
DST="/data/biol-micro-genomics/pkannan/TB/cond_gwas2"

for ab_dir in "$SRC"/*; do
  ab=$(basename "$ab_dir")

  # skip non-directories
  [ -d "$ab_dir" ] || continue

  echo "Processing $ab"
  mkdir -p "$DST/$ab"

  # copy required files if they exist
  cp -f "$ab_dir/annotated_kmers.txt"        "$DST/$ab/" 2>/dev/null || echo "  missing annotated_kmers.txt"
  cp -f "$ab_dir/$ab.tsv"                    "$DST/$ab/" 2>/dev/null || echo "  missing $ab.tsv"
  cp -f "$ab_dir/gene_hits.txt"              "$DST/$ab/" 2>/dev/null || echo "  missing gene_hits.txt"
  cp -f "$ab_dir/gene_hits.sorted.txt"       "$DST/$ab/" 2>/dev/null || echo "  missing gene_hits.sorted.txt"
done
