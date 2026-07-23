#!/bin/bash
#SBATCH --job-name=post_pyseer_check
#SBATCH --nodes=1
#SBATCH --partition=short,medium,long
#SBATCH --time=0-03:00:00
#SBATCH --error=/data/biol-micro-genomics/pkannan/logs/post_pyseer_check_%A.err
#SBATCH --output=/data/biol-micro-genomics/pkannan/logs/post_pyseer_check_%A.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=padmanabhan21@iisertvm.ac.in

set -euo pipefail

# -------------------------
# Config
# -------------------------
BASE="/data/biol-micro-genomics/pkannan/TB/cond_gwas2"
ENV_PATH="/data/biol-micro-genomics/pkannan/env/pyseer"

REFS="$BASE/references.txt"
REF_FNA="$(awk 'NR==1{print $1}' "$REFS")"
REF_GFF="$(awk 'NR==1{print $2}' "$REFS")"

COUNT_PATTERNS="$BASE/count_patterns.py"
REMAPPER_PY="$BASE/phandango_remapper.py"
SUMMARISE="$BASE/summarise_annotations.py"

# Only do these three (as per your plan)
ABLIST=(ethambutol linezolid moxifloxacin)

# Prefix used in the conditional run (must match your pyseer output)
RUN_PREFIX_SUFFIX="kmers_cond_plus_clof"

# -------------------------
# Activate conda env
# -------------------------
module purge || true
module load Anaconda3/2023.09-0 || true
source activate "$ENV_PATH"

# -------------------------
# Checks
# -------------------------
die(){ echo "ERROR: $*" >&2; exit 1; }

[[ -d "$BASE" ]] || die "BASE not found: $BASE"
[[ -s "$REFS" ]] || die "Missing/empty: $REFS"
[[ -s "$REF_FNA" ]] || die "Reference FNA missing (from references.txt): $REF_FNA"
[[ -s "$REF_GFF" ]] || die "Reference GFF missing (from references.txt): $REF_GFF"

[[ -s "$COUNT_PATTERNS" ]] || die "Missing: $COUNT_PATTERNS"
[[ -s "$REMAPPER_PY"    ]] || die "Missing: $REMAPPER_PY"
[[ -s "$SUMMARISE"      ]] || die "Missing: $SUMMARISE"

command -v phandango_mapper >/dev/null 2>&1 || die "phandango_mapper not found in PATH"
command -v annotate_hits_pyseer >/dev/null 2>&1 || die "annotate_hits_pyseer not found in PATH"

echo "Env active: $(python -V 2>&1)"

# -------------------------
# Main: iterate target antibiotics, process AB/check/
# -------------------------
cd "$BASE"

for AB in "${ABLIST[@]}"; do
  AB_DIR="${BASE}/${AB}"
  CHECK_DIR="${AB_DIR}/check"
  [[ -d "$CHECK_DIR" ]] || { echo "Skip $AB (no check dir): $CHECK_DIR"; continue; }

  echo "==== $AB (check) ===="
  cd "$CHECK_DIR" || { echo "Cannot cd into $CHECK_DIR"; continue; }

  RES="${AB}_${RUN_PREFIX_SUFFIX}_results.txt"
  PATS="${AB}_${RUN_PREFIX_SUFFIX}_patterns.txt"

  [[ -s "$RES"  ]] || { echo "Missing $RES";  cd "$BASE"; continue; }
  [[ -s "$PATS" ]] || { echo "Missing $PATS"; cd "$BASE"; continue; }

  # Threshold from patterns
  THR="$(python "$COUNT_PATTERNS" "$PATS" | awk '
    {for(i=1;i<=NF;i++) if($i ~ /^[0-9.]+E[-+][0-9]+$/){print $i; exit}}
  ')"
  [[ -n "${THR:-}" ]] || { echo "Could not parse threshold for $AB"; cd "$BASE"; continue; }
  echo "Threshold: $THR"

  # Filter significant kmers (p-value assumed column 4)
  SIG="significant_kmers.txt"
  cat <(head -n 1 "$RES") <(awk -v thr="$THR" 'NR>1 && $4 < thr {print $0}' "$RES") > "$SIG"

  mkdir -p plot

  # A) SIGNIFICANT kmers → plot + remap
  cp -f "$SIG" plot/
  ( cd plot && phandango_mapper "$(basename "$SIG")" "$REF_FNA" "${AB}_significant.plot" )

  python "$REMAPPER_PY" \
    --txt "plot/${AB}_significant.plot" \
    --gff "$REF_GFF" \
    --output "plot/${AB}_significant.remap.plot"

  # B) ALL kmers (full results) → plot + remap
  cp -f "$RES" plot/
  ( cd plot && phandango_mapper "$(basename "$RES")" "$REF_FNA" "${AB}_all.plot" )

  python "$REMAPPER_PY" \
    --txt "plot/${AB}_all.plot" \
    --gff "$REF_GFF" \
    --output "plot/${AB}_all.remap.plot"

  # C) Annotate + summarise (SIGNIFICANT)
  annotate_hits_pyseer "$SIG" "$REFS" annotated_kmers.txt --feature-type rRNA
  python "$SUMMARISE" annotated_kmers.txt > gene_hits.txt
  ( head -n 1 gene_hits.txt && tail -n +2 gene_hits.txt | sort -k3,3nr ) > gene_hits.sorted.txt

  # D) Annotate + summarise (ALL results)
  annotate_hits_pyseer "$RES" "$REFS" annotated_allkmers.txt --feature-type rRNA
  python "$SUMMARISE" annotated_allkmers.txt > gene_all_hits.txt
  ( head -n 1 gene_all_hits.txt && tail -n +2 gene_all_hits.txt | sort -k3,3nr ) > gene_all_hits.sorted.txt

  cd "$BASE"
done

echo "Done check-postprocessing"
