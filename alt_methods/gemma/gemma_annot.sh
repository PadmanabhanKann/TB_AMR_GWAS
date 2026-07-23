#!/bin/bash
#SBATCH --job-name=gemma_annot_all
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --partition=short,medium,long
#SBATCH --output=/data/biol-micro-genomics/pkannan/logs/gemma_annot_all.%j.out
#SBATCH --error=/data/biol-micro-genomics/pkannan/logs/gemma_annot_all.%j.err
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=padmanabhan21@iisertvm.ac.in 

set -euo pipefail

# ---------------- paths ----------------
BASE="/data/biol-micro-genomics/pkannan/TB/gemma"
SCRIPT="/data/biol-micro-genomics/pkannan/software/gemma-annotate/bin/gemma-annotate"
BED5="${BASE}/ref/ref.genes.bed5"
RES="${BASE}/results"

mkdir -p "${BASE}/logs"

[[ -x "$SCRIPT" ]] || { echo "ERROR: not executable: $SCRIPT"; exit 1; }
[[ -f "$BED5" ]]   || { echo "ERROR: BED5 missing: $BED5"; exit 1; }

echo "SCRIPT: $SCRIPT"
echo "BED5:   $BED5"
echo "BASE:   $BASE"
echo

# ---------- maxp thresholds (OLIGONUCLEOTIDE; SNP-level) ----------
# maxp = -log10(p)  →  p = 10^(-maxp)
declare -A MAXP
MAXP[amikacin]="7.38"
MAXP[kanamycin]="7.38"

MAXP[bedaquiline]="7.34"

MAXP[clofazimine]="7.34"
MAXP[levofloxacin]="7.34"

MAXP[delamanid]="7.36"
MAXP[ethionamide]="7.36"

MAXP[ethambutol]="7.32"

MAXP[isoniazid]="7.39"
MAXP[rifabutin]="7.39"

MAXP[linezolid]="7.33"

MAXP[moxifloxacin]="7.31"

MAXP[rifampicin]="7.37"

# ---------------- main loop ----------------
for drug in "${!MAXP[@]}"; do

  assoc="${RES}/${drug}/output/${drug}_lmm.assoc.txt"
  outjson="${RES}/${drug}/${drug}.annotated.json"
  outtsv="${RES}/${drug}/${drug}.annotated.tsv"
  maxp="${MAXP[$drug]}"

  if [[ ! -f "$assoc" ]]; then
    echo "SKIP: ${drug} (missing assoc: $assoc)"
    continue
  fi

  # compute p threshold from maxp
  pthr="$(python - <<PY
import math
print("{:.12e}".format(10**(-float("${maxp}"))))
PY
)"

  echo "========================================"
  echo "Drug: $drug"
  echo "maxp (-log10p): $maxp"
  echo "pthr (p):       $pthr"
  echo "Assoc: $assoc"
  echo "JSON:  $outjson"
  echo "TSV:   $outtsv"
  echo "========================================"

  # ---- GEMMA annotate ----
  "$SCRIPT" \
    --eval "p_lrt <= ${pthr}" \
    --value p_lrt \
    --bed "$BED5" \
    "$assoc" > "$outjson"

  # ---- JSON → TSV ----
  python - <<PY
import json
drug="${drug}"
inp="${outjson}"
out="${outtsv}"

recs=json.load(open(inp))
with open(out,"w") as f:
    f.write("drug\tchr\tpos\tsnp\tp_lrt\tlocus_tag\tgene\n")
    for r in recs:
        locus=r.get("anno","NA")
        gene=r.get("gene","NA")
        f.write(
            f"{drug}\t{r.get('chr')}\t{r.get('pos')}\t"
            f"{r.get('name')}\t{r.get('p_lrt')}\t"
            f"{locus}\t{gene}\n"
        )

print(f"{drug}: {len(recs)} hits → {out}")
PY

done

echo "DONE."
