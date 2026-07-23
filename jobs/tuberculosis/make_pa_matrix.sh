#!/bin/bash
#SBATCH --job-name=tb_pa_matrix
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=10:00:00
#SBATCH --partition=short,medium,long
#SBATCH --error=/data/biol-micro-genomics/pkannan/logs/pa_matrix_%j.err
#SBATCH --output=/data/biol-micro-genomics/pkannan/logs/pa_matrix_%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=padmanabhan21@iisertvm.ac.in
set -euo pipefail

TB_DIR="/data/biol-micro-genomics/pkannan/TB"
PYSEER2_DIR="${TB_DIR}/pyseer2"
UNITIG_GZ="${PYSEER2_DIR}/unitig.pyseer.gz"
OUT_BASE="${TB_DIR}/pa_matrix"

# activate your env (edit only if needed)
source activate /data/biol-micro-genomics/pkannan/env/pyseer

mkdir -p "${OUT_BASE}"

if [[ ! -f "${UNITIG_GZ}" ]]; then
  echo "ERROR: Cannot find ${UNITIG_GZ}"
  exit 1
fi

echo "Using unitigs: ${UNITIG_GZ}"
echo "Output base:  ${OUT_BASE}"
echo

for ab_dir in "${PYSEER2_DIR}"/*; do
  [[ -d "${ab_dir}" ]] || continue
  ab_name="$(basename "${ab_dir}")"

  # skip obvious non-drug dirs if present
  case "${ab_name}" in
    ref|logs) continue ;;
  esac

  echo "=============================="
  echo "Processing: ${ab_name}"
  echo "=============================="

  cd "${ab_dir}"

  # inputs expected per antibiotic
  if [[ ! -f "significant_kmers.txt" ]]; then
    echo "  [SKIP] missing significant_kmers.txt in ${ab_dir}"
    continue
  fi
  if [[ ! -f "annotated_flag.txt" ]]; then
    echo "  [SKIP] missing annotated_flag.txt in ${ab_dir}"
    continue
  fi

  # 1) get list of significant kmers (col1; skip header)
  tail -n +2 significant_kmers.txt | cut -f1 > sig_kmers.list

  # 2) extract presence lines for those kmers from unitig.pyseer.gz (order preserved)
  awk -F' \\| ' 'NR==FNR{rank[$1]=++c; next}
    ($1 in rank){print rank[$1]"\t"$0; seen[$1]=1}
    END{
      for(k in rank) if(!(k in seen)) print k > "missing_sig_kmers.list"
    }' \
    sig_kmers.list <(zcat "${UNITIG_GZ}") \
  | sort -n -k1,1 | cut -f2- > sigkmer_presence.ordered.pyseer

  gzip -f sigkmer_presence.ordered.pyseer

  # 3) unitig format -> presence/absence matrix (rows=kmers, cols=samples)
  python3 - << 'PY'
import pandas as pd
import gzip

inp = "sigkmer_presence.ordered.pyseer.gz"
out = "kmer_presence_ab.csv"

unitigs = []
sample_data = {}
all_samples = set()

with gzip.open(inp, 'rt') as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        parts = line.split(' | ')
        if len(parts) != 2:
            continue
        unitig = parts[0]
        samples = [s.split(':')[0] for s in parts[1].split()]
        unitigs.append(unitig)
        sample_data[unitig] = samples
        all_samples.update(samples)

all_samples = sorted(all_samples)

matrix = []
for u in unitigs:
    present = set(sample_data[u])
    matrix.append([1 if s in present else 0 for s in all_samples])

df = pd.DataFrame(matrix, index=unitigs, columns=all_samples)
df.to_csv(out, index=True)
print(f"Wrote {out} shape={df.shape}")
PY

  # 4) replace kmer IDs using annotated_flag.txt (col1=kmer, col8=annotation)
  python3 - << 'PY'
import csv

ANNOT = "annotated_flag.txt"
INCSV = "kmer_presence_ab.csv"
OUTCSV = "kmer_presence_ab_gene.csv"

k2g = {}
with open(ANNOT, newline="") as f:
    rdr = csv.reader(f, delimiter="\t")
    for row in rdr:
        if not row or len(row) < 8:
            continue
        kmer = row[0].strip()
        annot = row[7].strip()

        genes, seen = [], set()
        for tok in (t for t in annot.split(";") if t):
            if ":" in tok:   # drop contig:coords tokens
                continue
            if tok not in seen:
                seen.add(tok)
                genes.append(tok)

        if genes:
            k2g[kmer] = "|".join(genes)

with open(INCSV, newline="") as fi, open(OUTCSV, "w", newline="") as fo:
    rdr = csv.reader(fi)
    w = csv.writer(fo)
    first = True
    for row in rdr:
        if first:
            first = False
            w.writerow(row)
            continue
        if not row:
            w.writerow(row)
            continue
        kmer = row[0].strip()
        row[0] = k2g.get(kmer, kmer)
        w.writerow(row)

print(f"Wrote {OUTCSV}")
PY

# ---------------------------
# Step 5: transpose + append phenotype if present
# Final output: kmer_presence_ab_gene_T.csv
# ---------------------------
python3 - << 'PY'
import pandas as pd
import os

PA_IN  = "kmer_presence_ab_gene.csv"
PA_OUT = "kmer_presence_ab_gene_T.csv"
PHENO  = "phenotype.csv"

# transpose
df = pd.read_csv(PA_IN, index_col=0)
df_t = df.T
df_t.index.name = "id"
df_t.reset_index(inplace=True)

# if phenotype.csv exists, merge it
if os.path.isfile(PHENO):
    ph = pd.read_csv(PHENO)

    # enforce types
    ph["id"] = ph["id"].astype(str)
    df_t["id"] = df_t["id"].astype(str)

    # inner join ensures perfect alignment
    merged = ph.merge(df_t, on="id", how="inner")

    # reorder: id, phenotype, rest
    cols = ["id", "phenotype"] + [c for c in merged.columns if c not in ("id", "phenotype")]
    merged = merged[cols]

    merged.to_csv(PA_OUT, index=False)
    print(f"Wrote {PA_OUT} with phenotype column, shape={merged.shape}")

else:
    df_t.to_csv(PA_OUT, index=False)
    print(f"Wrote {PA_OUT} (no phenotype found), shape={df_t.shape}")
PY
  # 6) save to p/a_matrix/<antibiotic>/
  dest="${OUT_BASE}/${ab_name}"
  mkdir -p "${dest}"

  cp -f sig_kmers.list "${dest}/"
  cp -f missing_sig_kmers.list "${dest}/" 2>/dev/null || true
  cp -f sigkmer_presence.ordered.pyseer.gz "${dest}/"
  cp -f kmer_presence_ab.csv "${dest}/"
  cp -f kmer_presence_ab_gene.csv "${dest}/"
  cp -f kmer_presence_ab_gene_T.csv "${dest}/"

  echo "  [DONE] ${ab_name} -> ${dest}"
  echo
done

echo "All antibiotics processed."
