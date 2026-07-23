python3 - << 'PY'
import csv

ANNOT = "annotated_flag.txt"   # TSV: kmer in col1, annotation in col8
INCSV = "kmer_presence_ab.csv"      # CSV: kmer in col1
OUTCSV = "kmer_presence_ab_gene.csv"

# Build kmer -> gene(s) map
k2g = {}
with open(ANNOT, newline="") as f:
    rdr = csv.reader(f, delimiter="\t")
    for row in rdr:
        if not row or len(row) < 8:
            continue
        kmer = row[0].strip()
        annot = row[7].strip()
        # annotation looks like: NZ_...:start-end;gene;gene;...
        genes, seen = [], set()
        for tok in (t for t in annot.split(";") if t):
            if ":" in tok:              # skip the coordinate token
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
        if not row:
            w.writerow(row); continue
        if first and row[0].lower() in ("kmer","variant"):
            row[0] = "gene"
            first = False
        else:
            kmer = row[0].strip()
            gene = k2g.get(kmer)
            row[0] = gene if gene else kmer   # replace; keep kmer if no mapping
            # If you prefer NA instead, use: row[0] = gene if gene else "NA"
        w.writerow(row)

print(f"Done: wrote {OUTCSV}")
PY
