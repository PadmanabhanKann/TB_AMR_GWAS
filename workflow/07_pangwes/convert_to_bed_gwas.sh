awk -F'\t' '
{
  seq = $1;
  coords = $NF;

  n = split(coords, arr, ",");

  for (i = 1; i <= n; i++) {
    split(arr[i], a, ";");
    split(a[1], b, ":");
    split(b[2], c, "-");

    start = c[1];   # keep as-is (1-based)
    end   = c[2];

    print b[1], start, end, seq;
  }
}
' OFS="\t" annotated_all_flag_GWAS > gwas.bed
