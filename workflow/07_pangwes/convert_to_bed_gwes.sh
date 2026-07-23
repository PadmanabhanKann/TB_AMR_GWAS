awk -F'\t' '
{
  seq = $1;
  coords = $2;

  n = split(coords, arr, ",");

  for (i = 1; i <= n; i++) {
    split(arr[i], a, ";");
    split(a[1], b, ":");
    split(b[2], c, "-");

    start = c[1];   # no -1
    end   = c[2];

    print b[1], start, end, seq;
  }
}
' OFS="\t" annotated_flag_pangwes.txt > pangwes.bed
