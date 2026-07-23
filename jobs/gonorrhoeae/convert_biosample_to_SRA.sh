echo "BioSample,Run" > biosample_to_run.csv
for biosdir in downloads/SAMN*; do
  bios=$(basename "$biosdir")
  # handle SRR/ERR/DRR just in case
  for run in "$biosdir"/SRR* "$biosdir"/ERR* "$biosdir"/DRR*; do
    [ -d "$run" ] || continue
    echo "$bios,$(basename "$run")" >> biosample_to_run.csv
  done
done

awk -F, 'NR==FNR {              # pass 1: mapping
  if (NR>1) map[$1]=(map[$1]?map[$1]","$2:$2);  # BioSample -> comma list of runs
  next
}
FNR==1 { print "taxon\ttrait"; next }  # phenotype header
{
  bios=$1; mic=$2
  n=split(map[bios], arr, ",")
  for (i=1;i<=n;i++) if (arr[i]!="") print arr[i]"\t"mic
}' biosample_to_run.csv Azi_MIC.csv > phenotype_srr.tsv
