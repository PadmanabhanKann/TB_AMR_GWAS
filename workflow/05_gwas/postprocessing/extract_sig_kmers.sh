#sort by  p value(ost signgicant(lowest)first)
{ head -n1 significant_kmers.txt.;   tail -n +2 significant_kmers.txt. | sort -t $'\t' -k4,4g; } > significant_kmers.sorted_by_p.tsv

#sort by beta value(highest first)
{ head -n1 significant_kmers.txt; tail -n +2 significant_kmers.txt | sort -t $'\t' -k5,5gr; } > significant_kmers.sorted_by_beta_desc.tsv




tail -n +2 significant_kmers.sorted_by_beta_desc.tsv | cut -f1 > sig_kmers.list

tail -n +2 significant_kmers.sorted_by_p.tsv | cut -f1 > sig_kmers.list

awk -F' \\| ' 'NR==FNR{rank[$1]=++c; next}
  ($1 in rank){print rank[$1]"\t"$0; seen[$1]=1}
  END{for(k in rank) if(!(k in seen)) print k > "missing_sig_kmers.list"}' \
  sig_kmers.list unitig.pyseer \
| sort -n -k1,1 | cut -f2- > sigkmer_presence.ordered.pyseer
