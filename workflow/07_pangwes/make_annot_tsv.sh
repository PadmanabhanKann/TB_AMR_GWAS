awk -v RS='>' 'NR>1{split($1,h," "); id=h[1]; seq=""; for(i=2;i<=NF;i++) seq=seq $i; gsub(/\r/,"",seq); print seq "\t" id}' selected_unitigs.fa \
| LC_ALL=C sort \
| join -t $'\t' -1 1 -2 1 - <(LC_ALL=C sort annotated_flag.txt) \
| awk -F'\t' '
  function firstgene(s, n,i,f,a){ n=split(s,a,","); for(i=1;i<=n;i++){ split(a[i],f,";"); if(f[2]!="") return f[2] } return "NA" }
  { map[$2]=firstgene($3) } 
  END{
    print "U1\tU2\tlabel"
    while((getline l < "sig_pairs.txt")>0){
      if(l ~ /^awk\b/ || l ~ /^[[:space:]]*$/) continue
      n=split(l,f,"\t"); if(n<2 || f[1]=="U1") continue
      u1=f[1]; u2=f[2]
      g1=(u1 in map ? map[u1] : "NA")
      g2=(u2 in map ? map[u2] : "NA")
      printf "%s\t%s\t%s–%s\n", u1, u2, g1, g2
    }
  }' > annot.tsv
