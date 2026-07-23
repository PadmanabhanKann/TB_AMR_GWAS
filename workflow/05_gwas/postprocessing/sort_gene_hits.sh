(head -n1 gene.hits && tail -n +2 gene.hits | sort -t$'\t' -k3,3gr) > sorted_gene.hits
