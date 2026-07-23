#!/bin/bash
# Recorded invocation used to generate the annotated GWES plot for the
# gonorrhoeae 23S rRNA / azithromycin figure (n=1029 genomes).
./gwes_plot.r -i out.ud_sgg_0_based -n 1029 --highlight "$(paste -sd, rrna_unitigs.txt)" --annot annot.tsv -o pangwes_plot_annot.svg
