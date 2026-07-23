# For 16S rRNA (first copy)
efetch -db nuccore -id CP001050.1 -seq_start 1266912 -seq_stop 1268440 -strand 2 -format fasta > NCCP11945_16S_rRNA1.fasta

# For 23S rRNA (first copy)
efetch -db nuccore -id CP001050.1 -seq_start 1263410 -seq_stop 1266299 -strand 2 -format fasta > NCCP11945_23S_rRNA1.fasta

# For 5S rRNA (first copy)
efetch -db nuccore -id CP001050.1 -seq_start 1263199 -seq_stop 1263311 -strand 2 -format fasta > NCCP11945_5S_rRNA1.fasta

# 5S rRNA (NGK_rrna5s2)
efetch -db nuccore -id CP001050.1 -seq_start 1620687 -seq_stop 1620799 -strand 2 -format fasta > NCCP11945_rrna5s2.fasta

# 23S rRNA (NGK_rrna23s2)
efetch -db nuccore -id CP001050.1 -seq_start 1620898 -seq_stop 1623787 -strand 2 -format fasta > NCCP11945_rrna23s2.fasta

# 16S rRNA (NGK_rrna16s2)
efetch -db nuccore -id CP001050.1 -seq_start 1624400 -seq_stop 1625928 -strand 2 -format fasta > NCCP11945_rrna16s2.fasta


# 5S rRNA (NGK_rrna5s3)
efetch -db nuccore -id CP001050.1 -seq_start 1724586 -seq_stop 1724698 -strand 2 -format fasta > NCCP11945_rrna5s3.fasta

# 23S rRNA (NGK_rrna23s3)
efetch -db nuccore -id CP001050.1 -seq_start 1724797 -seq_stop 1727686 -strand 2 -format fasta > NCCP11945_rrna23s3.fasta

# 16S rRNA (NGK_rrna16s3)
efetch -db nuccore -id CP001050.1 -seq_start 1728299 -seq_stop 1729827 -strand 2 -format fasta > NCCP11945_rrna16s3.fasta

# 16S rRNA (NGK_rrna16s4)
efetch -db nuccore -id CP001050.1 -seq_start 1954347 -seq_stop 1955875 -strand 1 -format fasta > NCCP11945_rrna16s4.fasta

# 23S rRNA (NGK_rrna23s4)
efetch -db nuccore -id CP001050.1 -seq_start 1956488 -seq_stop 1959377 -strand 1 -format fasta > NCCP11945_rrna23s4.fasta

# 5S rRNA (NGK_rrna5s4)
efetch -db nuccore -id CP001050.1 -seq_start 1959476 -seq_stop 1959588 -strand 1 -format fasta > NCCP11945_rrna5s4.fasta
