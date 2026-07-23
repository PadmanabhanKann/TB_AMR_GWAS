DIR=/data/biol-micro-genomics/pkannan/sra/fastq

ls -1 $DIR/*_1.fastq | sort > R1.list
sed 's/_1\.fastq$/_2.fastq/' R1.list > R2.list
paste -d, R1.list R2.list > /data/biol-micro-genomics/pkannan/sra/spades/fixed_paired.txt
