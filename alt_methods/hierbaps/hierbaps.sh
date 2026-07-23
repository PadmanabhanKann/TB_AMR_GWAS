#!/bin/bash
#SBATCH --job-name=hierbaps
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=short,medium,long
#SBATCH --time=0-10:00:00
#SBATCH --error=/data/biol-micro-genomics/pkannan/logs/hierbaps_%j.err
#SBATCH --output=/data/biol-micro-genomics/pkannan/logs/hierbaps_%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=padmanabhan21@iisertvm.ac.in

set -euo pipefail

# Load conda module (ARC style)
module purge
module load Anaconda3/2023.09-0

# Activate your conda env that already contains R 4.3.3
source activate /data/biol-micro-genomics/pkannan/env/pyseer

# Go to alignment directory
cd /data/biol-micro-genomics/pkannan/TB/superdca/hierbaps

# Run hierBAPS and write Level-1 clusters only
Rscript --vanilla - <<'RSCRIPT'
library(ape)
library(rhierbaps)

aln_file <- "core_gene_alignment.aln"
dna <- read.dna(aln_file, format = "fasta")

res <- hierBAPS(dna, max.depth = 3, n.pops = 20)
df <- res$partition.df

sample_id <- rownames(df)
if (is.null(sample_id) || all(sample_id == "")) sample_id <- names(dna)

lvl1_col <- grep("^level\\.?1$|^level_?1$|^level\\.1$|^level1$",
                 colnames(df), ignore.case = TRUE, value = TRUE)
if (length(lvl1_col) == 0) {
  stop("Could not find Level 1 column. Columns are: ",
       paste(colnames(df), collapse = ", "))
}
lvl1_col <- lvl1_col[1]

out <- data.frame(
  lane_id = sample_id,
  hierBAPS_level1 = df[[lvl1_col]],
  stringsAsFactors = FALSE
)

write.csv(out, "hierBAPS_level1.csv", row.names = FALSE, quote = FALSE)
print(head(out))
cat("Wrote hierBAPS_level1.csv\n")
RSCRIPT
