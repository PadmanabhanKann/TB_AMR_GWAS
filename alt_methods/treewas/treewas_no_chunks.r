#!/usr/bin/env Rscript

################################################################################
# treeWAS Analysis Script - Corrected Version
# 
# Features:
# - No chunking (processes all k-mers at once)
# - ML (Maximum Likelihood) reconstruction for both SNPs and phenotype
# - Bonferroni correction for multiple testing
# - Optimized for 100k k-mers, 3080 isolates
#
# Usage: Rscript treewas_corrected.r <antibiotic>
################################################################################

suppressPackageStartupMessages({
  library(treeWAS)
  library(ape)
  library(data.table)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript treewas_corrected.r <antibiotic>")
}
ab <- args[1]

# Define file paths
BASE <- "/data/biol-micro-genomics/pkannan/TB/treewas"
TREE_FILE <- file.path(BASE, "core_gene_alignment.aln_iqtree.treefile")
KMER_CSV  <- file.path(BASE, "kmer_filtered.csv")

# Determine phenotype file (try two possible names)
PH1 <- file.path(BASE, ab, "phenotype.csv")
PH2 <- file.path(BASE, ab, paste0("log_", ab, "_MIC.csv"))
PH_FILE <- if (file.exists(PH1)) PH1 else PH2

if (!file.exists(PH_FILE)) {
  stop("Phenotype file not found: ", PH_FILE)
}

# Create output directory
OUT_DIR <- file.path(BASE, ab, "ml_bonferroni")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Set analysis parameters
PVALUE <- 0.01
N_SNPS_SIM <- 10000  # 10x minimum; increase if memory allows (e.g., 50000-100000)
TESTS <- c("terminal", "simultaneous", "subsequent")

cat("========================================\n")
cat("treeWAS Analysis Configuration\n")
cat("========================================\n")
cat("Antibiotic:", ab, "\n")
cat("P-value threshold:", PVALUE, "\n")
cat("P-value correction: Bonferroni\n")
cat("Reconstruction method: ML (Maximum Likelihood)\n")
cat("Number of simulated SNPs:", N_SNPS_SIM, "\n")
cat("Tests:", paste(TESTS, collapse = ", "), "\n")
cat("Output directory:", OUT_DIR, "\n")
cat("========================================\n\n")

# -----------------------------
# Load tree
# -----------------------------
cat("Loading phylogenetic tree...\n")
tr <- read.tree(TREE_FILE)
cat("  Tree loaded:", length(tr$tip.label), "tips\n\n")

# -----------------------------
# Load phenotype
# -----------------------------
cat("Loading phenotype data...\n")
ph <- fread(PH_FILE)
ph$id <- as.character(ph$id)

# Find numeric phenotype column
num_cols <- setdiff(names(ph)[sapply(ph, is.numeric)], "id")
if (length(num_cols) == 0) {
  stop("No numeric phenotype column found in: ", PH_FILE)
}
valcol <- num_cols[1]

phen <- setNames(as.numeric(ph[[valcol]]), ph$id)
cat("  Phenotype loaded:", length(phen), "isolates\n")
cat("  Phenotype column:", valcol, "\n")
cat("  Phenotype range:", round(min(phen, na.rm = TRUE), 3), "to", 
    round(max(phen, na.rm = TRUE), 3), "\n\n")

# -----------------------------
# Read full k-mer matrix
# -----------------------------
cat("Loading k-mer matrix (this may take several minutes)...\n")

# Read header to get column names
hdr <- names(fread(KMER_CSV, nrows = 0))
idcol <- hdr[1]
kcols <- hdr[-1]
cat("  Total k-mers:", length(kcols), "\n")

# Read full matrix
dt <- fread(KMER_CSV, showProgress = TRUE)
km_ids <- as.character(dt[[idcol]])
cat("  K-mer matrix loaded:", nrow(dt), "isolates x", length(kcols), "k-mers\n\n")

# -----------------------------
# Find common isolates
# -----------------------------
cat("Finding common isolates across tree, k-mer matrix, and phenotype...\n")
common <- Reduce(intersect, list(tr$tip.label, km_ids, names(phen)))

if (length(common) == 0) {
  stop("ERROR: No common isolates found among tree, k-mer matrix, and phenotype")
}

cat("  Common isolates:", length(common), "\n")
cat("  Dropped from tree:", length(tr$tip.label) - length(common), "\n")
cat("  Dropped from k-mer matrix:", length(km_ids) - length(common), "\n")
cat("  Dropped from phenotype:", length(phen) - length(common), "\n\n")

# -----------------------------
# Subset data to common isolates
# -----------------------------
cat("Subsetting data to common isolates...\n")

# Subset tree
tr2 <- keep.tip(tr, common)

# Subset phenotype
phen2 <- phen[tr2$tip.label]

# Subset and convert k-mer matrix
idx <- match(tr2$tip.label, dt[[idcol]])
if (any(is.na(idx))) {
  stop("ERROR: ID mismatch between tree tips and k-mer CSV sample column")
}

dt[[idcol]] <- NULL
geno <- as.matrix(dt[idx, , drop = FALSE])
storage.mode(geno) <- "numeric"
rownames(geno) <- tr2$tip.label

cat("  Final dataset:", nrow(geno), "isolates x", ncol(geno), "k-mers\n")

# Clear large objects to free memory
rm(dt, ph, phen, tr)
gc(verbose = FALSE)
cat("  Memory cleared\n\n")

# -----------------------------
# Run treeWAS
# -----------------------------
cat("========================================\n")
cat("Running treeWAS analysis...\n")
cat("This may take 30-60 minutes or longer\n")
cat("========================================\n\n")

start_time <- Sys.time()

out <- treeWAS(
  snps = geno,
  phen = phen2,
  tree = tr2,
  phen.type = "continuous",
  phen.reconstruction = "ML",
  snps.reconstruction = "ML",
  snps.sim.reconstruction = "ML",
  n.snps.sim = N_SNPS_SIM,
  test = TESTS,
  p.value = PVALUE,
  p.value.correct = "bonf",  # Bonferroni correction (NOT "bonferroni" or FALSE)
  plot.tree = TRUE,
  plot.manhattan = TRUE,
  plot.null.dist = TRUE,
  plot.dist = FALSE,
  filename.plot = file.path(OUT_DIR, "treewas_plots.pdf")
 # Set seed for reproducibility
)

end_time <- Sys.time()
elapsed <- difftime(end_time, start_time, units = "mins")

cat("\n========================================\n")
cat("treeWAS analysis completed in", round(elapsed, 2), "minutes\n")
cat("========================================\n\n")

# -----------------------------
# Save significant k-mers using built-in write.treeWAS
# -----------------------------
cat("Saving significant k-mers using write.treeWAS()...\n")

# Use the built-in treeWAS function to write results
# This ensures proper formatting and avoids manual extraction errors
results_file <- file.path(OUT_DIR, "treewas_results")
write.treeWAS(out, filename = results_file)

cat("Results saved to:", paste0(results_file, ".csv"), "\n")

# Print summary of findings
cat("\nSummary of significant findings:\n")
for (tname in TESTS) {
  n_sig <- 0
  if (!is.null(out[[tname]]$sig.snps) && NROW(out[[tname]]$sig.snps) > 0) {
    n_sig <- NROW(out[[tname]]$sig.snps)
  }
  cat("  ", tname, "test:", n_sig, "significant k-mers\n")
}

# -----------------------------
# Save full treeWAS output object
# -----------------------------
cat("\nSaving full treeWAS output object...\n")
rds_file <- file.path(OUT_DIR, "treewas_output.rds")
saveRDS(out, rds_file)

# -----------------------------
# Generate summary report
# -----------------------------
cat("\nGenerating summary report...\n")

summary_file <- file.path(OUT_DIR, "analysis_summary.txt")
sink(summary_file)

cat("========================================\n")
cat("treeWAS Analysis Summary\n")
cat("========================================\n\n")

cat("Analysis Parameters:\n")
cat("  Antibiotic:", ab, "\n")
cat("  Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("  P-value threshold:", PVALUE, "\n")
cat("  P-value correction: Bonferroni\n")
cat("  Reconstruction method: ML\n")
cat("  Number of simulated SNPs:", N_SNPS_SIM, "\n")
cat("  Tests:", paste(TESTS, collapse = ", "), "\n\n")

cat("Input Data:\n")
cat("  Tree file:", TREE_FILE, "\n")
cat("  K-mer file:", KMER_CSV, "\n")
cat("  Phenotype file:", PH_FILE, "\n")
cat("  Phenotype column:", valcol, "\n\n")

cat("Dataset Dimensions:\n")
cat("  Number of isolates:", nrow(geno), "\n")
cat("  Number of k-mers:", ncol(geno), "\n\n")

cat("Phenotype Statistics:\n")
cat("  Min:", round(min(phen2, na.rm = TRUE), 3), "\n")
cat("  Max:", round(max(phen2, na.rm = TRUE), 3), "\n")
cat("  Mean:", round(mean(phen2, na.rm = TRUE), 3), "\n")
cat("  Median:", round(median(phen2, na.rm = TRUE), 3), "\n")
cat("  SD:", round(sd(phen2, na.rm = TRUE), 3), "\n\n")

cat("Results:\n")
for (tname in TESTS) {
  n_sig <- 0
  if (!is.null(out[[tname]]$sig.snps) && NROW(out[[tname]]$sig.snps) > 0) {
    n_sig <- NROW(out[[tname]]$sig.snps)
  }
  cat("  ", tname, "test:", n_sig, "significant k-mers\n")
}

cat("\nComputation Time:\n")
cat("  Elapsed time:", round(elapsed, 2), "minutes\n")

cat("\nOutput Files:\n")
cat("  Output directory:", OUT_DIR, "\n")
cat("  Results CSV:", paste0(results_file, ".csv"), "\n")
cat("  RDS object:", rds_file, "\n")
cat("  This summary:", summary_file, "\n")

cat("\n========================================\n")

sink()

cat("Summary report saved to:", summary_file, "\n")

# -----------------------------
# Final message
# -----------------------------
cat("\n========================================\n")
cat("Analysis Complete!\n")
cat("========================================\n")
cat("Output directory:", OUT_DIR, "\n")
cat("Results file:", paste0(results_file, ".csv"), "\n")
cat("Full output object:", rds_file, "\n")
cat("Summary report:", summary_file, "\n")
cat("========================================\n")
