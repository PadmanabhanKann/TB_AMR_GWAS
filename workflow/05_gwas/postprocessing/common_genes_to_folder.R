#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default=NULL) {
  if (!(flag %in% args)) return(default)
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  args[i + 1]
}

K <- as.integer(get_arg("--k", "5"))
outdir <- get_arg("--outdir", "common_gene_analysis")
pattern <- "gene_hits\\.sorted\\.txt$"

if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# Find all gene_hits.sorted.txt anywhere under current dir
files <- list.files(".", pattern=pattern, recursive=TRUE, full.names=TRUE)

# Keep only those that end with /gene_hits.sorted.txt (allow optional leading ./)
files <- files[grepl("(^|\\./).+/gene_hits\\.sorted\\.txt$", files)]

if (length(files) == 0) {
  stop("No gene_hits.sorted.txt files found under current directory. Are you in cond_gwas2/?")
}

# Drug = directory name containing the file
drug <- sub("^\\./", "", files)
drug <- sub("/gene_hits\\.sorted\\.txt$", "", drug)

# Read gene lists
gene_lists <- lapply(files, function(f) {
  df <- read.delim(f, header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names=FALSE)
  if (!("gene" %in% colnames(df))) stop("Missing 'gene' column in: ", f)
  unique(df$gene)
})
names(gene_lists) <- drug

all_genes <- sort(unique(unlist(gene_lists)))

presence <- matrix(
  0L,
  nrow = length(all_genes),
  ncol = length(gene_lists),
  dimnames = list(all_genes, names(gene_lists))
)

for (d in names(gene_lists)) presence[gene_lists[[d]], d] <- 1L

n_drugs <- rowSums(presence)
gene_counts <- data.frame(gene=rownames(presence), n_drugs=as.integer(n_drugs), stringsAsFactors=FALSE)
gene_counts <- gene_counts[order(-gene_counts$n_drugs, gene_counts$gene), ]

genes_all <- gene_counts$gene[gene_counts$n_drugs == ncol(presence)]
genes_geK <- gene_counts$gene[gene_counts$n_drugs >= K]

write.table(presence, file=file.path(outdir, "gene_presence_matrix.tsv"),
            sep="\t", quote=FALSE, col.names=NA)

write.table(gene_counts, file=file.path(outdir, "gene_counts.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

writeLines(genes_all, file.path(outdir, "genes_common_all.txt"))
writeLines(genes_geK, file.path(outdir, paste0("genes_common_ge", K, ".txt")))

summary_file <- file.path(outdir, "summary.txt")
cat("Common gene analysis summary\n",
    "============================\n",
    "Total drugs analysed: ", ncol(presence), "\n",
    "Total unique genes: ", length(all_genes), "\n",
    "Genes common to ALL drugs: ", length(genes_all), "\n",
    "Genes present in >= ", K, " drugs: ", length(genes_geK), "\n\n",
    "Files used:\n",
    paste(files, collapse="\n"), "\n\n",
    "Top 20 genes by occurrence:\n",
    sep="", file=summary_file)

write.table(head(gene_counts, 20), file=summary_file, sep="\t", quote=FALSE,
            row.names=FALSE, col.names=TRUE, append=TRUE)

cat("✔ Results written to:", outdir, "\n")
