#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
})

base_dir <- getwd()

## central folder for all bubbleplots
all_bubble_dir <- file.path(base_dir, "bubbleplots")
dir.create(all_bubble_dir, showWarnings = FALSE)

ab_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)

for (ab_dir in ab_dirs) {

  ab_name <- basename(ab_dir)
  gene_hits_file <- file.path(ab_dir, "gene_hits.txt")

  if (!file.exists(gene_hits_file)) next

  outdir <- file.path(ab_dir, "bubbleplot")
  dir.create(outdir, showWarnings = FALSE)

  message("Processing: ", ab_name)

  gene_hits <- read.table(
    gene_hits_file,
    header = TRUE,
    stringsAsFactors = FALSE,
    sep = "\t",
    quote = "",
    comment.char = ""
  )

  ## Drop any gene containing "cds"
  gene_hits <- gene_hits[!grepl("cds", gene_hits$gene, ignore.case = TRUE), ]

  req_cols <- c("gene", "avg_beta", "maxp", "avg_maf", "hits")
  missing <- setdiff(req_cols, colnames(gene_hits))
  if (length(missing) > 0) {
    stop(paste("Missing columns in", ab_name, ":", paste(missing, collapse = ", ")))
  }

  ## Drug-specific labels (only these will be annotated)
  label_genes <- "esxQ"  # in all

  if (ab_name == "clofazimine") {
    label_genes <- c(label_genes, "metK")
  } else if (ab_name == "moxifloxacin") {
    label_genes <- c(label_genes, "gyrA", "gyrA1")
  } else if (ab_name == "ethambutol") {
    label_genes <- c(label_genes, "embB", "katG", "rpoB")
  } else if (ab_name == "linezolid") {
    label_genes <- c(label_genes, "rpoB")
  }

  ## Create label column: only marked genes get text, others are NA
  gene_hits$label <- ifelse(gene_hits$gene %in% label_genes, gene_hits$gene, NA)

  p <- ggplot(
    gene_hits,
    aes(x = avg_beta, y = maxp, colour = avg_maf, size = hits)
  ) +
    geom_point(alpha = 0.5) +
    geom_text_repel(
      aes(label = label),
      show.legend = FALSE,
      colour = "black",
      max.overlaps = Inf
    ) +
    scale_size("Number of k-mers", range = c(1, 10)) +
    scale_colour_gradient("Average MAF") +
    theme_bw(base_size = 14) +
    ggtitle(paste(ab_name, "resistance")) +
    xlab("Average effect size") +
    ylab("Maximum -log10(p-value)")

  png_file <- file.path(outdir, paste0(ab_name, "_bubbleplot.png"))
  pdf_file <- file.path(outdir, paste0(ab_name, "_bubbleplot.pdf"))

  ggsave(png_file, p, width = 10, height = 7, dpi = 300)
  ggsave(pdf_file, p, width = 10, height = 7)

  ## copy PNG to central folder
  file.copy(
    from = png_file,
    to   = file.path(all_bubble_dir, paste0(ab_name, "_bubbleplot.png")),
    overwrite = TRUE
  )
}
