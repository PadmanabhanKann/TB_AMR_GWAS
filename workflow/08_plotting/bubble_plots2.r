#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
})

base_dir <- getwd()

## If esxQ is missing, label top N rows instead
top_n_to_label_if_no_esxQ <- 5

## central folder for all bubbleplots
all_bubble_dir <- file.path(base_dir, "bubbleplots")
dir.create(all_bubble_dir, showWarnings = FALSE)

ab_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)

for (ab_dir in ab_dirs) {

  ab_name <- basename(ab_dir)

  ## Use sorted file
  gene_hits_file <- file.path(ab_dir, "gene_hits.sorted.txt")
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

  ## =========================
  ## LABEL LOGIC:
  ## - If esxQ exists: label everything from top down to esxQ (inclusive)
  ## - If esxQ missing: label top N rows
  ## =========================
  gene_hits$label <- NA

  idx_esxQ <- which(gene_hits$gene == "esxQ")[1]

  if (!is.na(idx_esxQ)) {
    gene_hits$label[1:idx_esxQ] <- gene_hits$gene[1:idx_esxQ]
  } else {
    top_n <- min(top_n_to_label_if_no_esxQ, nrow(gene_hits))
    gene_hits$label[1:top_n] <- gene_hits$gene[1:top_n]
    warning(paste("esxQ not found in", ab_name, "- labeling top", top_n, "rows instead."))
  }

  ## =========================
  ## PLOT
  ## =========================
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
