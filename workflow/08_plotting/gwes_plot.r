#!/usr/bin/env Rscript

# Help output
print_usage <- function() {
  cat("Usage: ./gwes_plot.r -i input_links_file -n number_of_assemblies [other options]\n")
  cat("Options:\n")
  cat("  -i, --input            Input links file (required)\n")
  cat("  -n, --number           Number of assemblies in the analysis (required)\n")
  cat("  -o, --output           Output plot path (default: ./pangwes_plot.svg)\n")
  cat("  --device               Output device: svg|png|pdf (default: svg)\n")
  cat("  -l, --ld-dist          LD distance cutoff (default: 50000)\n")
  cat("  --no-deps              Generate a basic plot without any R package dependencies (slower)\n")
  cat("  --highlight            Comma-separated list of unitig IDs to highlight (e.g., 102028,102030)\n")
  cat("  --annot                TSV/CSV with columns: U1 U2 label   (annotate red points with arrow + text)\n")
}

# Parse args
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  options <- list(
    input = NULL,
    output = "pangwes_plot.svg",
    device = "svg",
    number = NULL,
    basic_plot_only = FALSE,
    ld_dist = 50000,
    highlight = NULL,
    annot = NULL
  )
  i <- 1
  while (i <= length(args)) {
    if (args[i] == "-i" || args[i] == "--input") {
      options$input <- args[i + 1]; i <- i + 1
    } else if (args[i] == "-o" || args[i] == "--output") {
      options$output <- args[i + 1]; i <- i + 1
    } else if (args[i] == "--device") {
      options$device <- tolower(args[i + 1]); i <- i + 1
    } else if (args[i] == "-n" || args[i] == "--number") {
      options$number <- as.integer(args[i + 1])
      if (is.na(options$number)) stop("Number of assemblies must be an integer", call. = FALSE)
      i <- i + 1
    } else if (args[i] == "--no-deps") {
      options$basic_plot_only <- TRUE
    } else if (args[i] == "-l" || args[i] == "--ld-dist") {
      options$ld_dist <- as.integer(args[i + 1]); i <- i + 1
    } else if (args[i] == "--highlight") {
      options$highlight <- args[i + 1]; i <- i + 1
    } else if (args[i] == "--annot") {
      options$annot <- args[i + 1]; i <- i + 1
    } else if (args[i] == "-h" || args[i] == "--help") {
      print_usage(); quit(status = 0)
    } else {
      stop(paste("Unknown input option:", args[i], ". Use --help for options."), call. = FALSE)
    }
    i <- i + 1
  }
  if (is.null(options$input) || is.null(options$number)) {
    stop("Input links file and the number of assemblies must be provided, use --help for options.", call. = FALSE)
  }
  if (!options$device %in% c("svg","png","pdf")) {
    stop("Unsupported --device. Use one of: svg, png, pdf", call. = FALSE)
  }
  return(options)
}

# Parse 9-column input
parse_input <- function(input) {
  cat("Parsing input...\n")
  if (ncol(input) != 9) stop("Input links file must have 9 columns", call. = FALSE)
  if (nrow(input) < 100) stop("Input has less than 100 rows, cannot generate plot!", call. = FALSE)
  colnames(input) <- c("U1","U2","sep","ARACNE","MI","count","d_M2","d_min","d_max")

  n_K_disconnected <- length(which(input$sep == -1))
  if (length(n_K_disconnected) == 0) n_K_disconnected <- 0

  cat("Filtering connected links...")
  connected <- input[input$sep != -1, ]
  if (nrow(connected) == 0) stop("Cannot generate plot, no connected links found in the input links file!")
  connected$rsd <- (sqrt(connected$d_M2 / (num_assemblies - 1))) / connected$sep

  idx_to_keep <- which(connected$count >= ceiling(0.05 * num_assemblies) & connected$rsd <= 1)
  if (length(idx_to_keep) < 100) stop("Cannot generate plot, less than 100 links remain after filtering!")

  filtered_data <- connected[idx_to_keep, ]
  q <- stats::quantile(filtered_data$MI, c(0.25, 0.75))
  cat("Computing outlier thresholds...")
  outlier_thresh <- q[[2]] + 1.5 * (q[[2]] - q[[1]])
  extreme_outlier_thresh  <- q[[2]] + 3 * (q[[2]] - q[[1]])

  return(list(
    filtered_data = filtered_data,
    n_K_disconnected = n_K_disconnected,
    n_K_connected = floor(nrow(connected) / 1000),
    outlier_thresh = outlier_thresh,
    extreme_outlier_thresh = extreme_outlier_thresh
  ))
}

# Check and install packages
check_and_install_packages <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      cat(paste("Package", pkg, "is not installed. Installing now...\n"))
      install.packages(pkg, repos = "http://cran.us.r-project.org")
    }
  }
}

# Robust saver (svg -> svglite if available; else base svg)
save_plot <- function(p, path, device = "svg", width_in = 25, height_in = 12.5, pointsize = 12) {
  dev <- tolower(device)
  if (dev == "svg") {
    if (requireNamespace("svglite", quietly = TRUE)) {
      svglite::svglite(path, width = width_in, height = height_in)
      print(p); grDevices::dev.off()
    } else {
      message("svglite not found; using base grDevices::svg()")
      grDevices::svg(path, width = width_in, height = height_in, pointsize = pointsize)
      print(p); grDevices::dev.off()
    }
  } else if (dev == "pdf") {
    grDevices::pdf(path, width = width_in, height = height_in, useDingbats = FALSE)
    print(p); grDevices::dev.off()
  } else if (dev == "png") {
    # use ~133 dpi
    grDevices::png(path, width = as.integer(133 * width_in), height = as.integer(133 * height_in), pointsize = pointsize)
    print(p); grDevices::dev.off()
  } else {
    stop("Unknown device: ", dev)
  }
}

# ---- Main ----
opts <- parse_args()

highlight_unitigs <- NULL
if (!is.null(opts$highlight)) {
  highlight_unitigs <- as.numeric(strsplit(opts$highlight, ",")[[1]])
  cat("Highlighting unitigs:", paste(highlight_unitigs, collapse = ", "), "\n")
}

input_path <- opts$input
if (!file.exists(input_path)) stop("Input file does not exist", call. = FALSE)

output_path <- opts$output
if (file.exists(output_path)) stop("Plot already exists, please delete before running again", call. = FALSE)

num_assemblies <- opts$number
basic_plot_only <- opts$basic_plot_only

# Echo primary parameters
cat("Input file:", input_path, "\n")
cat("Output file:", output_path, "\n")
cat("Number of assemblies:", num_assemblies, "\n")
cat("Device:", opts$device, "\n")

## Packages (skip in no-deps mode)
if (!basic_plot_only) {
  required_packages <- c("data.table", "ggplot2", "ggthemes", "hexbin", "ggrastr", "ggrepel")
  check_and_install_packages(required_packages)
}

## Read input
cat("Reading input file...")
time_reading_start <- proc.time()
if (basic_plot_only) {
  input <- utils::read.csv(input_path, header = FALSE, sep = " ")
} else {
  library(data.table)
  library(ggplot2)
  library(ggrastr)
  library(ggthemes)
  library(hexbin)
  library(ggrepel)
  library(grid)
  input <- data.table::fread(input_path, header = FALSE, sep = " ")
}
time_reading_end <- proc.time()
time_reading <- (time_reading_end - time_reading_start)[[3]]
cat(paste0("Done in ", time_reading, "s\n"))

input_list <- parse_input(input)
connected <- input_list$filtered_data
outlier_threshold <- input_list$outlier_thresh
extreme_outlier_threshold <- input_list$extreme_outlier_thresh

# Build highlight rows
highlight_direct <- highlight_indirect <- NULL
if (!is.null(highlight_unitigs) && length(highlight_unitigs) >= 1) {
  in_pair <- (connected$U1 %in% highlight_unitigs) | (connected$U2 %in% highlight_unitigs)
  highlight_direct   <- connected[in_pair & connected$ARACNE == 1, , drop = FALSE]
  highlight_indirect <- connected[in_pair & connected$ARACNE == 0, , drop = FALSE]
  cat("Highlight stats:\n")
  cat("  Direct (ARACNE=1):   ", nrow(highlight_direct),   "edges\n")
  cat("  Indirect (ARACNE=0): ", nrow(highlight_indirect), "edges\n")
}

# ---- Build annotation points (optional; only for ARACNE=1/red) ----
annot_points <- NULL
if (!is.null(opts$annot)) {
  cat("Reading annotation file:", opts$annot, "\n")
  if (!file.exists(opts$annot)) stop("Annotation file not found: ", opts$annot, call. = FALSE)

  read_any <- function(p) {
    if (requireNamespace("data.table", quietly = TRUE)) {
      data.table::fread(p, header = TRUE)
    } else {
      utils::read.csv(p, header = TRUE, sep = "", check.names = FALSE)
    }
  }
  ann <- read_any(opts$annot)

  need_cols <- c("U1","U2","label")
  if (!all(need_cols %in% names(ann))) {
    stop("Annotation file must have columns: U1, U2, label", call. = FALSE)
  }
  ann$U1 <- as.numeric(ann$U1)
  ann$U2 <- as.numeric(ann$U2)

  red <- connected[connected$ARACNE == 1, ]

  key <- function(u1, u2) paste(pmin(u1,u2), pmax(u1,u2), sep=":")
  red$key <- key(red$U1, red$U2)
  ann$key <- key(ann$U1, ann$U2)

  annot_points <- merge(
    ann[, c("U1","U2","label","key")],
    red[, c("U1","U2","sep","MI","key")],
    by = "key", all.x = TRUE, suffixes = c("_ann","_red")
  )
  annot_points <- annot_points[!is.na(annot_points$sep), ]
  if (nrow(annot_points) > 0) {
    cat("Will annotate", nrow(annot_points), "red points.\n")
  } else {
    cat("No matching red points found to annotate.\n")
  }
}

# Fixed plot parameters
plot_width <- 2400   # px (for base path scaling); ~25 in @96 dpi
plot_height <- 1200
color_ld <- "black"
color_outlier <- grDevices::rgb(165, 0, 38, maxColorValue = 255)
color_extreme_outlier <- grDevices::rgb(215, 48, 39, maxColorValue = 255)
color_direct <- grDevices::rgb(0, 115, 190, maxColorValue = 255)    # blue
color_indirect <- grDevices::rgb(192, 192, 192, maxColorValue = 255) # gray

## Generate plot
cat("Generating plot... \n")
time_plotting_start <- proc.time()

if (basic_plot_only) {
  # ----- Base graphics -----
  plot_pointsize <- 16
  plot_symbol <- 19
  cex_direct <- 0.2
  cex_indirect <- 0.1
  cex_legend <- 1.2

  min_mi <- min(connected[, 5], na.rm = TRUE)
  max_mi <- max(connected[, 5], na.rm = TRUE)
  max_distance <- max(connected[, 3], na.rm = TRUE)
  exponent <- round(log10(max_distance)) - 1

  # Open device by flag
  if (opts$device == "svg") {
    grDevices::svg(output_path, width = plot_width/96, height = plot_height/96, pointsize = plot_pointsize)
  } else if (opts$device == "png") {
    grDevices::png(output_path, width = plot_width, height = plot_height, pointsize = plot_pointsize)
  } else if (opts$device == "pdf") {
    grDevices::pdf(output_path, width = plot_width/96, height = plot_height/96, pointsize = plot_pointsize, useDingbats = FALSE)
  }

  plot(connected[!connected[, 4], 3], connected[!connected[, 4], 5], col = color_indirect, type = "p", pch = plot_symbol, cex = cex_indirect,
       xlim = c(0, max_distance), ylim = c(min_mi, max_mi), xaxs = "i", yaxs = "i",
       xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
  points(connected[as.logical(connected[, 4]), 3], connected[as.logical(connected[, 4]), 5], col = color_direct, pch = plot_symbol, cex = cex_direct)
  axis(1, at = seq(0, max_distance, 10^exponent), tick = FALSE, labels = seq(0, max_distance / 10^exponent), line = -0.8)
  axis(2, at = seq(0.05, 1, 0.05), labels = FALSE, tcl = -0.5)
  axis(2, at = seq(0.1, 1, 0.1), labels = seq(0.1, 1, 0.1), las = 1, tcl = -0.5)
  title(xlab = "Distance between unitigs (bp)", line = 1.2)
  title(xlab = substitute(x10^exp, list(exp = exponent)), line = 1.4, adj = 1)
  title(ylab = "Mutual Information", line = 2.5)

  # Outlier thresholds
  if (outlier_threshold > 0) {
    segments(0, outlier_threshold, max_distance, outlier_threshold, col = color_outlier, lty = 2, lwd = 2)
    text(0, outlier_threshold, "*", col = color_outlier, pos = 2, offset = 0.2, cex = 1, xpd = NA)
  }
  if (extreme_outlier_threshold > 0) {
    segments(0, extreme_outlier_threshold, max_distance, extreme_outlier_threshold, col = color_extreme_outlier, lty = 2, lwd = 2)
    text(0, extreme_outlier_threshold, "**", col = color_extreme_outlier, pos = 2, offset = 0.2, cex = 1, xpd = NA)
  }

  # LD line
  if (opts$ld_dist > 0) { segments(opts$ld_dist, min_mi, opts$ld_dist, 1, col = color_ld, lty = 2, lwd = 2) }

  # Highlights
  if (!is.null(highlight_direct) && nrow(highlight_direct) > 0) {
    points(highlight_direct$sep, highlight_direct$MI, col = "red", pch = 19, cex = 1.6)
  }
  if (!is.null(highlight_indirect) && nrow(highlight_indirect) > 0) {
    points(highlight_indirect$sep, highlight_indirect$MI, col = "orange", pch = 17, cex = 1.4)
  }

  # Annotations (arrows + text)
  if (!is.null(annot_points) && nrow(annot_points) > 0) {
    dx <- 0.02 * max_distance
    dy <- 0.02 * (max_mi - min_mi)
    for (k in seq_len(nrow(annot_points))) {
      x0 <- annot_points$sep[k]; y0 <- annot_points$MI[k]
      x1 <- x0 + dx; y1 <- y0 + dy
      arrows(x0, y0, x1, y1, length = 0.08, angle = 25, col = "red", lwd = 1)
      text(x1, y1, labels = annot_points$label[k], pos = 4, cex = 0.9)
    }
  }

  # Legend
  legend("topright",
         legend = c("Indirect", "Direct", "rRNA (ARACNE = 1)", "rRNA (ARACNE ≠ 1)"),
         col    = c(color_indirect, color_direct, "red", "orange"),
         pch    = c(19, 19, 19, 17),
         pt.cex = c(0.8, 0.8, 1.2, 1.1),
         bty = "n", cex = cex_legend)

  grDevices::dev.off()

} else {
  # ---- ggplot2 path (legend for grey/blue + bigger red/orange) ----
  bg_size <- 0.25       # background dot size (rasterized)
  hi_size_red <- 2.2    # highlight circle (red)
  hi_size_orange <- 2.0 # highlight triangle (orange)

  # Split background
  bg_indirect <- connected[connected$ARACNE != 1, ]
  bg_direct   <- connected[connected$ARACNE == 1, ]

  # Legend-enabled highlight DF
  legend_df <- NULL
  if (!is.null(highlight_direct) && nrow(highlight_direct) > 0) {
    hd <- highlight_direct; hd$rrna_type <- "rRNA (ARACNE = 1)"
    legend_df <- rbind(legend_df, hd)
  }
  if (!is.null(highlight_indirect) && nrow(highlight_indirect) > 0) {
    hi <- highlight_indirect; hi$rrna_type <- "rRNA (ARACNE ≠ 1)"
    legend_df <- rbind(legend_df, hi)
  }

  p1 <- ggplot()

  # Rasterized background WITH legend (map color to string labels)
  p1 <- p1 +
    ggrastr::geom_point_rast(
      data = bg_indirect,
      aes(x = sep, y = MI, color = "Indirect"),
      size = bg_size, alpha = 0.25, shape = 16, show.legend = TRUE
    ) +
    ggrastr::geom_point_rast(
      data = bg_direct,
      aes(x = sep, y = MI, color = "Direct"),
      size = bg_size, alpha = 0.6, shape = 16, show.legend = TRUE
    )

  # Vector highlights (bigger)
  if (!is.null(legend_df) && nrow(legend_df) > 0) {
    p1 <- p1 +
      geom_point(
        data = subset(legend_df, rrna_type == "rRNA (ARACNE = 1)"),
        aes(x = sep, y = MI, color = rrna_type, shape = rrna_type),
        size = hi_size_red, stroke = 0.2, show.legend = TRUE
      ) +
      geom_point(
        data = subset(legend_df, rrna_type == "rRNA (ARACNE ≠ 1)"),
        aes(x = sep, y = MI, color = rrna_type, shape = rrna_type),
        size = hi_size_orange, stroke = 0.2, show.legend = TRUE
      )
  }

  # Annotations (arrow + label)
  if (!is.null(annot_points) && nrow(annot_points) > 0) {
    p1 <- p1 + geom_point(
      data = annot_points, aes(x = sep, y = MI),
      color = "red", shape = 16, size = hi_size_red + 0.4, stroke = 0
    )
    p1 <- p1 + ggrepel::geom_label_repel(
      data = annot_points,
      aes(x = sep, y = MI, label = label),
      min.segment.length = 0,
      seed = 1,
      box.padding = 0.25,
      point.padding = 0.3,
      label.size = 0.15,
      label.r = unit(0.1, "lines"),
      size = 2.8,
      segment.size = 0.3,
      segment.alpha = 0.9,
      arrow = arrow(length = unit(0.02, "npc")),
      fill = "white",
      color = "black",
      show.legend = FALSE
    )
  }

  # Unified legend + theming
  p1 <- p1 +
    scale_color_manual(
      name = NULL,
      values = c(
        "Indirect" = color_indirect,
        "Direct"   = color_direct,
        "rRNA (ARACNE = 1)" = "red",
        "rRNA (ARACNE ≠ 1)" = "orange"
      ),
      breaks = c("Indirect", "Direct", "rRNA (ARACNE = 1)", "rRNA (ARACNE ≠ 1)")
    ) +
    scale_shape_manual(
      values = c("rRNA (ARACNE = 1)" = 16, "rRNA (ARACNE ≠ 1)" = 17),
      breaks = c("rRNA (ARACNE = 1)", "rRNA (ARACNE ≠ 1)")
    ) +
    guides(shape = "none") +
    scale_x_continuous(breaks = function(x) pretty(x, n = 10)) +
    geom_hline(yintercept = outlier_threshold,         col = '#a50026', linetype = 'dashed', linewidth = 0.3) +
    geom_hline(yintercept = extreme_outlier_threshold, col = '#d73027', linetype = 'dashed', linewidth = 0.3) +
    geom_vline(xintercept = c(opts$ld_dist), col = 'black', linetype = 'dashed', linewidth = 0.3) +
    ggthemes::theme_clean(base_size = 8) +
    xlab("Graph distance") + ylab("Mutual information") +
    theme(
      legend.position = "right",
      legend.background = element_blank(),
      legend.box.background = element_blank()
    )

  # Save (~25x12.5 in is your original px at 96dpi)
  save_plot(p1, output_path, device = opts$device,
            width_in = plot_width/96, height_in = plot_height/96, pointsize = 12)
}

time_plotting_end <- proc.time()
time_plotting <- (time_plotting_end - time_plotting_start)[[3]]
cat(paste0("All done! Plotting time: ", time_plotting, "s\n"))
