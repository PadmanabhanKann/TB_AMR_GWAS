#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(svglite)
library(optparse)

options(bitmapType = "cairo")

# ----------------------------
# Command line options
# ----------------------------
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Input LD file (PLINK .ld). Must be tab-delimited with header.", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output directory", metavar = "character"),
  make_option(c("-m", "--maxdist"), type = "double", default = NA,
              help = "Max distance (bp) to plot. Default = use full range in data.", metavar = "number"),
  make_option(c("-g", "--genome_len"), type = "double", default = NA,
              help = "Genome length (bp) for circular distance. If provided, uses circular distance; otherwise direct distance.",
              metavar = "number"),
  make_option(c("-b", "--bin_short"), type = "double", default = 1000,
              help = "Bin size (bp) for distances up to 100 kb (default 1000).", metavar = "number"),
  make_option(c("-B", "--bin_long"), type = "double", default = 10000,
              help = "Bin size (bp) for distances > 100 kb (default 10000).", metavar = "number")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input) || is.null(opt$output)) {
  print_help(opt_parser)
  stop("Both input file (-i) and output directory (-o) must be supplied.", call. = FALSE)
}

input_file <- opt$input
output_dir <- opt$output
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# Load data
# ----------------------------
message("Starting data loading...")
ld_data <- read.table(input_file, header = TRUE, sep = "", stringsAsFactors = FALSE, comment.char = "")

required_cols <- c("BP_A", "BP_B", "R2")
missing <- setdiff(required_cols, colnames(ld_data))
if (length(missing) > 0) {
  stop(paste("Input file is missing required columns:", paste(missing, collapse = ", ")), call. = FALSE)
}

# Ensure numeric
ld_data$BP_A <- as.numeric(ld_data$BP_A)
ld_data$BP_B <- as.numeric(ld_data$BP_B)
ld_data$R2   <- as.numeric(ld_data$R2)

# Drop rows with NA essentials
ld_data <- ld_data %>% filter(!is.na(BP_A) & !is.na(BP_B) & !is.na(R2))

# ----------------------------
# Distance calculation
# ----------------------------
message("Calculating distances...")
direct_dist <- abs(ld_data$BP_B - ld_data$BP_A)

if (!is.na(opt$genome_len) && opt$genome_len > 0) {
  GENOME_LEN <- opt$genome_len
  ld_data$distance <- pmin(direct_dist, GENOME_LEN - direct_dist)
  message(paste0("Using circular distance with genome length = ", GENOME_LEN, " bp"))
} else {
  ld_data$distance <- direct_dist
  message("Using direct (non-circular) distance")
}

if (!is.na(opt$maxdist) && opt$maxdist > 0) {
  ld_data <- ld_data %>% filter(distance <= opt$maxdist)
  message(paste0("Filtering to max distance = ", opt$maxdist, " bp"))
}

# ----------------------------
# Binning / aggregation
# ----------------------------
BIN_SHORT <- opt$bin_short
BIN_LONG  <- opt$bin_long

message("Aggregating data with adaptive binning...")
ld_data_summary <- ld_data %>%
  mutate(distance_interval = case_when(
    distance <= 100000 ~ floor(distance / BIN_SHORT) * BIN_SHORT,
    TRUE               ~ floor(distance / BIN_LONG)  * BIN_LONG
  )) %>%
  group_by(distance_interval) %>%
  summarise(
    average_R2 = mean(R2, na.rm = TRUE),
    n_pairs = dplyr::n(),
    .groups = "drop"
  ) %>%
  arrange(distance_interval)

message("Aggregation completed.")

write.csv(ld_data_summary, file.path(output_dir, "distance.tsv"), row.names = FALSE)
ld_data_zoom <- ld_data_summary %>% filter(distance_interval <= 100000)
write.csv(ld_data_zoom, file.path(output_dir, "distance_zoom.tsv"), row.names = FALSE)

# ----------------------------
# Summary stats
# ----------------------------
average_R2_after_100kb <- ld_data_summary %>%
  filter(distance_interval > 100000) %>%
  summarise(avg = mean(average_R2, na.rm = TRUE)) %>%
  pull(avg)

cleaned <- ld_data_summary %>%
  filter(!is.na(average_R2) & average_R2 > 0 & !is.na(distance_interval) & distance_interval > 0)

num_rows <- nrow(cleaned)
trend_n <- max(10, floor(num_rows / 10))
trend_data <- cleaned[1:trend_n, ]

log_model <- lm(average_R2 ~ log10(distance_interval), data = trend_data)
slope <- coef(log_model)[2]
intercept <- coef(log_model)[1]

recommended_LD_length <- NA
if (!is.na(average_R2_after_100kb) && is.finite(average_R2_after_100kb) &&
    is.finite(slope) && slope != 0 && is.finite(intercept)) {
  recommended_LD_length <- 10^((average_R2_after_100kb - intercept) / slope)
}

stats_file <- file.path(output_dir, "statistic.txt")
writeLines(c(
  paste("Average R2 after 100,000 bp:", average_R2_after_100kb),
  paste("Slope of log10 trend-line (fit on first 10% bins):", slope),
  paste("Y-intercept of log10 trend-line:", intercept),
  paste("Recommended LD length (heuristic):", recommended_LD_length),
  paste("Binning: <=100 kb:", BIN_SHORT, "bp; >100 kb:", BIN_LONG, "bp"),
  paste("Max distance filter used:", ifelse(is.na(opt$maxdist), "None (full range)", opt$maxdist)),
  paste("Circular distance:", ifelse(is.na(opt$genome_len), "No", paste0("Yes (L=", opt$genome_len, ")")))
), con = stats_file)

writeLines(paste("Recommended LD length (heuristic):", recommended_LD_length),
           con = file.path(output_dir, "recommended_ld_length.txt"))

message("Statistics saved.")

# ----------------------------
# esxQ reference line
# ----------------------------
esxQ_midpoint <- (3376490 + 3376852) / 2  # 3,376,671 bp

# ----------------------------
# X-axis formatting helpers
# ----------------------------
max_x <- max(cleaned$distance_interval, na.rm = TRUE)

# Sane tick spacing
step_linear <- if (max_x <= 200000) {
  20000
} else if (max_x <= 2000000) {
  200000
} else {
  500000
}
custom_breaks <- seq(0, max_x, by = step_linear)

# k / M label formatter — no scientific notation
fmt_bp <- function(x) {
  ifelse(x >= 1e6,
         paste0(formatC(x / 1e6, format = "f", digits = 1), "M"),
         ifelse(x >= 1e3,
                paste0(formatC(x / 1e3, format = "f", digits = 0), "k"),
                as.character(x)))
}

# Log-axis breaks: one per order of magnitude, labelled with fmt_bp
min_x_log      <- min(cleaned$distance_interval[cleaned$distance_interval > 0], na.rm = TRUE)
min_pow        <- floor(log10(min_x_log))
max_pow        <- ceiling(log10(max_x))
custom_log_breaks <- 10^(min_pow:max_pow)

# ----------------------------
# Plots
# ----------------------------
message("Starting plotting...")

# ---- Linear x ----
p <- ggplot(cleaned, aes(x = distance_interval, y = average_R2)) +
  geom_point(size = 1.2, alpha = 0.7) +
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE) +
  geom_vline(xintercept = esxQ_midpoint,
             linetype = "dashed", linewidth = 1, color = "red") +
  annotate("text",
           x = esxQ_midpoint,
           y = max(cleaned$average_R2, na.rm = TRUE),
           label = "esxQ", color = "red",
           angle = 90, vjust = -0.5, hjust = 1, size = 5) +
  scale_x_continuous(
    breaks = custom_breaks,
    labels = fmt_bp,
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(
    breaks = seq(0, 1.0, by = 0.1),
    limits = c(0, NA)
  ) +
  xlab("Distance (bp)") +
  ylab("Average R²") +
  ggtitle("LD Decay") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# ---- Log x ----
p2 <- ggplot(cleaned, aes(x = distance_interval, y = average_R2)) +
  geom_point(size = 1.2, alpha = 0.7) +
  geom_smooth(data = trend_data, method = "lm", formula = y ~ log10(x), se = FALSE) +
  geom_vline(xintercept = esxQ_midpoint,
             linetype = "dashed", linewidth = 1, color = "red") +
  annotate("text",
           x = esxQ_midpoint,
           y = max(cleaned$average_R2, na.rm = TRUE),
           label = "esxQ", color = "red",
           angle = 90, vjust = -0.5, hjust = 1, size = 5) +
  scale_x_log10(
    breaks = custom_log_breaks,
    labels = fmt_bp
  ) +
  scale_y_continuous(
    breaks = seq(0, 1.0, by = 0.1),
    limits = c(0, NA)
  ) +
  xlab("Distance (bp, log10 scale)") +
  ylab("Average R²") +
  ggtitle("LD Decay (log x)") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Save
ggsave(file.path(output_dir, "ld_decay.png"),     plot = p,  width = 12, height = 6, device = "png")
ggsave(file.path(output_dir, "ld_decay_log.png"), plot = p2, width = 12, height = 6, device = "png")
ggsave(file.path(output_dir, "ld_decay.svg"),     plot = p,  width = 12, height = 6, device = "svg")
ggsave(file.path(output_dir, "ld_decay_log.svg"), plot = p2, width = 12, height = 6, device = "svg")

message("Plots saved successfully:")
message(paste0(" - ", file.path(output_dir, "ld_decay.png")))
message(paste0(" - ", file.path(output_dir, "ld_decay_log.png")))
message(paste0(" - ", file.path(output_dir, "ld_decay.svg")))
message(paste0(" - ", file.path(output_dir, "ld_decay_log.svg")))
