#!/usr/bin/env Rscript
# =============================================================================
# Standard LD Decay Plot Around esxQ
# =============================================================================
# Input:  PLINK .ld file (pairwise LD from --r2)
# Output: LD decay plot (r² vs distance from focal gene)
#
# Usage:
#   Rscript esxQ_LD_decay_standard.R
#
# PLINK command to generate input:
#   plink --bfile <input> --r2 --ld-window-r2 0 \
#         --ld-window 999999 --ld-window-kb 1000 --out ld_results_sampled
# =============================================================================

library(readr)
library(dplyr)
library(ggplot2)

# -----------------------------------------------------------------------------
# 1. Load PLINK LD output
# -----------------------------------------------------------------------------
ld <- read_table("ld_results.ld", col_names = TRUE) %>%
  mutate(
    BP_A = as.numeric(BP_A),
    BP_B = as.numeric(BP_B),
    R2   = as.numeric(R2)
  ) %>%
  filter(!is.na(R2))

# -----------------------------------------------------------------------------
# 2. Define focal gene coordinates
# -----------------------------------------------------------------------------
esxq_start <- 3376490
esxq_end   <- 3376852
esxq_mid   <- (esxq_start + esxq_end) / 2

# -----------------------------------------------------------------------------
# 3. Filter: keep pairs where at least one SNP is within the gene region
#    Then compute distance from the OTHER SNP to esxQ midpoint
# -----------------------------------------------------------------------------
gene_window <- 500  # SNPs within this distance of esxQ midpoint are "in the gene"

ld_esxq <- ld %>%
  filter(
    abs(BP_A - esxq_mid) <= gene_window | abs(BP_B - esxq_mid) <= gene_window
  ) %>%
  mutate(
    # The "other" SNP position (the one NOT in esxQ)
    other_pos = ifelse(abs(BP_A - esxq_mid) <= gene_window, BP_B, BP_A),
    # Distance from esxQ midpoint (in bp)
    dist_bp   = abs(other_pos - esxq_mid),
    # Distance in kb for readability
    dist_kb   = dist_bp / 1000
  )

cat("Pairs involving esxQ SNPs:", nrow(ld_esxq), "\n")
cat("Max distance (kb):", max(ld_esxq$dist_kb), "\n")

# -----------------------------------------------------------------------------
# 4. LD Decay Plot — r² vs distance from esxQ
# -----------------------------------------------------------------------------
p_decay <- ggplot(ld_esxq, aes(x = dist_kb, y = R2)) +
  geom_point(size = 0.6, alpha = 0.3, color = "grey30") +
  geom_smooth(
    method = "loess",
    span   = 0.3,
    se     = TRUE,
    color  = "#d95f5f",
    fill   = "#d95f5f",
    alpha  = 0.15,
    linewidth = 1
  ) +
  geom_hline(yintercept = 0.2, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  annotate("text", x = max(ld_esxq$dist_kb) * 0.95, y = 0.22,
           label = "r² = 0.2", hjust = 1, size = 3, color = "grey40") +
  labs(
    x = "Distance from esxQ (kb)",
    y = expression(r^2),
    title = "LD Decay Around esxQ"
  ) +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.2)) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "italic", size = 13)
  )

print(p_decay)
ggsave("esxQ_LD_decay_full.png", p_decay, width = 7, height = 4, dpi = 300)

# -----------------------------------------------------------------------------
# 5. Regional LD Plot — r² vs absolute genomic position (centered on esxQ)
# -----------------------------------------------------------------------------
p_regional <- ggplot(ld_esxq, aes(x = other_pos / 1000, y = R2)) +
  geom_point(size = 0.6, alpha = 0.3, color = "grey30") +
  geom_smooth(
    method = "loess",
    span   = 0.3,
    se     = TRUE,
    color  = "#3575b5",
    fill   = "#3575b5",
    alpha  = 0.15,
    linewidth = 1
  ) +
  # Mark the gene region
  annotate("rect",
           xmin = esxq_start / 1000, xmax = esxq_end / 1000,
           ymin = -0.03, ymax = 1.05,
           fill = "#d95f5f", alpha = 0.1) +
  geom_vline(xintercept = esxq_mid / 1000, color = "#d95f5f",
             linetype = "dashed", linewidth = 0.6) +
  annotate("text", x = esxq_mid / 1000, y = 1.02,
           label = "esxQ", fontface = "italic", size = 4, color = "#d95f5f") +
  geom_hline(yintercept = 0.2, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  labs(
    x = "Genomic position (kb)",
    y = expression(r^2),
    title = "LD With esxQ Across the Region"
  ) +
  scale_y_continuous(limits = c(-0.03, 1.05), breaks = seq(0, 1, 0.2)) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "italic", size = 13)
  )

print(p_regional)
ggsave("esxQ_LD_regional_full.png", p_regional, width = 8, height = 4, dpi = 300)

# -----------------------------------------------------------------------------
# 6. Summary statistics
# -----------------------------------------------------------------------------
cat("\n--- LD Summary Around esxQ ---\n")
cat("Mean r²:", round(mean(ld_esxq$R2), 4), "\n")
cat("Median r²:", round(median(ld_esxq$R2), 4), "\n")

# Mean r² in distance bins
bins <- ld_esxq %>%
  mutate(bin_kb = cut(dist_kb,
                      breaks = c(0, 5, 10, 25, 50, 100, 250, 500, Inf),
                      right = FALSE,
                      labels = c("0-5", "5-10", "10-25", "25-50",
                                 "50-100", "100-250", "250-500", "500+"))) %>%
  group_by(bin_kb) %>%
  summarise(
    n       = n(),
    mean_r2 = round(mean(R2), 4),
    .groups = "drop"
  )

cat("\nMean r² by distance bin:\n")
print(as.data.frame(bins))
# =============================================================================
# Extract SNPs at the ~600-800 kb peak in LD with esxQ
# Run AFTER your main plot script so ld_esxq is already in memory
# =============================================================================
 
# Define the peak window (adjust based on your plot)
peak_min_kb <- 500
peak_max_kb <- 1000
r2_cutoff   <- 0.2   # only keep SNPs above this r² threshold
 
peak_snps <- ld_esxq %>%
  filter(dist_kb >= peak_min_kb & dist_kb <= peak_max_kb & R2 >= r2_cutoff) %>%
  arrange(desc(R2)) %>%
  select(SNP_A, BP_A, SNP_B, BP_B, R2, dist_kb, other_pos)
 
cat("SNPs in the", peak_min_kb, "-", peak_max_kb, "kb peak with r² >=", r2_cutoff, ":\n")
cat("Count:", nrow(peak_snps), "\n\n")
print(as.data.frame(peak_snps))
 
# Get unique SNP IDs at the peak (the "other" SNP, not the esxQ one)
peak_ids <- peak_snps %>%
  distinct(other_pos) %>%
  mutate(snp_id = case_when(
    other_pos %in% peak_snps$BP_A ~ peak_snps$SNP_A[match(other_pos, peak_snps$BP_A)],
    other_pos %in% peak_snps$BP_B ~ peak_snps$SNP_B[match(other_pos, peak_snps$BP_B)],
    TRUE ~ NA_character_
  ))
 
cat("\nUnique SNP positions at the peak:\n")
print(as.data.frame(peak_ids))
 
# Save to file
write.csv(peak_snps, "esxQ_peak_500_1000kb_snps.csv", row.names = FALSE)
cat("\nSaved to esxQ_peak_500_1000kb_snps.csv\n")
