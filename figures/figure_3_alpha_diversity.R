# =============================================================================
# Hospital Sink Disinfection Study - Figure 3: Alpha Diversity + PCA
#
# Authors: Megan Stroup Hill, Claude Sonnet 4
# Description: Generate Figure 3 showing alpha diversity and Aitchison PCA
# Input: alpha_diversity_data.RData, processed_data.RData
# Output: Figure 3 (alpha diversity boxplots + PCA per location)
# =============================================================================

source("analysis/00_load_libraries.R")

cat("Generating Figure 3: Alpha Diversity + PCA...\n")

# =============================================================================
# Load Data
# =============================================================================
if (!file.exists("analysis/alpha_diversity_data.RData")) {
  stop("Error: Run analysis/02_alpha_diversity.R first to generate required data")
}

# processed_data.RData holds meta_clean without diversity metrics, so load it
# first and let alpha_diversity_data.RData overwrite meta_clean.
load("analysis/processed_data.RData")
load("analysis/alpha_diversity_data.RData")

# =============================================================================
# Panel A: Alpha Diversity Boxplots with Significance
# =============================================================================
cat("Creating Panel A: Alpha diversity boxplots (rarefied to", RAREFY_DEPTH, "reads)...\n")

# Wilcoxon tests with BH correction, from analysis/02_alpha_diversity.R
alpha_significance <- alpha_results %>%
  select(Sample_Type, Metric, p_adjusted, significance_adj)

metric_labels <- c("Observed" = "Observed genera", "Shannon" = "Shannon")

# Prepare data for plotting (samples below the rarefaction depth are NA)
meta_long <- meta_clean %>%
  select(SeqID, Sample_Type, Treatment, Observed, Shannon) %>%
  pivot_longer(cols = c(Observed, Shannon), names_to = "Metric", values_to = "Value") %>%
  filter(!is.na(Value))

meta_long$Metric <- factor(metric_labels[meta_long$Metric], levels = metric_labels)
meta_long$Treatment <- factor(meta_long$Treatment, levels = c("Treated", "Untreated"))

# Calculate bracket positions
bracket_pos <- meta_long %>%
  group_by(Metric, Sample_Type) %>%
  summarise(ymax = max(Value, na.rm = TRUE), .groups = "drop")

# Merge significance with bracket positions
alpha_significance$Metric <- factor(metric_labels[alpha_significance$Metric], levels = metric_labels)
alpha_significance$Sample_Type <- factor(alpha_significance$Sample_Type, levels = SAMPLE_TYPE_LABELS)

bracket_data <- merge(alpha_significance, bracket_pos, by = c("Metric", "Sample_Type"))
bracket_data$y <- bracket_data$ymax * 1.08
bracket_data$ystar <- bracket_data$ymax * 1.12

panel_A <- ggplot(meta_long, aes(x = Sample_Type, y = Value, fill = Treatment)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3, linewidth = 0.4,
               position = position_dodge(width = 0.75), width = 0.65) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = TREATMENT_COLORS) +
  geom_segment(data = bracket_data,
               aes(x = as.numeric(Sample_Type) - 0.2, xend = as.numeric(Sample_Type) + 0.2,
                   y = y, yend = y), inherit.aes = FALSE, linewidth = 0.4) +
  geom_segment(data = bracket_data,
               aes(x = as.numeric(Sample_Type) - 0.2, xend = as.numeric(Sample_Type) - 0.2,
                   y = y, yend = y - (ymax * 0.02)), inherit.aes = FALSE, linewidth = 0.4) +
  geom_segment(data = bracket_data,
               aes(x = as.numeric(Sample_Type) + 0.2, xend = as.numeric(Sample_Type) + 0.2,
                   y = y, yend = y - (ymax * 0.02)), inherit.aes = FALSE, linewidth = 0.4) +
  geom_text(data = bracket_data,
            aes(x = as.numeric(Sample_Type), y = ystar, label = significance_adj),
            inherit.aes = FALSE, size = 4, vjust = 0) +
  labs(x = "Sample Location", y = NULL, fill = "Treatment", title = "A   Alpha Diversity") +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(size = 10),
    plot.title = element_text(face = "bold", size = 12),
    plot.margin = margin(5, 10, 5, 10)
  )

# =============================================================================
# Panel B: Aitchison PCA per Location
# =============================================================================
cat("Creating Panel B: Aitchison PCA per location...\n")

# Create CLR-transformed data if not available
if (!exists("comm_clr")) {
  pseudocount <- min(comm_rel[comm_rel > 0]) * 0.5
  comm_clr <- log(comm_rel + pseudocount)
  comm_clr <- comm_clr - rowMeans(comm_clr)
}

# Create PCA plots for each location
pca_plots <- list()
location_labels <- c("TP" = "Tail Piece", "PT" = "P-Trap", "TOP" = "Countertop")

for (st in SAMPLE_TYPE_LEVELS) {
  idx <- meta_clean$Sample_Type == SAMPLE_TYPE_LABELS[which(SAMPLE_TYPE_LEVELS == st)]
  sub_comm <- comm_rel[idx, ]
  sub_meta <- meta_clean[idx, ]

  # CLR transform
  sub_clr <- log(sub_comm + min(sub_comm[sub_comm > 0]) * 0.5)
  sub_clr <- sub_clr - rowMeans(sub_clr)

  # PCA on CLR-transformed data
  sub_pca <- prcomp(sub_clr, center = TRUE, scale. = FALSE)
  sub_eig <- round(100 * sub_pca$sdev[1:2]^2 / sum(sub_pca$sdev^2), 1)

  sub_meta$PC1 <- sub_pca$x[, 1]
  sub_meta$PC2 <- sub_pca$x[, 2]
  sub_meta$Treatment <- factor(sub_meta$Treatment, levels = c("Treated", "Untreated"))

  p <- ggplot(sub_meta, aes(x = PC1, y = PC2, color = Treatment)) +
    geom_point(alpha = 0.45, size = 1.3) +
    scale_color_manual(values = TREATMENT_COLORS) +
    stat_ellipse(level = 0.95, linewidth = 0.6) +
    labs(
      x = paste0("PC1 (", sub_eig[1], "%)"),
      y = paste0("PC2 (", sub_eig[2], "%)"),
      title = location_labels[st]
    ) +
    theme_bw(base_size = 10) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
      panel.grid.minor = element_blank(),
      plot.margin = margin(5, 5, 5, 5)
    )

  pca_plots[[st]] <- p
}

# Combine PCA plots in correct order
row_B <- pca_plots[["TP"]] + pca_plots[["PT"]] + pca_plots[["TOP"]] +
  plot_layout(ncol = 3)

# Add panel label
panel_B <- wrap_elements(row_B) +
  labs(title = "B   Aitchison PCA") +
  theme(plot.title = element_text(face = "bold", size = 12))

# =============================================================================
# Create Shared Legend
# =============================================================================
legend_plot <- ggplot(data.frame(
  Treatment = factor(c("Treated", "Untreated"), levels = c("Treated", "Untreated")),
  x = 1:2, y = 1:2),
  aes(x = x, y = y, fill = Treatment, color = Treatment)) +
  geom_point(size = 3, shape = 22) +
  scale_fill_manual(values = TREATMENT_COLORS) +
  scale_color_manual(values = TREATMENT_COLORS) +
  theme_void() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10))

shared_legend <- get_legend(legend_plot)

# =============================================================================
# Combine All Components
# =============================================================================
cat("Combining panels...\n")

combined_figure <- (panel_A / panel_B / wrap_elements(shared_legend)) +
  plot_layout(heights = c(4, 3.5, 0.5)) +
  plot_annotation(tag_levels = list(c("", "", "")))

# =============================================================================
# Save Figure
# =============================================================================
cat("Saving Figure 3...\n")

# Create figures directory if it doesn't exist
dir.create("figures/output", showWarnings = FALSE, recursive = TRUE)

# Save as PDF and PNG
ggsave("figures/output/Figure_3_alpha_diversity.pdf", combined_figure,
       width = 11, height = 9, device = "pdf")
ggsave("figures/output/Figure_3_alpha_diversity.png", combined_figure,
       width = 11, height = 9, dpi = 300)

cat("✓ Figure 3 saved to figures/output/\n")

# =============================================================================
# Figure Summary
# =============================================================================
cat("\n=== FIGURE 3 SUMMARY ===\n")
cat("Panel A: Observed genera and Shannon (rarefied to", RAREFY_DEPTH, "reads) with significance brackets\n")
cat("Panel B: Aitchison PCA per sample location with 95% confidence ellipses\n")
cat("Ordering: Tail Piece → P-Trap → Countertop\n")
cat("Statistical testing: Wilcoxon rank-sum with Benjamini-Hochberg correction\n")
cat("PCA: Performed on CLR-transformed data (appropriate for compositional data)\n")

# Print significance results
cat("\nSignificance results:\n")
print(alpha_significance[, c("Sample_Type", "Metric", "p_adjusted", "significance_adj")])

cat("\n✓ Figure 3 generation completed!\n")