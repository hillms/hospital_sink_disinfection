# =============================================================================
# Hospital Sink Disinfection Study - Figure 5: Differential Abundance
#
# Authors: Megan Stroup Hill, Claude Sonnet 4
# Description: Generate Figure 5 showing differential abundance volcano plots
# Input: processed_data.RData
# Output: Figure 5 (volcano plots with site labels)
# =============================================================================

source("analysis/00_load_libraries.R")

cat("Generating Figure 5: Differential Abundance...\n")

# =============================================================================
# Load Data
# =============================================================================
load("analysis/processed_data.RData")

# Check if ALDEx2 is available
if (!requireNamespace("ALDEx2", quietly = TRUE)) {
  stop("ALDEx2 package required. Install with: BiocManager::install('ALDEx2')")
}

library(ALDEx2)

# =============================================================================
# Prepare Data for ALDEx2
# =============================================================================
cat("Preparing data for ALDEx2 analysis...\n")

# ALDEx2 requires integer counts - convert from relative abundances
# We'll need to reconstruct approximate counts (this is an approximation)
mean_lib_size <- 10000  # Assume a standard library size for conversion
genus_counts <- round(comm_rel * mean_lib_size)

# Remove zero-sum samples
keep_samples <- rowSums(genus_counts) > 0
genus_counts <- genus_counts[keep_samples, ]
meta_analysis <- meta_clean[keep_samples, ]

cat("✓ Data prepared:", nrow(genus_counts), "samples ×", ncol(genus_counts), "genera\n")

# =============================================================================
# Run ALDEx2 for Each Location
# =============================================================================
cat("Running ALDEx2 differential abundance analysis...\n")

aldex_results <- list()

# Analyze Tail Piece and P-Trap (where we expect the strongest effects)
for (st in c("TP", "PT")) {
  st_label <- SAMPLE_TYPE_LABELS[which(SAMPLE_TYPE_LEVELS == st)]
  cat("Analyzing", st_label, "samples...\n")

  # Subset data for this sample type
  sub_meta <- meta_analysis %>% filter(Sample_Type == st_label)
  sub_counts <- t(genus_counts[sub_meta$SeqID, ])

  # Filter low-abundance genera (>=5% prevalence)
  prevalence <- rowSums(sub_counts > 0) / ncol(sub_counts)
  sub_counts <- sub_counts[prevalence >= 0.05, ]

  cat("  Testing", nrow(sub_counts), "genera (>=5% prevalence)\n")

  if (nrow(sub_counts) >= 10 && length(unique(sub_meta$Treatment)) == 2) {
    # Create conditions vector
    conditions <- sub_meta$Treatment

    # Run ALDEx2
    set.seed(42)
    aldex_out <- aldex(sub_counts, conditions, mc.samples = 128, test = "t",
                       effect = TRUE, verbose = FALSE)

    # Add metadata
    aldex_out$genus <- rownames(aldex_out)
    aldex_out$location <- st
    aldex_out$location_label <- st_label

    # Identify ESKAPEE pathogens
    aldex_out$is_eskapee <- aldex_out$genus %in% ESKAPEE_GENERA

    # Significance categories
    aldex_out$significant <- (aldex_out$wi.eBH < 0.05 | aldex_out$we.eBH < 0.05)
    aldex_out$direction <- ifelse(aldex_out$effect < 0, "Treated-enriched", "Untreated-enriched")

    aldex_results[[st]] <- aldex_out

    cat("  ✓ Significant genera (BH < 0.05):", sum(aldex_out$significant), "\n")
  } else {
    cat("  ⚠ Insufficient data for", st_label, "\n")
  }
}

if (length(aldex_results) == 0) {
  stop("No successful ALDEx2 analyses. Check data quality.")
}

# =============================================================================
# Create Volcano Plot Function
# =============================================================================
create_volcano_plot <- function(data, site_label) {
  # Determine axis limits for consistency
  max_effect <- max(abs(data$effect), na.rm = TRUE) * 1.1
  max_log_p <- max(-log10(pmin(data$we.eBH, data$wi.eBH)), na.rm = TRUE) * 1.1

  # Create significance and color categories
  data$plot_category <- "Not significant"
  data$plot_category[data$significant & data$effect < 0] <- "Treated-enriched"
  data$plot_category[data$significant & data$effect > 0] <- "Untreated-enriched"

  # Color palette
  colors <- c("Not significant" = "grey70",
              "Treated-enriched" = TREATMENT_COLORS["Treated"],
              "Untreated-enriched" = TREATMENT_COLORS["Untreated"])

  # Base plot with site label
  p <- ggplot(data, aes(x = effect, y = -log10(pmin(we.eBH, wi.eBH)))) +
    geom_point(aes(color = plot_category, size = is_eskapee, alpha = significant)) +
    scale_color_manual(values = colors, name = "Differential abundance") +
    scale_size_manual(values = c("FALSE" = 1.5, "TRUE" = 2.5), guide = "none") +
    scale_alpha_manual(values = c("FALSE" = 0.5, "TRUE" = 0.9), guide = "none") +

    # Reference lines
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", alpha = 0.7) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", alpha = 0.7) +

    # Axes and labels with site identification
    xlim(-max_effect, max_effect) +
    ylim(0, max_log_p) +
    labs(
      x = "Effect size (ALDEx2)",
      y = expression(-log[10]("BH-adjusted P-value")),
      title = site_label
    ) +

    # Theme
    theme_bw(base_size = 11) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = 9),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12)
    )

  # Add text labels for significant genera
  label_data <- data %>%
    filter(significant) %>%
    mutate(
      priority = case_when(
        is_eskapee ~ 4,  # Highest priority - ESKAPEE pathogens
        effect < -0.1 ~ 3,  # High priority - treated-enriched
        abs(effect) > 0.3 ~ 2,  # High effect size untreated-enriched
        TRUE ~ 1  # Other significant
      )
    ) %>%
    arrange(desc(priority), desc(abs(effect))) %>%
    slice_head(n = 15)

  if (nrow(label_data) > 0) {
    # Split labels by side for better positioning
    left_labels <- label_data %>% filter(effect < 0)  # Treated-enriched
    right_labels <- label_data %>% filter(effect > 0)  # Untreated-enriched

    # Add labels for treated-enriched (left side)
    if (nrow(left_labels) > 0) {
      p <- p +
        geom_text_repel(
          data = left_labels,
          aes(label = genus),
          size = 3,
          fontface = ifelse(left_labels$is_eskapee, "bold", "italic"),
          color = "black",
          bg.color = "white",
          bg.r = 0.1,
          box.padding = 0.3,
          point.padding = 0.3,
          max.overlaps = 10,
          xlim = c(-Inf, -0.02),
          seed = 42
        )
    }

    # Add labels for untreated-enriched (right side)
    if (nrow(right_labels) > 0) {
      p <- p +
        geom_text_repel(
          data = right_labels,
          aes(label = genus),
          size = 3,
          fontface = ifelse(right_labels$is_eskapee, "bold", "italic"),
          color = "black",
          bg.color = "white",
          bg.r = 0.1,
          box.padding = 0.3,
          point.padding = 0.3,
          max.overlaps = 10,
          xlim = c(0.02, Inf),
          seed = 42
        )
    }
  }

  return(p)
}

# =============================================================================
# Create Individual Volcano Plots
# =============================================================================
cat("Creating volcano plots...\n")

individual_plots <- list()

# Create plots in Tail Piece → P-Trap order
for (st in c("TP", "PT")) {
  if (st %in% names(aldex_results)) {
    st_label <- SAMPLE_TYPE_LABELS[which(SAMPLE_TYPE_LEVELS == st)]
    individual_plots[[st]] <- create_volcano_plot(aldex_results[[st]], st_label)
  }
}

# =============================================================================
# Combine Plots
# =============================================================================
cat("Combining volcano plots...\n")

if (length(individual_plots) == 2) {
  # Both plots available - combine in Tail Piece → P-Trap order
  combined_plot <- individual_plots[["TP"]] + individual_plots[["PT"]] +
    plot_layout(ncol = 2, guides = "collect") &
    theme(legend.position = "bottom")
} else if (length(individual_plots) == 1) {
  # Only one plot available
  combined_plot <- individual_plots[[1]]
} else {
  stop("No successful volcano plots generated")
}

# Add interpretation subtitle
final_plot <- combined_plot +
  plot_annotation(
    subtitle = "Effect size > 0: enriched in untreated samples; Effect size < 0: enriched in treated samples\nBold labels: ESKAPEE pathogens; Dashed lines: significance thresholds",
    theme = theme(
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey30")
    )
  )

# =============================================================================
# Save Figure
# =============================================================================
cat("Saving Figure 5...\n")

# Create figures directory if it doesn't exist
dir.create("figures/output", showWarnings = FALSE, recursive = TRUE)

# Save as PDF and PNG
ggsave("figures/output/Figure_5_differential_abundance.pdf", final_plot,
       width = 12, height = 8, device = "pdf")
ggsave("figures/output/Figure_5_differential_abundance.png", final_plot,
       width = 12, height = 8, dpi = 300)

cat("✓ Figure 5 saved to figures/output/\n")

# =============================================================================
# Generate Summary Statistics
# =============================================================================
cat("Generating differential abundance summary...\n")

summary_results <- data.frame()

for (st in names(aldex_results)) {
  data <- aldex_results[[st]]
  st_label <- SAMPLE_TYPE_LABELS[which(SAMPLE_TYPE_LEVELS == st)]

  # Overall statistics
  total_tested <- nrow(data)
  significant <- sum(data$significant)
  treated_enriched <- sum(data$significant & data$effect < 0)
  untreated_enriched <- sum(data$significant & data$effect > 0)
  eskapee_significant <- sum(data$significant & data$is_eskapee)

  summary_results <- rbind(summary_results, data.frame(
    Location = st_label,
    Total_Tested = total_tested,
    Significant = significant,
    Treated_Enriched = treated_enriched,
    Untreated_Enriched = untreated_enriched,
    ESKAPEE_Significant = eskapee_significant
  ))

  # Top enriched in each direction
  cat("\n", st_label, "results:\n")
  cat("  Total genera tested:", total_tested, "\n")
  cat("  Significant genera:", significant, "\n")
  cat("  Treated-enriched:", treated_enriched, "\n")
  cat("  Untreated-enriched:", untreated_enriched, "\n")
  cat("  ESKAPEE pathogens significant:", eskapee_significant, "\n")

  if (treated_enriched > 0) {
    top_treated <- data %>%
      filter(significant, effect < 0) %>%
      arrange(effect) %>%
      slice_head(n = 3)

    cat("  Top treated-enriched:\n")
    for (i in 1:nrow(top_treated)) {
      cat(sprintf("    %s (effect = %.3f)\n",
                  top_treated$genus[i], top_treated$effect[i]))
    }
  }

  if (untreated_enriched > 0) {
    top_untreated <- data %>%
      filter(significant, effect > 0) %>%
      arrange(desc(effect)) %>%
      slice_head(n = 3)

    cat("  Top untreated-enriched:\n")
    for (i in 1:nrow(top_untreated)) {
      cat(sprintf("    %s (effect = %.3f)\n",
                  top_untreated$genus[i], top_untreated$effect[i]))
    }
  }
}

# Save summary
write.csv(summary_results, "analysis/differential_abundance_summary.csv", row.names = FALSE)

# =============================================================================
# Figure Summary
# =============================================================================
cat("\n=== FIGURE 5 SUMMARY ===\n")
cat("Dual volcano plots showing differential abundance analysis\n")
cat("Method: ALDEx2 with 128 Monte Carlo samples\n")
cat("Comparison: Treated vs Untreated for each sample location\n")
cat("Significance: Benjamini-Hochberg corrected p-values < 0.05\n")
cat("Effect size: ALDEx2 effect size (log2-like scale)\n")
cat("Site labels: Clear identification of Tail Piece and P-Trap\n")

cat("\nKey findings:\n")
print(summary_results)

cat("\nInterpretation:\n")
cat("- Left side (effect < 0): Genera enriched in treated samples\n")
cat("- Right side (effect > 0): Genera enriched in untreated samples\n")
cat("- Bold labels: ESKAPEE pathogens (clinical concern)\n")
cat("- Point size: Larger for ESKAPEE pathogens\n")

cat("\n✓ Figure 5 generation completed!\n")