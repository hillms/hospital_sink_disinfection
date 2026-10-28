# =============================================================================
# Hospital Sink Disinfection Study - Figure 5: Differential Abundance
#
# Authors: Megan Stroup Hill, Claude
# Description: Volcano plots of ALDEx2 results (raw genus counts) for the
#              tail piece and P-trap
# Input: analysis/aldex2_results_raw_counts.csv (04_differential_abundance.R)
# Output: Figure 5 (volcano plots with site labels)
#
# Labels: ESKAPEE genera (bold, enlarged points), genera discussed in the text,
# and the 5 largest significant effects on each side. The key to the plot
# (effect-size sign, bold labels, dashed lines) is in the figure legend.
# =============================================================================

source("analysis/00_load_libraries.R")

cat("Generating Figure 5: Differential Abundance...\n")

if (!file.exists("analysis/aldex2_results_raw_counts.csv")) {
  stop("Error: Run analysis/04_differential_abundance.R first to generate required data")
}
results <- read.csv("analysis/aldex2_results_raw_counts.csv")

# =============================================================================
# Volcano Plot Function
# =============================================================================
create_volcano_plot <- function(data, site_label) {
  # Determine axis limits for consistency
  max_effect <- max(abs(data$effect), na.rm = TRUE) * 1.35
  max_log_p <- max(-log10(pmin(data$we.eBH, data$wi.eBH)), na.rm = TRUE) * 1.1

  # Create significance and color categories
  data$plot_category <- "Not significant"
  data$plot_category[data$significant & data$effect < 0] <- "Treated-enriched"
  data$plot_category[data$significant & data$effect > 0] <- "Untreated-enriched"

  # Color palette
  colors <- c("Not significant" = "grey70",
              "Treated-enriched" = "#457B9D",
              "Untreated-enriched" = "#D62828")

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
      title = site_label  # Add site label as title
    ) +

    # Theme
    theme_bw(base_size = 11) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = 9),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12)  # Center and bold title
    )

  # Add text labels for significant genera
  # Enhanced labeling: prioritize ESKAPEE pathogens and more treated-enriched taxa
  # Label ESKAPEE genera, genera discussed in the text and the 5 largest effects per side
  text_g <- c("Afipia", "Bradyrhizobium", "Pelomonas", "Sphingomonas", "Comamonas")
  sig <- data %>% filter(significant)
  top <- bind_rows(sig %>% filter(effect < 0) %>% arrange(effect) %>% slice_head(n = 5),
                   sig %>% filter(effect > 0) %>% arrange(desc(effect)) %>% slice_head(n = 5))
  label_data <- sig %>%
    filter(is_eskapee | genus %in% text_g | genus %in% top$genus)

  if (nrow(label_data) > 0) {
    # Split labels by side for better positioning
    left_labels <- label_data %>% filter(effect < 0)  # Treated-enriched (left side)
    right_labels <- label_data %>% filter(effect > 0)  # Untreated-enriched (right side)

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
          max.overlaps = 30,
          xlim = c(-Inf, -0.02),  # Keep labels on left side
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
          max.overlaps = 30,
          xlim = c(0.02, Inf),  # Keep labels on right side
          seed = 42
        )
    }
  }

  return(p)
}

# =============================================================================
# Create and Combine Plots (Tail Piece first)
# =============================================================================
tp_plot <- create_volcano_plot(results %>% filter(Sample_Type == "TP"), "Tail Piece")
pt_plot <- create_volcano_plot(results %>% filter(Sample_Type == "PT"), "P-Trap")

final_plot <- tp_plot + pt_plot +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

# =============================================================================
# Save Figure
# =============================================================================
dir.create("figures/output", showWarnings = FALSE, recursive = TRUE)
ggsave("figures/output/Figure_5_differential_abundance.pdf", final_plot,
       width = 12, height = 8, device = "pdf")
ggsave("figures/output/Figure_5_differential_abundance.png", final_plot,
       width = 12, height = 8, dpi = 300)

cat("✓ Figure 5 saved to figures/output/\n")
print(results %>%
        filter(Sample_Type != "TOP") %>%
        group_by(Sample_Type) %>%
        summarise(tested = n(),
                  treated_enriched = sum(significant & effect < 0),
                  untreated_enriched = sum(significant & effect > 0),
                  significant = sum(significant)))
