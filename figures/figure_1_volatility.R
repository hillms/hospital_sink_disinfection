# =============================================================================
# Hospital Sink Disinfection Study - Figure 1: Community Volatility
#
# Authors: Megan Stroup Hill, Claude Sonnet 4
# Description: Generate Figure 1 showing community volatility analysis
# Input: volatility_analysis.RData
# Output: Figure 1 (boxplot + trajectory panels)
# =============================================================================

source("analysis/00_load_libraries.R")

cat("Generating Figure 1: Community Volatility...\n")

# =============================================================================
# Load Data
# =============================================================================
if (!file.exists("analysis/volatility_analysis.RData")) {
  stop("Error: Run analysis/05_community_volatility.R first to generate required data")
}

load("analysis/volatility_analysis.RData")

# =============================================================================
# Panel A: Volatility Boxplots
# =============================================================================
cat("Creating Panel A: Volatility boxplots...\n")

panel_A <- ggplot(volatility_df,
                  aes(x = Treatment, y = mean_volatility,
                      fill = Treatment, color = Treatment)) +
  geom_boxplot(alpha = 0.35, outlier.shape = NA, linewidth = 0.6) +
  geom_jitter(width = 0.18, size = 1.8, alpha = 0.7) +
  facet_wrap(~ Sample_Type, ncol = 3) +
  scale_fill_manual(values = TREATMENT_COLORS) +
  scale_color_manual(values = TREATMENT_COLORS) +
  labs(
    x = NULL,
    y = "Mean Aitchison distance\n(consecutive weeks)",
    title = "A   Community Volatility"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey92"),
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 10),
    plot.title = element_text(face = "bold", size = 12),
    plot.margin = margin(20, 20, 20, 40, "pt")
  )

# =============================================================================
# Panel B: Temporal Trajectories
# =============================================================================
cat("Creating Panel B: Temporal trajectories...\n")

# Calculate weekly averages and standard errors
weekly_trajectory <- step_data %>%
  group_by(Sample_Type, Treatment, Week_to) %>%
  summarise(
    mean_distance = mean(distance, na.rm = TRUE),
    se_distance = sd(distance, na.rm = TRUE) / sqrt(n()),
    n_obs = n(),
    .groups = "drop"
  ) %>%
  filter(n_obs >= 3)  # Require at least 3 observations

panel_B <- ggplot(weekly_trajectory,
                  aes(x = Week_to, y = mean_distance,
                      color = Treatment, fill = Treatment)) +
  geom_ribbon(aes(ymin = pmax(mean_distance - se_distance, 0),
                  ymax = mean_distance + se_distance),
              alpha = 0.18, color = NA) +
  geom_line(linewidth = 0.75) +
  geom_point(size = 1.4) +
  facet_wrap(~ Sample_Type, ncol = 3) +
  scale_color_manual(values = TREATMENT_COLORS) +
  scale_fill_manual(values = TREATMENT_COLORS) +
  scale_x_continuous(breaks = seq(0, 35, by = 5)) +
  labs(
    x = "Week",
    y = "Mean Aitchison distance (± SE)",
    title = "B   Temporal Dynamics",
    color = "Treatment",
    fill = "Treatment"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey92"),
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 12),
    plot.margin = margin(20, 20, 40, 40, "pt"),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10)
  )

# =============================================================================
# Combine Panels
# =============================================================================
cat("Combining panels...\n")

combined_figure <- panel_A / panel_B +
  plot_layout(heights = c(1, 1.2))

# =============================================================================
# Save Figure
# =============================================================================
cat("Saving Figure 1...\n")

# Create figures directory if it doesn't exist
dir.create("figures/output", showWarnings = FALSE, recursive = TRUE)

# Save as PDF and PNG
ggsave("figures/output/Figure_1_volatility.pdf", combined_figure,
       width = 10, height = 8, device = "pdf")
ggsave("figures/output/Figure_1_volatility.png", combined_figure,
       width = 10, height = 8, dpi = 300)

cat("✓ Figure 1 saved to figures/output/\n")

# =============================================================================
# Generate Figure Summary
# =============================================================================
cat("Generating figure summary...\n")

# Load statistics if available
if (exists("volatility_stats")) {
  figure_stats <- volatility_stats %>%
    mutate(
      interpretation = case_when(
        significance %in% c("****", "***", "**", "*") ~
          paste("Significant", ifelse(treated_mean > untreated_mean, "increase", "decrease")),
        TRUE ~ "No significant difference"
      )
    ) %>%
    select(Sample_Type, interpretation, p_value)

  print("Statistical significance by sample type:")
  print(figure_stats)
}

# Weekly significance summary if available
if (exists("weekly_summary")) {
  print("Weekly volatility patterns:")
  print(weekly_summary)
}

cat("\n=== FIGURE 1 SUMMARY ===\n")
cat("Panel A: Community volatility boxplots by treatment and sample type\n")
cat("Panel B: Temporal trajectories of weekly volatility over 35 weeks\n")
cat("Color scheme: Treated (blue), Untreated (red)\n")
cat("Ordering: Tail Piece → P-Trap → Countertop\n")
cat("Metric: Aitchison distance (Euclidean distance on CLR-transformed data)\n")
cat("Interpretation: Week-to-week compositional turnover\n")

cat("\n✓ Figure 1 generation completed!\n")