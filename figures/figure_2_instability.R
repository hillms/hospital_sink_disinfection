# =============================================================================
# Hospital Sink Disinfection Study - Figure 2: Community Instability
#
# Authors: Megan Stroup Hill, Claude Sonnet 4
# Description: Generate Figure 2 showing community instability analysis
# Input: instability_analysis.RData
# Output: Figure 2 (ACF + PERMANOVA R² + joint scatter)
# =============================================================================

source("analysis/00_load_libraries.R")

cat("Generating Figure 2: Community Instability...\n")

# =============================================================================
# Load Data
# =============================================================================
if (!file.exists("analysis/instability_analysis.RData")) {
  stop("Error: Run analysis/06_community_instability.R first to generate required data")
}

load("analysis/instability_analysis.RData")

# =============================================================================
# Panel A: Lag-1 ACF Boxplots
# =============================================================================
cat("Creating Panel A: Lag-1 autocorrelation boxplots...\n")

# Calculate significance for each sample type
acf_significance <- acf_df %>%
  group_by(Sample_Type) %>%
  summarise(
    p_value = if(length(unique(Treatment)) == 2) {
      wilcox.test(lag1_acf ~ Treatment, exact = FALSE)$p.value
    } else {
      NA
    },
    significance = format_p_value(p_value),
    .groups = "drop"
  ) %>%
  mutate(
    y_pos = case_when(
      Sample_Type == "Tail Piece" ~ 0.72,
      Sample_Type == "P-Trap" ~ 0.57,
      Sample_Type == "Countertop" ~ 0.92
    )
  )

panel_A <- ggplot(acf_df,
                  aes(x = Treatment, y = lag1_acf, fill = Treatment, color = Treatment)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.45) +
  geom_boxplot(alpha = 0.30, outlier.shape = NA, linewidth = 0.55, width = 0.55) +
  geom_jitter(width = 0.14, size = 1.9, alpha = 0.75) +
  geom_text(data = acf_significance,
            aes(x = 1.5, y = y_pos, label = significance),
            inherit.aes = FALSE, size = 4.5, fontface = "bold") +
  facet_wrap(~ Sample_Type, ncol = 3) +
  scale_fill_manual(values = TREATMENT_COLORS) +
  scale_color_manual(values = TREATMENT_COLORS) +
  scale_y_continuous(limits = c(-0.65, 0.80)) +
  labs(
    x = NULL,
    y = "Lag-1 ACF",
    title = "A   Short-term predictability (Lag-1 ACF)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey92"),
    strip.text = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 11)
  )

# =============================================================================
# Panel B: PERMANOVA R² Boxplots
# =============================================================================
cat("Creating Panel B: PERMANOVA R² boxplots...\n")

# Calculate significance for PERMANOVA R²
r2_significance <- permanova_df %>%
  group_by(Sample_Type) %>%
  summarise(
    p_value = if(length(unique(Treatment)) == 2) {
      wilcox.test(R2_week ~ Treatment, exact = FALSE)$p.value
    } else {
      NA
    },
    significance = format_p_value(p_value),
    .groups = "drop"
  ) %>%
  mutate(
    y_pos = case_when(
      Sample_Type == "Tail Piece" ~ 0.185,
      Sample_Type == "P-Trap" ~ 0.225,
      Sample_Type == "Countertop" ~ 0.235
    )
  )

panel_B <- ggplot(permanova_df,
                  aes(x = Treatment, y = R2_week, fill = Treatment, color = Treatment)) +
  geom_boxplot(alpha = 0.30, outlier.shape = NA, linewidth = 0.55, width = 0.55) +
  geom_jitter(aes(shape = P_week < 0.05), width = 0.14, size = 1.9, alpha = 0.75) +
  scale_shape_manual(values = c("TRUE" = 19, "FALSE" = 1),
                     guide = guide_legend(title = "Week P < 0.05",
                                          override.aes = list(color = "grey30", size = 2))) +
  geom_text(data = r2_significance,
            aes(x = 1.5, y = y_pos, label = significance),
            inherit.aes = FALSE, size = 4.5, fontface = "bold") +
  facet_wrap(~ Sample_Type, ncol = 3) +
  scale_fill_manual(values = TREATMENT_COLORS) +
  scale_color_manual(values = TREATMENT_COLORS) +
  labs(
    x = NULL,
    y = expression(R^2 ~ "(Week)"),
    title = expression(bold("B   Directional drift (PERMANOVA ") * bolditalic(R^2) * bold(")"))
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "right",
    strip.background = element_rect(fill = "grey92"),
    strip.text = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 11)
  ) +
  guides(fill = "none", color = "none")

# =============================================================================
# Panel C: Joint Scatter Plot with Quadrants
# =============================================================================
cat("Creating Panel C: Joint ACF vs R² scatter...\n")

# Calculate centroids
centroids <- joint_data %>%
  group_by(Sample_Type, Treatment) %>%
  summarise(mean_R2 = mean(R2_week), mean_acf = mean(lag1_acf), .groups = "drop")

# Define quadrant background data
x_cut <- 0.085
y_cut <- 0.0

quad_df <- data.frame(
  Sample_Type = rep(factor(SAMPLE_TYPE_LABELS, levels = SAMPLE_TYPE_LABELS), each = 4),
  xmin = rep(c(-Inf, x_cut, -Inf, x_cut), 3),
  xmax = rep(c(x_cut, Inf, x_cut, Inf), 3),
  ymin = rep(c(y_cut, y_cut, -Inf, -Inf), 3),
  ymax = rep(c(Inf, Inf, y_cut, y_cut), 3),
  fill = rep(c("#d4eac8", "#c8dcea", "#f5e6cc", "#f5cccc"), 3),
  label = rep(c("Stable\n(high ACF, low R²)",
                "Predictable\ndrift",
                "Erratic\ninstability",
                "Directional\nbut noisy"), 3),
  lx = rep(c(-Inf + 0.001, x_cut + 0.001, -Inf + 0.001, x_cut + 0.001), 3),
  ly = rep(c(Inf, Inf, -Inf, -Inf), 3),
  hjust = rep(c(0, 0, 0, 0), 3),
  vjust = rep(c(1, 1, 0, 0), 3)
)

panel_C <- ggplot(joint_data, aes(x = R2_week, y = lag1_acf, color = Treatment, fill = Treatment)) +
  # Quadrant shading
  geom_rect(data = quad_df,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
            inherit.aes = FALSE, alpha = 0.18) +
  scale_fill_identity() +
  # Reference lines
  geom_vline(xintercept = x_cut, linetype = "dotted", color = "grey45", linewidth = 0.45) +
  geom_hline(yintercept = y_cut, linetype = "dotted", color = "grey45", linewidth = 0.45) +
  # Confidence ellipses
  stat_ellipse(aes(color = Treatment), level = 0.80, linewidth = 0.7, linetype = "solid", show.legend = FALSE) +
  # Individual room points with borders
  geom_point(size = 2.3, alpha = 0.80, shape = 21, aes(fill = Treatment), color = "white", stroke = 0.4) +
  # Group centroids
  geom_point(data = centroids, aes(x = mean_R2, y = mean_acf, fill = Treatment),
             shape = 23, size = 4, color = "black", stroke = 0.8, inherit.aes = FALSE) +
  # Quadrant labels
  geom_label(data = quad_df, aes(x = lx, y = ly, label = label, hjust = hjust, vjust = vjust),
             inherit.aes = FALSE, fill = NA, label.size = 0, color = "grey35",
             size = 2.4, lineheight = 0.9, label.padding = unit(0.1, "lines")) +
  facet_wrap(~ Sample_Type, ncol = 3) +
  scale_color_manual(values = TREATMENT_COLORS) +
  scale_fill_manual(values = TREATMENT_COLORS,
                    guide = guide_legend(override.aes = list(shape = 21, size = 3, color = "black", stroke = 0.5))) +
  scale_x_continuous(labels = number_format(accuracy = 0.01)) +
  labs(
    x = expression(R^2 ~ "(Week) — directional drift"),
    y = "Lag-1 ACF — short-term predictability",
    title = "C   Joint instability landscape per room",
    fill = "Treatment", color = "Treatment"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey92"),
    strip.text = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 11)
  )

# =============================================================================
# Combine All Panels
# =============================================================================
cat("Combining panels...\n")

combined_figure <- gridExtra::arrangeGrob(
  panel_A, panel_B, panel_C,
  ncol = 1,
  heights = c(1, 1, 1.2)
)

# =============================================================================
# Save Figure
# =============================================================================
cat("Saving Figure 2...\n")

# Create figures directory if it doesn't exist
dir.create("figures/output", showWarnings = FALSE, recursive = TRUE)

# Save as PDF and PNG
ggsave("figures/output/Figure_2_instability.pdf", combined_figure,
       width = 10, height = 13, device = "pdf")
ggsave("figures/output/Figure_2_instability.png", combined_figure,
       width = 10, height = 13, dpi = 300)

cat("✓ Figure 2 saved to figures/output/\n")

# =============================================================================
# Figure Summary
# =============================================================================
cat("\n=== FIGURE 2 SUMMARY ===\n")
cat("Panel A: Lag-1 autocorrelation boxplots (short-term predictability)\n")
cat("Panel B: PERMANOVA R² boxplots (directional drift over time)\n")
cat("Panel C: Joint scatter plot showing instability landscape\n")
cat("Key finding: Treatment creates 'controlled instability' signature\n")
cat("Quadrants: Stable | Predictable drift | Erratic | Directional+noisy\n")
cat("Treatment signature: Low ACF + Low R² = Controlled instability\n")

cat("\n✓ Figure 2 generation completed!\n")