# =============================================================================
# Hospital Sink Disinfection Study - Figure 4: Top Taxa Stacked Bar Chart
#
# Authors: Megan Stroup Hill, Claude Sonnet 4
# Description: Generate Figure 4 showing top 15 genera composition
# Input: processed_data.RData
# Output: Figure 4 (stacked bar chart by treatment and location)
# =============================================================================

source("analysis/00_load_libraries.R")

cat("Generating Figure 4: Top Taxa Composition...\n")

# =============================================================================
# Load Data
# =============================================================================
load("analysis/processed_data.RData")

# =============================================================================
# Calculate Mean Relative Abundances by Group
# =============================================================================
cat("Calculating mean relative abundances by group...\n")

group_data <- data.frame()

# Calculate for each location and treatment combination
for (loc in SAMPLE_TYPE_LEVELS) {
  for (trt in c("Treated", "Untreated")) {
    subset_idx <- meta_clean$Sample_Type == SAMPLE_TYPE_LABELS[which(SAMPLE_TYPE_LEVELS == loc)] &
                  meta_clean$Treatment == trt

    if (sum(subset_idx) > 0) {
      subset_comm <- comm_rel[subset_idx, , drop = FALSE]
      mean_abundance <- colMeans(subset_comm)

      group_data <- rbind(group_data, data.frame(
        Location = loc,
        Treatment = trt,
        Genus = names(mean_abundance),
        Relative_Abundance = mean_abundance,
        stringsAsFactors = FALSE
      ))
    }
  }
}

# Create location labels
location_labels <- c("TP" = "Tail Piece", "PT" = "P-Trap", "TOP" = "Countertop")
group_data$Location_Label <- location_labels[group_data$Location]

# Create combined group labels in correct order
group_data$Group <- paste(group_data$Location_Label, group_data$Treatment, sep = "_")

cat("✓ Group abundances calculated for", length(unique(group_data$Group)), "groups\n")

# =============================================================================
# Identify Top Taxa and Prepare for Plotting
# =============================================================================
cat("Identifying top taxa for plotting...\n")

# Get top 15 taxa by overall abundance
top_taxa <- group_data %>%
  group_by(Genus) %>%
  summarise(Total_Abundance = sum(Relative_Abundance), .groups = "drop") %>%
  arrange(desc(Total_Abundance)) %>%
  slice_head(n = 15) %>%
  pull(Genus)

cat("✓ Top 15 genera identified:", paste(head(top_taxa, 5), collapse = ", "), "...\n")

# Prepare data for stacking with correct ordering
plot_data <- group_data %>%
  mutate(
    Genus_Cat = ifelse(Genus %in% top_taxa, Genus, "Other"),
    # Correct factor ordering: Tail Piece → P-Trap → Countertop
    Group = factor(Group, levels = c("Tail Piece_Treated", "Tail Piece_Untreated",
                                   "P-Trap_Treated", "P-Trap_Untreated",
                                   "Countertop_Treated", "Countertop_Untreated"))
  ) %>%
  group_by(Group, Genus_Cat) %>%
  summarise(Relative_Abundance = sum(Relative_Abundance), .groups = "drop")

# Order taxa by abundance for consistent coloring
taxa_order <- plot_data %>%
  group_by(Genus_Cat) %>%
  summarise(Total = sum(Relative_Abundance), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  pull(Genus_Cat)

plot_data$Genus_Cat <- factor(plot_data$Genus_Cat, levels = taxa_order)

cat("✓ Plot data prepared with", length(levels(plot_data$Genus_Cat)), "taxa categories\n")

# =============================================================================
# Create Color Palette
# =============================================================================
cat("Creating color palette...\n")

n_taxa <- length(unique(plot_data$Genus_Cat))

# Create a diverse color palette
set.seed(42)  # For reproducible colors

# Use a combination of color palettes for better distinction
colors1 <- brewer.pal(min(11, n_taxa), "Spectral")
if (n_taxa > 11) {
  colors2 <- brewer.pal(min(8, n_taxa - 11), "Set2")
  if (n_taxa > 19) {
    colors3 <- viridis(n_taxa - 19, option = "plasma")
    color_palette <- c(colors1, colors2, colors3)
  } else {
    color_palette <- c(colors1, colors2)
  }
} else {
  color_palette <- colors1
}

# Make sure "Other" is grey if present
if ("Other" %in% levels(plot_data$Genus_Cat)) {
  other_idx <- which(levels(plot_data$Genus_Cat) == "Other")
  color_palette[other_idx] <- "grey80"
}

names(color_palette) <- levels(plot_data$Genus_Cat)

# =============================================================================
# Create the Stacked Bar Plot
# =============================================================================
cat("Creating stacked bar plot...\n")

stacked_plot <- ggplot(plot_data, aes(x = Group, y = Relative_Abundance, fill = Genus_Cat)) +
  geom_bar(stat = "identity", color = "white", size = 0.1) +
  scale_fill_manual(values = color_palette, name = "Genus") +
  scale_y_continuous(labels = scales::percent_format(scale = 100),
                     expand = c(0, 0)) +
  labs(
    x = NULL,
    y = "Relative Abundance (%)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    # Clean theme
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line.y = element_line(color = "black", size = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 11, margin = margin(r = 10)),
    legend.position = "right",
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10, face = "bold"),
    legend.key.size = unit(0.4, "cm"),
    legend.margin = margin(l = 10),
    plot.margin = margin(10, 10, 10, 10),
    strip.background = element_blank(),
    strip.text = element_text(size = 10, face = "bold")
  ) +
  guides(fill = guide_legend(ncol = 1, reverse = TRUE))

# =============================================================================
# Save Figure
# =============================================================================
cat("Saving Figure 4...\n")

# Create figures directory if it doesn't exist
dir.create("figures/output", showWarnings = FALSE, recursive = TRUE)

# Save as PDF and PNG
ggsave("figures/output/Figure_4_top_taxa.pdf", stacked_plot,
       width = 10, height = 6, device = "pdf")
ggsave("figures/output/Figure_4_top_taxa.png", stacked_plot,
       width = 10, height = 6, dpi = 300)

cat("✓ Figure 4 saved to figures/output/\n")

# =============================================================================
# Generate Summary Statistics
# =============================================================================
cat("Generating summary statistics...\n")

# Calculate summary by group
summary_stats <- plot_data %>%
  group_by(Group) %>%
  summarise(
    N_Taxa = n(),
    Total_Abundance = sum(Relative_Abundance),
    .groups = "drop"
  )

cat("Data summary by group:\n")
print(summary_stats)

# Top taxa information
top_taxa_info <- plot_data %>%
  filter(Genus_Cat %in% top_taxa) %>%
  group_by(Genus_Cat) %>%
  summarise(
    Total_Abundance = sum(Relative_Abundance),
    Max_Group = Group[which.max(Relative_Abundance)],
    Max_Abundance = max(Relative_Abundance),
    .groups = "drop"
  ) %>%
  arrange(desc(Total_Abundance))

cat("\nTop 10 most abundant genera:\n")
print(head(top_taxa_info, 10))

# =============================================================================
# Figure Summary
# =============================================================================
cat("\n=== FIGURE 4 SUMMARY ===\n")
cat("Stacked bar chart showing relative abundance of top 15 genera\n")
cat("Grouping: Tail Piece → P-Trap → Countertop, each with Treated/Untreated\n")
cat("Color coding: Distinct colors for each genus, grey for 'Other'\n")
cat("Data: Mean relative abundances across samples within each group\n")
cat("Total taxa represented:", n_taxa, "(", length(top_taxa), "named +",
    ifelse("Other" %in% levels(plot_data$Genus_Cat), "Other", "no Other"), ")\n")

# Treatment-specific patterns
treatment_patterns <- plot_data %>%
  separate(Group, into = c("Location", "Treatment"), sep = "_", remove = FALSE) %>%
  group_by(Location, Genus_Cat) %>%
  summarise(
    Treated_Abundance = sum(Relative_Abundance[Treatment == "Treated"]),
    Untreated_Abundance = sum(Relative_Abundance[Treatment == "Untreated"]),
    .groups = "drop"
  ) %>%
  mutate(
    Treatment_Ratio = Treated_Abundance / (Untreated_Abundance + 1e-6),  # Avoid division by zero
    Enrichment = case_when(
      Treatment_Ratio > 1.5 ~ "Treated-enriched",
      Treatment_Ratio < 0.67 ~ "Untreated-enriched",
      TRUE ~ "Similar"
    )
  ) %>%
  filter(Genus_Cat != "Other")

cat("\nTreatment enrichment patterns (top genera):\n")
enrichment_summary <- treatment_patterns %>%
  count(Location, Enrichment) %>%
  pivot_wider(names_from = Enrichment, values_from = n, values_fill = 0)

print(enrichment_summary)

cat("\n✓ Figure 4 generation completed!\n")