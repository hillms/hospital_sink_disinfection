# =============================================================================
# Hospital Sink Disinfection Study - Alpha Diversity Analysis
#
# Authors: Megan Stroup Hill, Claude Sonnet 4
# Description: Calculate and analyze alpha diversity metrics
# Input: processed_data.RData
# Output: Alpha diversity metrics and statistical comparisons
# =============================================================================

source("analysis/00_load_libraries.R")

cat("Starting alpha diversity analysis...\n")

# =============================================================================
# Load Processed Data
# =============================================================================
load("analysis/processed_data.RData")

# =============================================================================
# Calculate Alpha Diversity Metrics
# =============================================================================
cat("Calculating alpha diversity metrics...\n")

# Shannon diversity
shannon <- diversity(comm, index = "shannon")

# Chao1 richness estimator (requires integer counts)
comm_int <- round(comm)
storage.mode(comm_int) <- "integer"
chao1_raw <- estimateR(comm_int)
chao1 <- chao1_raw["S.chao1", ]

# Observed richness
observed <- rowSums(comm > 0)

# Simpson diversity
simpson <- diversity(comm, index = "simpson")

# Add to metadata
meta_clean$Shannon <- shannon[meta_clean$SeqID]
meta_clean$Chao1 <- chao1[meta_clean$SeqID]
meta_clean$Observed <- observed[meta_clean$SeqID]
meta_clean$Simpson <- simpson[meta_clean$SeqID]

cat("✓ Alpha diversity metrics calculated\n")

# =============================================================================
# Statistical Testing by Sample Type
# =============================================================================
cat("Performing statistical tests by sample type...\n")

# Function to format p-values
format_p <- function(p) {
  case_when(
    p < 0.0001 ~ "****",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

# Test each metric by sample type
metrics <- c("Shannon", "Chao1", "Observed", "Simpson")
alpha_results <- data.frame()

for (metric in metrics) {
  for (st in SAMPLE_TYPE_LEVELS) {
    sub <- meta_clean %>% filter(Sample_Type == SAMPLE_TYPE_LABELS[which(SAMPLE_TYPE_LEVELS == st)])

    if (nrow(sub) > 0 && length(unique(sub$Treatment)) == 2) {
      # Wilcoxon rank-sum test
      test_result <- wilcox.test(sub[[metric]] ~ sub$Treatment)

      # Effect size (Cohen's d approximation)
      treated_vals <- sub[[metric]][sub$Treatment == "Treated"]
      untreated_vals <- sub[[metric]][sub$Treatment == "Untreated"]

      pooled_sd <- sqrt(((length(treated_vals) - 1) * var(treated_vals) +
                        (length(untreated_vals) - 1) * var(untreated_vals)) /
                       (length(treated_vals) + length(untreated_vals) - 2))

      cohens_d <- (mean(treated_vals) - mean(untreated_vals)) / pooled_sd

      # Store results
      alpha_results <- rbind(alpha_results, data.frame(
        Sample_Type = st,
        Metric = metric,
        p_value = test_result$p.value,
        significance = format_p(test_result$p.value),
        cohens_d = cohens_d,
        treated_mean = mean(treated_vals),
        untreated_mean = mean(untreated_vals),
        treated_n = length(treated_vals),
        untreated_n = length(untreated_vals)
      ))
    }
  }
}

# Apply multiple testing correction
alpha_results$p_adjusted <- p.adjust(alpha_results$p_value, method = "BH")
alpha_results$significance_adj <- format_p(alpha_results$p_adjusted)

print(alpha_results)

# =============================================================================
# Summary Statistics
# =============================================================================
cat("\nGenerating summary statistics...\n")

alpha_summary <- meta_clean %>%
  group_by(Sample_Type, Treatment) %>%
  summarise(
    across(c(Shannon, Chao1, Observed, Simpson),
           list(mean = mean, sd = sd, median = median),
           na.rm = TRUE, .names = "{.col}_{.fn}"),
    n = n(),
    .groups = "drop"
  )

print(alpha_summary)

# =============================================================================
# Effect Size Analysis
# =============================================================================
cat("\nAnalyzing effect sizes...\n")

# Calculate percent reduction for significant comparisons
effect_summary <- alpha_results %>%
  filter(p_adjusted < 0.05) %>%
  mutate(
    percent_change = ((treated_mean - untreated_mean) / untreated_mean) * 100,
    interpretation = case_when(
      abs(cohens_d) < 0.2 ~ "Small effect",
      abs(cohens_d) < 0.5 ~ "Medium effect",
      abs(cohens_d) < 0.8 ~ "Large effect",
      TRUE ~ "Very large effect"
    )
  ) %>%
  select(Sample_Type, Metric, percent_change, cohens_d, interpretation,
         p_adjusted, significance_adj)

if (nrow(effect_summary) > 0) {
  print(effect_summary)
} else {
  cat("No significant differences detected after multiple testing correction.\n")
}

# =============================================================================
# ESKAPEE Pathogen Analysis
# =============================================================================
cat("\nAnalyzing ESKAPEE pathogen richness...\n")

# Identify ESKAPEE genera in dataset
eskapee_in_data <- intersect(ESKAPEE_GENERA, colnames(comm))
cat("ESKAPEE genera detected:", length(eskapee_in_data), "of", length(ESKAPEE_GENERA), "\n")
cat("Detected:", paste(eskapee_in_data, collapse = ", "), "\n")

if (length(eskapee_in_data) > 0) {
  # Calculate ESKAPEE richness and abundance
  eskapee_richness <- rowSums(comm[, eskapee_in_data] > 0)
  eskapee_abundance <- rowSums(comm_rel[, eskapee_in_data])

  meta_clean$ESKAPEE_richness <- eskapee_richness
  meta_clean$ESKAPEE_abundance <- eskapee_abundance

  # Test ESKAPEE patterns
  eskapee_results <- data.frame()

  for (st in SAMPLE_TYPE_LEVELS) {
    sub <- meta_clean %>% filter(Sample_Type == SAMPLE_TYPE_LABELS[which(SAMPLE_TYPE_LEVELS == st)])

    if (nrow(sub) > 0 && length(unique(sub$Treatment)) == 2) {
      # Test richness
      rich_test <- wilcox.test(sub$ESKAPEE_richness ~ sub$Treatment)
      # Test abundance
      abund_test <- wilcox.test(sub$ESKAPEE_abundance ~ sub$Treatment)

      eskapee_results <- rbind(eskapee_results, data.frame(
        Sample_Type = st,
        Metric = c("ESKAPEE_richness", "ESKAPEE_abundance"),
        p_value = c(rich_test$p.value, abund_test$p.value),
        significance = c(format_p(rich_test$p.value), format_p(abund_test$p.value))
      ))
    }
  }

  print(eskapee_results)
}

# =============================================================================
# Save Results
# =============================================================================
cat("\nSaving alpha diversity results...\n")

# Save statistical results
write.csv(alpha_results, "analysis/alpha_diversity_statistics.csv", row.names = FALSE)
write.csv(alpha_summary, "analysis/alpha_diversity_summary.csv", row.names = FALSE)

if (exists("effect_summary") && nrow(effect_summary) > 0) {
  write.csv(effect_summary, "analysis/alpha_diversity_effects.csv", row.names = FALSE)
}

if (exists("eskapee_results")) {
  write.csv(eskapee_results, "analysis/eskapee_alpha_diversity.csv", row.names = FALSE)
}

# Save updated metadata with alpha diversity metrics
save(meta_clean, comm, comm_rel, genus_table, alpha_results, alpha_summary,
     file = "analysis/alpha_diversity_data.RData")

cat("✓ Results saved to analysis/ directory\n")

# =============================================================================
# Key Findings Summary
# =============================================================================
cat("\n=== ALPHA DIVERSITY ANALYSIS SUMMARY ===\n")

# Significant results
sig_results <- alpha_results %>% filter(p_adjusted < 0.05)

if (nrow(sig_results) > 0) {
  cat("Significant treatment effects detected:\n")
  for (i in 1:nrow(sig_results)) {
    result <- sig_results[i, ]
    direction <- ifelse(result$cohens_d < 0, "reduced", "increased")
    cat(sprintf("  %s %s in %s: %.1f%% %s (p = %.3f)\n",
                result$Metric,
                direction,
                result$Sample_Type,
                abs(((result$treated_mean - result$untreated_mean) / result$untreated_mean) * 100),
                direction,
                result$p_adjusted))
  }
} else {
  cat("No significant treatment effects detected after multiple testing correction.\n")
}

# Sample type patterns
cat("\nSample type patterns:\n")
for (st_label in SAMPLE_TYPE_LABELS) {
  st_results <- alpha_results %>% filter(Sample_Type == SAMPLE_TYPE_LEVELS[which(SAMPLE_TYPE_LABELS == st_label)])
  sig_count <- sum(st_results$p_adjusted < 0.05)
  cat(sprintf("  %s: %d of %d metrics significantly affected\n", st_label, sig_count, nrow(st_results)))
}

cat("\n✓ Alpha diversity analysis completed!\n")