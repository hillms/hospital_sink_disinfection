# =============================================================================
# Hospital Sink Disinfection Study - Alpha Diversity Analysis
#
# Authors: Megan Stroup Hill, Claude
# Description: Observed genera and Shannon diversity on raw genus counts
#              rarefied to 1,000 reads per sample
# Input: processed_data.RData
# Output: Alpha diversity metrics and statistical comparisons
#
# Read depth is strongly confounded with treatment (median tail-piece reads
# 1,459 treated vs 67,339 untreated), so richness from unrarefied or
# depth-normalized counts largely tracks depth. DADA2 output also has no
# singletons, so Chao1 equals observed richness; observed genera is reported.
# Samples with < 1,000 reads are excluded from these tests.
# =============================================================================

source("analysis/00_load_libraries.R")

cat("Starting alpha diversity analysis...\n")

RAREFY_DEPTH <- 1000

# =============================================================================
# Load Processed Data
# =============================================================================
load("analysis/processed_data.RData")

# =============================================================================
# Rarefy and Calculate Alpha Diversity Metrics
# =============================================================================
cat("Rarefying raw genus counts to", RAREFY_DEPTH, "reads per sample...\n")

depth <- rowSums(comm_raw)
ok <- depth >= RAREFY_DEPTH
set.seed(42)
comm_rare <- rrarefy(comm_raw[ok, ], RAREFY_DEPTH)

meta_clean$Depth <- depth[meta_clean$SeqID]
meta_clean$Observed <- rowSums(comm_rare > 0)[meta_clean$SeqID]
meta_clean$Shannon <- vegan::diversity(comm_rare, index = "shannon")[meta_clean$SeqID]

cat("✓ Samples rarefied:", sum(ok), "of", length(ok), "\n")
print(meta_clean %>%
        group_by(Sample_Type, Treatment) %>%
        summarise(n = n(), median_reads = median(Depth),
                  n_below_depth = sum(Depth < RAREFY_DEPTH), .groups = "drop"))

# =============================================================================
# Statistical Testing by Sample Type
# =============================================================================
cat("\nPerforming Wilcoxon tests by sample type (BH across 6 tests)...\n")

metrics <- c("Observed", "Shannon")
alpha_results <- data.frame()

for (metric in metrics) {
  for (st in SAMPLE_TYPE_LABELS) {
    sub <- meta_clean %>% filter(Sample_Type == st, !is.na(.data[[metric]]))
    test_result <- wilcox.test(sub[[metric]] ~ sub$Treatment)
    treated_vals <- sub[[metric]][sub$Treatment == "Treated"]
    untreated_vals <- sub[[metric]][sub$Treatment == "Untreated"]

    alpha_results <- rbind(alpha_results, data.frame(
      Sample_Type = st,
      Metric = metric,
      p_value = test_result$p.value,
      treated_mean = mean(treated_vals),
      untreated_mean = mean(untreated_vals),
      percent_change = 100 * (mean(treated_vals) - mean(untreated_vals)) / mean(untreated_vals),
      treated_n = length(treated_vals),
      untreated_n = length(untreated_vals)
    ))
  }
}

alpha_results$p_adjusted <- p.adjust(alpha_results$p_value, method = "BH")
alpha_results$significance_adj <- format_p_value(alpha_results$p_adjusted)

print(alpha_results)

# =============================================================================
# Save Results
# =============================================================================
cat("\nSaving alpha diversity results...\n")

write.csv(alpha_results, "analysis/alpha_diversity_statistics.csv", row.names = FALSE)
write.csv(meta_clean %>% select(SeqID, Sample_Type, Treatment, Depth, Observed, Shannon),
          "analysis/alpha_diversity_per_sample.csv", row.names = FALSE)

save(meta_clean, alpha_results, RAREFY_DEPTH, file = "analysis/alpha_diversity_data.RData")

cat("✓ Results saved to analysis/ directory\n")

# =============================================================================
# Key Findings Summary
# =============================================================================
cat("\n=== ALPHA DIVERSITY ANALYSIS SUMMARY ===\n")
for (i in seq_len(nrow(alpha_results))) {
  r <- alpha_results[i, ]
  cat(sprintf("  %s %s: treated %.1f vs untreated %.1f (%+.0f%%), BH P = %.2g %s\n",
              r$Sample_Type, r$Metric, r$treated_mean, r$untreated_mean,
              r$percent_change, r$p_adjusted, r$significance_adj))
}

cat("\n✓ Alpha diversity analysis completed!\n")
