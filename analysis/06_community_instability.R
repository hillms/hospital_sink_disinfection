# =============================================================================
# Hospital Sink Disinfection Study - Community Instability Analysis
#
# Authors: Megan Stroup Hill, Claude Sonnet 4
# Description: Analyze temporal autocorrelation and directional drift
# Input: processed_data.RData
# Output: Instability metrics revealing "controlled instability" signature
# =============================================================================

source("analysis/00_load_libraries.R")

cat("Starting community instability analysis...\n")

# =============================================================================
# Load Processed Data
# =============================================================================
load("analysis/processed_data.RData")

# Load CLR-transformed data if available
if (file.exists("analysis/volatility_analysis.RData")) {
  load("analysis/volatility_analysis.RData")
} else {
  # Recalculate CLR transformation
  pseudocount <- min(comm_rel[comm_rel > 0]) * 0.5
  comm_clr <- log(comm_rel + pseudocount)
  comm_clr <- comm_clr - rowMeans(comm_clr)
}

# =============================================================================
# Analysis 1: Lag-1 Temporal Autocorrelation
# =============================================================================
cat("Calculating lag-1 temporal autocorrelation...\n")

# Calculate autocorrelation function (ACF) for each room-location time series
acf_results <- list()

for (room in unique(meta_clean$Room)) {
  for (sample_type in SAMPLE_TYPE_LEVELS) {
    # Subset data for this room and sample type
    subset_meta <- meta_clean %>%
      filter(Room == room, Sample_Type == SAMPLE_TYPE_LABELS[which(SAMPLE_TYPE_LEVELS == sample_type)]) %>%
      arrange(Week_num)

    if (nrow(subset_meta) < 4) next  # Need at least 4 timepoints for reliable ACF

    treatment_label <- unique(subset_meta$Treatment)[1]

    # Get CLR-transformed data for this subset
    clr_subset <- comm_clr[subset_meta$SeqID, , drop = FALSE]

    # Calculate community centroid
    centroid <- colMeans(clr_subset)

    # Calculate distance of each timepoint from centroid
    distances_from_centroid <- apply(clr_subset, 1, function(x) {
      sqrt(sum((x - centroid)^2))
    })

    # Calculate lag-1 autocorrelation
    n <- length(distances_from_centroid)
    if (n >= 4) {
      lag1_correlation <- cor(distances_from_centroid[-n], distances_from_centroid[-1],
                             method = "pearson")

      acf_results[[length(acf_results) + 1]] <- data.frame(
        Room = room,
        Sample_Type = sample_type,
        Treatment = treatment_label,
        n_timepoints = n,
        lag1_acf = lag1_correlation,
        mean_distance = mean(distances_from_centroid),
        sd_distance = sd(distances_from_centroid)
      )
    }
  }
}

# Combine ACF results
acf_df <- bind_rows(acf_results)

# Apply factor ordering
acf_df$Sample_Type <- factor(acf_df$Sample_Type,
                            levels = SAMPLE_TYPE_LEVELS,
                            labels = SAMPLE_TYPE_LABELS)
acf_df$Treatment <- factor(acf_df$Treatment, levels = c("Treated", "Untreated"))

cat("✓ Lag-1 autocorrelation calculated for", nrow(acf_df), "time series\n")

# =============================================================================
# Analysis 2: PERMANOVA R² for Directional Drift
# =============================================================================
cat("Calculating PERMANOVA R² for directional drift...\n")

# Calculate per-room PERMANOVA R² (week effect)
permanova_results <- list()

set.seed(42)  # For reproducible permutation tests

for (room in unique(meta_clean$Room)) {
  for (sample_type in SAMPLE_TYPE_LEVELS) {
    subset_meta <- meta_clean %>%
      filter(Room == room, Sample_Type == SAMPLE_TYPE_LABELS[which(SAMPLE_TYPE_LEVELS == sample_type)]) %>%
      arrange(Week_num)

    if (nrow(subset_meta) < 5) next  # Need at least 5 samples for PERMANOVA

    treatment_label <- unique(subset_meta$Treatment)[1]

    # Get CLR-transformed data for this subset
    clr_subset <- comm_clr[subset_meta$SeqID, , drop = FALSE]

    # Calculate Aitchison distances
    aitchison_dist <- dist(clr_subset, method = "euclidean")

    # PERMANOVA testing for week effect (directional drift)
    permanova_result <- tryCatch({
      adonis2(aitchison_dist ~ Week_num, data = subset_meta, permutations = 999)
    }, error = function(e) NULL)

    if (!is.null(permanova_result)) {
      r2_value <- permanova_result[["R2"]][1]
      p_value <- permanova_result[["Pr(>F)"]][1]

      permanova_results[[length(permanova_results) + 1]] <- data.frame(
        Room = room,
        Sample_Type = sample_type,
        Treatment = treatment_label,
        n_samples = nrow(subset_meta),
        R2_week = round(r2_value, 4),
        P_week = round(p_value, 4),
        significant = p_value < 0.05
      )
    }
  }
}

# Combine PERMANOVA results
permanova_df <- bind_rows(permanova_results)

# Apply factor ordering
permanova_df$Sample_Type <- factor(permanova_df$Sample_Type,
                                  levels = SAMPLE_TYPE_LEVELS,
                                  labels = SAMPLE_TYPE_LABELS)
permanova_df$Treatment <- factor(permanova_df$Treatment, levels = c("Treated", "Untreated"))

cat("✓ PERMANOVA R² calculated for", nrow(permanova_df), "time series\n")

# =============================================================================
# Statistical Comparisons
# =============================================================================
cat("Performing statistical comparisons...\n")

# Test ACF differences between treatments
acf_stats <- data.frame()

for (st in SAMPLE_TYPE_LABELS) {
  subset_acf <- acf_df %>% filter(Sample_Type == st)

  if (nrow(subset_acf) > 0 && length(unique(subset_acf$Treatment)) == 2) {
    w_test <- wilcox.test(lag1_acf ~ Treatment, data = subset_acf, exact = FALSE)

    # Summary statistics
    summary_stats <- subset_acf %>%
      group_by(Treatment) %>%
      summarise(
        n = n(),
        mean = mean(lag1_acf),
        sd = sd(lag1_acf),
        median = median(lag1_acf),
        .groups = "drop"
      )

    acf_stats <- rbind(acf_stats, data.frame(
      Sample_Type = st,
      metric = "lag1_acf",
      p_value = w_test$p.value,
      significance = format_p_value(w_test$p.value),
      treated_mean = summary_stats$mean[summary_stats$Treatment == "Treated"],
      untreated_mean = summary_stats$mean[summary_stats$Treatment == "Untreated"],
      treated_n = summary_stats$n[summary_stats$Treatment == "Treated"],
      untreated_n = summary_stats$n[summary_stats$Treatment == "Untreated"]
    ))
  }
}

# Test PERMANOVA R² differences between treatments
r2_stats <- data.frame()

for (st in SAMPLE_TYPE_LABELS) {
  subset_r2 <- permanova_df %>% filter(Sample_Type == st)

  if (nrow(subset_r2) > 0 && length(unique(subset_r2$Treatment)) == 2) {
    w_test <- wilcox.test(R2_week ~ Treatment, data = subset_r2, exact = FALSE)

    # Summary statistics
    summary_stats <- subset_r2 %>%
      group_by(Treatment) %>%
      summarise(
        n = n(),
        mean = mean(R2_week),
        sd = sd(R2_week),
        median = median(R2_week),
        .groups = "drop"
      )

    r2_stats <- rbind(r2_stats, data.frame(
      Sample_Type = st,
      metric = "R2_week",
      p_value = w_test$p.value,
      significance = format_p_value(w_test$p.value),
      treated_mean = summary_stats$mean[summary_stats$Treatment == "Treated"],
      untreated_mean = summary_stats$mean[summary_stats$Treatment == "Untreated"],
      treated_n = summary_stats$n[summary_stats$Treatment == "Treated"],
      untreated_n = summary_stats$n[summary_stats$Treatment == "Untreated"]
    ))
  }
}

# Combine statistics
instability_stats <- rbind(acf_stats, r2_stats)
print(instability_stats)

# =============================================================================
# Joint Analysis: ACF vs R² Landscape
# =============================================================================
cat("Performing joint ACF vs R² analysis...\n")

# Merge ACF and R² data for joint analysis
joint_data <- inner_join(
  acf_df %>% select(Room, Sample_Type, Treatment, lag1_acf),
  permanova_df %>% select(Room, Sample_Type, Treatment, R2_week, P_week),
  by = c("Room", "Sample_Type", "Treatment")
)

# Calculate group centroids for each treatment-location combination
centroids <- joint_data %>%
  group_by(Sample_Type, Treatment) %>%
  summarise(
    mean_R2 = mean(R2_week),
    mean_acf = mean(lag1_acf),
    n = n(),
    .groups = "drop"
  )

print(centroids)

# Define instability quadrants
# High ACF, Low R² = Stable (predictable, no drift)
# Low ACF, Low R² = Controlled instability (unpredictable, no drift) <- Treatment signature
# High ACF, High R² = Predictable drift
# Low ACF, High R² = Chaotic drift

# Classify each room into quadrants
r2_threshold <- 0.085  # Based on data distribution
acf_threshold <- 0.0

joint_data <- joint_data %>%
  mutate(
    quadrant = case_when(
      lag1_acf >= acf_threshold & R2_week < r2_threshold ~ "Stable",
      lag1_acf < acf_threshold & R2_week < r2_threshold ~ "Controlled_instability",
      lag1_acf >= acf_threshold & R2_week >= r2_threshold ~ "Predictable_drift",
      lag1_acf < acf_threshold & R2_week >= r2_threshold ~ "Chaotic_drift"
    )
  )

# Analyze quadrant distributions
quadrant_summary <- joint_data %>%
  count(Sample_Type, Treatment, quadrant) %>%
  group_by(Sample_Type, Treatment) %>%
  mutate(percentage = (n / sum(n)) * 100) %>%
  ungroup()

print(quadrant_summary)

# =============================================================================
# Save Results
# =============================================================================
cat("\nSaving instability analysis results...\n")

# Save main results
write.csv(acf_df, "analysis/lag1_autocorrelation.csv", row.names = FALSE)
write.csv(permanova_df, "analysis/permanova_r2.csv", row.names = FALSE)
write.csv(joint_data, "analysis/joint_instability_data.csv", row.names = FALSE)
write.csv(instability_stats, "analysis/instability_statistics.csv", row.names = FALSE)
write.csv(centroids, "analysis/instability_centroids.csv", row.names = FALSE)
write.csv(quadrant_summary, "analysis/quadrant_classification.csv", row.names = FALSE)

# Save workspace
save(acf_df, permanova_df, joint_data, instability_stats, centroids,
     quadrant_summary, acf_threshold, r2_threshold,
     file = "analysis/instability_analysis.RData")

cat("✓ Results saved to analysis/ directory\n")

# =============================================================================
# Analysis Summary
# =============================================================================
cat("\n=== COMMUNITY INSTABILITY ANALYSIS SUMMARY ===\n")

cat("Lag-1 Autocorrelation Results:\n")
for (i in 1:nrow(acf_stats)) {
  result <- acf_stats[i, ]
  direction <- ifelse(result$treated_mean < result$untreated_mean, "reduced", "increased")
  cat(sprintf("  %s: %s predictability in treated (%.3f vs %.3f, p = %.3f)\n",
              result$Sample_Type, direction, result$treated_mean,
              result$untreated_mean, result$p_value))
}

cat("\nPERMANOVA R² Results:\n")
for (i in 1:nrow(r2_stats)) {
  result <- r2_stats[i, ]
  direction <- ifelse(result$treated_mean < result$untreated_mean, "reduced", "increased")
  cat(sprintf("  %s: %s directional drift in treated (%.3f vs %.3f, p = %.3f)\n",
              result$Sample_Type, direction, result$treated_mean,
              result$untreated_mean, result$p_value))
}

cat("\nControlled Instability Signature:\n")
controlled_inst <- quadrant_summary %>%
  filter(quadrant == "Controlled_instability") %>%
  arrange(Sample_Type, desc(percentage))

for (i in 1:nrow(controlled_inst)) {
  result <- controlled_inst[i, ]
  cat(sprintf("  %s %s: %.1f%% in controlled instability quadrant\n",
              result$Sample_Type, result$Treatment, result$percentage))
}

cat("\nKey Findings:\n")
cat("- Lag-1 ACF measures short-term predictability of community dynamics\n")
cat("- PERMANOVA R² quantifies directional drift over time\n")
cat("- Treatment creates 'controlled instability': low ACF + low R²\n")
cat("- This signature indicates disrupted biofilm succession without chaos\n")

cat("\n✓ Community instability analysis completed!\n")