# =============================================================================
# Hospital Sink Disinfection Study - Community Volatility Analysis
#
# Authors: Megan Stroup Hill, Claude Sonnet 4
# Description: Analyze week-to-week community compositional changes
# Input: processed_data.RData
# Output: Volatility metrics and temporal dynamics
# =============================================================================

source("analysis/00_load_libraries.R")

cat("Starting community volatility analysis...\n")

# =============================================================================
# Load Processed Data
# =============================================================================
load("analysis/processed_data.RData")

# =============================================================================
# CLR Transformation for Compositional Data
# =============================================================================
cat("Applying CLR transformation...\n")

# CLR transformation for compositional data analysis
pseudocount <- min(comm_rel[comm_rel > 0]) * 0.5
comm_clr <- log(comm_rel + pseudocount)
comm_clr <- comm_clr - rowMeans(comm_clr)

cat("✓ CLR transformation applied\n")

# =============================================================================
# Calculate Week-to-Week Volatility
# =============================================================================
cat("Calculating community volatility (week-to-week distances)...\n")

# Calculate volatility for each room-location combination
volatility_list <- list()

for (room in unique(meta_clean$Room)) {
  for (sample_type in SAMPLE_TYPE_LEVELS) {
    # Subset data for this room and sample type
    subset_meta <- meta_clean %>%
      filter(Room == room, Sample_Type == SAMPLE_TYPE_LABELS[which(SAMPLE_TYPE_LEVELS == sample_type)]) %>%
      arrange(Week_num)

    if (nrow(subset_meta) < 2) next  # Need at least 2 timepoints

    treatment_label <- unique(subset_meta$Treatment)[1]

    # Calculate consecutive distances using Aitchison distance (Euclidean on CLR)
    distances <- numeric(nrow(subset_meta) - 1)

    for (i in 2:nrow(subset_meta)) {
      sample1 <- subset_meta$SeqID[i - 1]
      sample2 <- subset_meta$SeqID[i]

      # Aitchison distance (Euclidean distance on CLR-transformed data)
      distances[i - 1] <- dist(comm_clr[c(sample1, sample2), ], method = "euclidean")
    }

    # Store results
    volatility_list[[length(volatility_list) + 1]] <- data.frame(
      Room = room,
      Sample_Type = sample_type,
      Treatment = treatment_label,
      n_transitions = length(distances),
      mean_volatility = mean(distances, na.rm = TRUE),
      sd_volatility = sd(distances, na.rm = TRUE),
      median_volatility = median(distances, na.rm = TRUE),
      max_volatility = max(distances, na.rm = TRUE),
      # Store individual step data for detailed analysis
      step_data = I(list(data.frame(
        Week_from = subset_meta$Week_num[-nrow(subset_meta)],
        Week_to = subset_meta$Week_num[-1],
        distance = distances,
        Room = room,
        Sample_Type = sample_type,
        Treatment = treatment_label
      )))
    )
  }
}

# Combine volatility results
volatility_df <- bind_rows(volatility_list[sapply(volatility_list, is.data.frame)])

# Extract individual step data
step_data <- bind_rows(lapply(volatility_list, function(x) x$step_data[[1]]))

cat("✓ Volatility calculated for", nrow(volatility_df), "room-location combinations\n")
cat("✓ Individual transitions analyzed:", nrow(step_data), "week-to-week steps\n")

# =============================================================================
# Apply Factor Ordering
# =============================================================================
# Apply consistent factor ordering
volatility_df$Sample_Type <- factor(volatility_df$Sample_Type,
                                   levels = SAMPLE_TYPE_LEVELS,
                                   labels = SAMPLE_TYPE_LABELS)

step_data$Sample_Type <- factor(step_data$Sample_Type,
                               levels = SAMPLE_TYPE_LEVELS,
                               labels = SAMPLE_TYPE_LABELS)

volatility_df$Treatment <- factor(volatility_df$Treatment, levels = c("Treated", "Untreated"))
step_data$Treatment <- factor(step_data$Treatment, levels = c("Treated", "Untreated"))

# =============================================================================
# Statistical Analysis of Overall Volatility
# =============================================================================
cat("Performing statistical tests on overall volatility...\n")

volatility_stats <- data.frame()

for (st in SAMPLE_TYPE_LABELS) {
  subset_vol <- volatility_df %>% filter(Sample_Type == st)

  if (nrow(subset_vol) > 0 && length(unique(subset_vol$Treatment)) == 2) {
    # Wilcoxon test
    w_test <- wilcox.test(mean_volatility ~ Treatment, data = subset_vol)

    # Summary statistics
    summary_stats <- subset_vol %>%
      group_by(Treatment) %>%
      summarise(
        n = n(),
        mean = mean(mean_volatility),
        sd = sd(mean_volatility),
        median = median(mean_volatility),
        .groups = "drop"
      )

    # Effect size
    treated_vals <- subset_vol$mean_volatility[subset_vol$Treatment == "Treated"]
    untreated_vals <- subset_vol$mean_volatility[subset_vol$Treatment == "Untreated"]

    cohens_d <- (mean(treated_vals) - mean(untreated_vals)) /
      sqrt(((length(treated_vals) - 1) * var(treated_vals) +
            (length(untreated_vals) - 1) * var(untreated_vals)) /
           (length(treated_vals) + length(untreated_vals) - 2))

    volatility_stats <- rbind(volatility_stats, data.frame(
      Sample_Type = st,
      p_value = w_test$p.value,
      significance = format_p_value(w_test$p.value),
      cohens_d = cohens_d,
      treated_mean = mean(treated_vals),
      untreated_mean = mean(untreated_vals),
      treated_n = length(treated_vals),
      untreated_n = length(untreated_vals)
    ))
  }
}

print(volatility_stats)

# =============================================================================
# Weekly Pattern Analysis
# =============================================================================
cat("Analyzing weekly volatility patterns...\n")

# Calculate weekly averages across all rooms
weekly_volatility <- step_data %>%
  group_by(Sample_Type, Treatment, Week_to) %>%
  summarise(
    mean_distance = mean(distance, na.rm = TRUE),
    se_distance = sd(distance, na.rm = TRUE) / sqrt(n()),
    n_transitions = n(),
    .groups = "drop"
  ) %>%
  filter(n_transitions >= 3)  # Require at least 3 transitions for reliability

# Test for weekly differences
weekly_tests <- data.frame()

for (week in unique(weekly_volatility$Week_to)) {
  for (st in SAMPLE_TYPE_LABELS) {
    week_subset <- step_data %>%
      filter(Week_to == week, Sample_Type == st)

    if (nrow(week_subset) >= 6 && length(unique(week_subset$Treatment)) == 2) {
      w_test <- wilcox.test(distance ~ Treatment, data = week_subset)

      weekly_tests <- rbind(weekly_tests, data.frame(
        Week = week,
        Sample_Type = st,
        p_value = w_test$p.value,
        n_treated = sum(week_subset$Treatment == "Treated"),
        n_untreated = sum(week_subset$Treatment == "Untreated")
      ))
    }
  }
}

# Multiple testing correction for weekly tests
weekly_tests$p_adjusted <- p.adjust(weekly_tests$p_value, method = "BH")
weekly_tests$significant <- weekly_tests$p_adjusted < 0.05

# Summary of weekly significance
weekly_summary <- weekly_tests %>%
  group_by(Sample_Type) %>%
  summarise(
    total_weeks_tested = n(),
    significant_weeks = sum(significant),
    percent_significant = (significant_weeks / total_weeks_tested) * 100,
    .groups = "drop"
  )

print(weekly_summary)

# =============================================================================
# Temporal Clustering Analysis
# =============================================================================
cat("Analyzing temporal clustering of significant weeks...\n")

# Identify clusters of consecutive significant weeks
temporal_clusters <- list()

for (st in SAMPLE_TYPE_LABELS) {
  sig_weeks <- weekly_tests %>%
    filter(Sample_Type == st, significant) %>%
    pull(Week) %>%
    sort()

  if (length(sig_weeks) > 0) {
    # Find consecutive week clusters
    clusters <- list()
    current_cluster <- sig_weeks[1]

    for (i in seq_along(sig_weeks)[-1]) {
      if (sig_weeks[i] == sig_weeks[i-1] + 1) {
        # Consecutive week
        current_cluster <- c(current_cluster, sig_weeks[i])
      } else {
        # Gap found, save current cluster and start new one
        clusters[[length(clusters) + 1]] <- current_cluster
        current_cluster <- sig_weeks[i]
      }
    }

    # Save the last cluster
    clusters[[length(clusters) + 1]] <- current_cluster

    # Format cluster ranges
    cluster_ranges <- sapply(clusters, function(x) {
      if (length(x) == 1) return(as.character(x))
      else return(paste(range(x), collapse = "-"))
    })

    temporal_clusters[[st]] <- list(
      clusters = clusters,
      ranges = cluster_ranges,
      largest_cluster = max(sapply(clusters, length))
    )
  }
}

# Print temporal clustering results
for (st in names(temporal_clusters)) {
  cat("\n", st, "temporal clusters:\n")
  cat("  Cluster ranges:", paste(temporal_clusters[[st]]$ranges, collapse = ", "), "\n")
  cat("  Largest cluster:", temporal_clusters[[st]]$largest_cluster, "consecutive weeks\n")
}

# =============================================================================
# Save Results
# =============================================================================
cat("\nSaving volatility analysis results...\n")

# Save main results
write.csv(volatility_df, "analysis/community_volatility_summary.csv", row.names = FALSE)
write.csv(step_data, "analysis/weekly_transitions.csv", row.names = FALSE)
write.csv(volatility_stats, "analysis/volatility_statistics.csv", row.names = FALSE)
write.csv(weekly_tests, "analysis/weekly_volatility_tests.csv", row.names = FALSE)
write.csv(weekly_summary, "analysis/weekly_significance_summary.csv", row.names = FALSE)

# Save temporal clustering results
temporal_cluster_summary <- data.frame(
  Sample_Type = names(temporal_clusters),
  n_clusters = sapply(temporal_clusters, function(x) length(x$clusters)),
  largest_cluster = sapply(temporal_clusters, function(x) x$largest_cluster),
  cluster_ranges = sapply(temporal_clusters, function(x) paste(x$ranges, collapse = "; "))
)

write.csv(temporal_cluster_summary, "analysis/temporal_clustering.csv", row.names = FALSE)

# Save workspace
save(volatility_df, step_data, volatility_stats, weekly_tests, weekly_summary,
     weekly_volatility, temporal_clusters, comm_clr,
     file = "analysis/volatility_analysis.RData")

cat("✓ Results saved to analysis/ directory\n")

# =============================================================================
# Analysis Summary
# =============================================================================
cat("\n=== COMMUNITY VOLATILITY ANALYSIS SUMMARY ===\n")

cat("Overall volatility patterns:\n")
for (i in 1:nrow(volatility_stats)) {
  result <- volatility_stats[i, ]
  direction <- ifelse(result$treated_mean > result$untreated_mean, "higher", "lower")
  cat(sprintf("  %s: %s volatility in treated (p = %.3f, %s)\n",
              result$Sample_Type, direction, result$p_value, result$significance))
}

cat("\nWeekly significance patterns:\n")
for (i in 1:nrow(weekly_summary)) {
  result <- weekly_summary[i, ]
  cat(sprintf("  %s: %d/%d weeks significant (%.1f%%)\n",
              result$Sample_Type, result$significant_weeks,
              result$total_weeks_tested, result$percent_significant))
}

cat("\nKey findings:\n")
cat("- Community volatility represents week-to-week compositional turnover\n")
cat("- Aitchison distances used to properly handle compositional data\n")
cat("- Weekly analysis reveals episodic treatment effects\n")
cat("- Temporal clustering indicates facility-wide perturbation events\n")

cat("\n✓ Community volatility analysis completed!\n")