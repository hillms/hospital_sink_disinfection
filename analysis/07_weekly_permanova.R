# =============================================================================
# Hospital Sink Disinfection Study - Weekly Treatment PERMANOVA (Table S3)
#
# Authors: Megan Stroup Hill, Claude
# Description: For each week and location, PERMANOVA of Bray-Curtis distances
#              with treatment as the predictor
# Input: processed_data.RData
# Output: analysis/weekly_permanova_pvalues.csv (long),
#         analysis/weekly_permanova_summary.csv (significant weeks per location)
#
# Weeks are tested only when both treatments have samples (>= 4 samples in
# total), so the number of weeks tested differs by location (30 TP, 28 PT,
# 30 TOP). Significant-week counts are computed from the results, not set by hand.
# =============================================================================

source("analysis/00_load_libraries.R")

cat("Starting weekly PERMANOVA...\n")

load("analysis/processed_data.RData")

weekly <- list()
set.seed(42)
for (week in sort(unique(meta_clean$Week_num))) {
  for (st in SAMPLE_TYPE_LABELS) {
    wm <- meta_clean %>% filter(Week_num == week, Sample_Type == st)
    if (nrow(wm) < 4 || length(unique(wm$Treatment)) < 2) next
    wc <- comm[wm$SeqID, , drop = FALSE]
    wm <- wm[rowSums(wc) > 0, ]
    wc <- wc[wm$SeqID, , drop = FALSE]
    if (nrow(wc) < 4 || length(unique(wm$Treatment)) < 2) next

    res <- tryCatch(adonis2(vegdist(wc, method = "bray") ~ Treatment, data = wm,
                            permutations = 999, by = "margin"),
                    error = function(e) NULL)
    if (is.null(res)) next
    weekly[[length(weekly) + 1]] <- data.frame(
      Week = week, Sample_Type = st, R2 = res$R2[1], P_value = res$`Pr(>F)`[1],
      n_treated = sum(wm$Treatment == "Treated"),
      n_untreated = sum(wm$Treatment == "Untreated"))
  }
}
weekly <- bind_rows(weekly)

summary_weeks <- weekly %>%
  group_by(Sample_Type) %>%
  summarise(weeks_tested = n(), significant_weeks = sum(P_value < 0.05),
            first_week_P = P_value[which.min(Week)], .groups = "drop")

write.csv(weekly, "analysis/weekly_permanova_pvalues.csv", row.names = FALSE)
write.csv(summary_weeks, "analysis/weekly_permanova_summary.csv", row.names = FALSE)

cat("\n=== WEEKLY PERMANOVA SUMMARY ===\n")
print(summary_weeks)
cat("\n✓ Weekly PERMANOVA completed!\n")
