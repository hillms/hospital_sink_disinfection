# =============================================================================
# Hospital Sink Disinfection Study - Differential Abundance (ALDEx2)
#
# Authors: Megan Stroup Hill, Claude
# Description: ALDEx2 on raw (unrarefied) genus read counts, treated vs
#              untreated, separately for each sample location
# Input: processed_data.RData
# Output: analysis/aldex2_results_raw_counts.csv, analysis/da_shared_genera.csv
#
# Settings: genera present in >= 5% of samples at a location; CLR with 128
# Monte Carlo instances; Welch t-test and Wilcoxon; significant if either
# BH-adjusted P < 0.05. Negative effect = enriched in treated samples.
# DESeq2 and edgeR were also tested; see 04b_da_calibration.R for why only
# ALDEx2 is reported.
# =============================================================================

source("analysis/00_load_libraries.R")

cat("Starting differential abundance analysis (ALDEx2, raw counts)...\n")

load("analysis/processed_data.RData")

# =============================================================================
# ALDEx2 per location
# =============================================================================
# Genus names that contain an ESKAPEE pathogen (Table 1 *, bold in Figure 5)
ESKAPEE_CANONICAL <- c("Enterococcus", "Staphylococcus", "Klebsiella", "Acinetobacter",
                       "Pseudomonas", "Enterobacter", "Escherichia")

run_aldex <- function(st, treatment = meta_clean$Treatment) {
  idx <- meta_clean$Sample_Type == SAMPLE_TYPE_LABELS[SAMPLE_TYPE_LEVELS == st]
  sc <- t(comm_raw[idx, ])
  sc <- sc[rowSums(sc > 0) / ncol(sc) >= 0.05, ]
  set.seed(42)
  a <- aldex(sc, as.character(treatment[idx]), mc.samples = 128, test = "t",
             effect = TRUE, verbose = FALSE)
  data.frame(genus = rownames(a), Sample_Type = st, a, row.names = NULL)
}

aldex_results <- bind_rows(lapply(SAMPLE_TYPE_LEVELS, function(st) {
  r <- run_aldex(st)
  cat(sprintf("  %s: %d genera tested, %d significant\n", st, nrow(r),
              sum(r$wi.eBH < 0.05 | r$we.eBH < 0.05)))
  r
}))

aldex_results <- aldex_results %>%
  mutate(significant = wi.eBH < 0.05 | we.eBH < 0.05,
         direction = ifelse(effect < 0, "Treated-enriched", "Untreated-enriched"),
         is_eskapee = genus %in% ESKAPEE_CANONICAL)

# =============================================================================
# Genera shifted in the same direction in tail pieces and P-traps
# =============================================================================
shared <- aldex_results %>%
  filter(Sample_Type %in% c("TP", "PT")) %>%
  select(genus, Sample_Type, effect, significant) %>%
  pivot_wider(names_from = Sample_Type, values_from = c(effect, significant)) %>%
  filter(significant_TP %in% TRUE, significant_PT %in% TRUE,
         sign(effect_TP) == sign(effect_PT))

cat("\nGenera significant in both tail pieces and P-traps, same direction:", nrow(shared), "\n")
cat("  (compare with room-shuffled labels in 04b_da_calibration.R)\n")

# =============================================================================
# Save Results
# =============================================================================
write.csv(aldex_results, "analysis/aldex2_results_raw_counts.csv", row.names = FALSE)
write.csv(shared, "analysis/da_shared_genera.csv", row.names = FALSE)
save(aldex_results, shared, run_aldex, file = "analysis/differential_abundance.RData")

cat("\n=== DIFFERENTIAL ABUNDANCE SUMMARY ===\n")
print(aldex_results %>%
        group_by(Sample_Type) %>%
        summarise(tested = n(),
                  treated_enriched = sum(significant & effect < 0),
                  untreated_enriched = sum(significant & effect > 0),
                  significant = sum(significant)))
cat("\nESKAPEE genera:\n")
print(aldex_results %>% filter(is_eskapee, Sample_Type != "TOP") %>%
        transmute(Sample_Type, genus, effect = round(effect, 3),
                  BH_P = signif(pmin(wi.eBH, we.eBH), 2), significant))

cat("\n✓ Differential abundance analysis completed!\n")
