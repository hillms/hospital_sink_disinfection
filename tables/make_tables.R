# =============================================================================
# Hospital Sink Disinfection Study - Manuscript Tables
#
# Authors: Megan Stroup Hill, Claude
# Description: Build Table 1 and Tables S1-S3 from the analysis outputs
# Input: analysis/aldex2_results_raw_counts.csv   (04_differential_abundance.R)
#        result/decontam/decontam_scores_all_ASVs.tsv (preprocessing/02_decontam.R)
#        analysis/da_calibration_*.csv            (04b_da_calibration.R)
#        analysis/weekly_permanova_pvalues.csv    (07_weekly_permanova.R)
# Output: tables/output/Table_1.csv, Table_S1.csv, Table_S2.csv, Table_S3.csv
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(data.table)
})

cat("Building manuscript tables...\n")
dir.create("tables/output", showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# Table 1: significant genera, ALDEx2 on raw counts (* = contains an ESKAPEE pathogen)
# =============================================================================
res <- read.csv("analysis/aldex2_results_raw_counts.csv")
table1 <- res %>%
  filter(significant, Sample_Type %in% c("PT", "TP")) %>%
  mutate(Location = factor(c(PT = "P-Trap", TP = "Tail Piece")[Sample_Type],
                           levels = c("P-Trap", "Tail Piece")),
         Genus = paste0(genus, ifelse(is_eskapee, " *", ""))) %>%
  arrange(Location, effect) %>%
  transmute(Location, Genus, Effect = round(effect, 3),
            `BH(wi)` = signif(wi.eBH, 3), `BH(we)` = signif(we.eBH, 3),
            diff.btw = round(diff.btw, 3),
            Direction = ifelse(effect < 0, "TREATED", "UNTREATED"))
fwrite(table1, "tables/output/Table_1.csv")
cat("✓ Table 1:", nrow(table1), "genera\n")

# =============================================================================
# Table S1: ASVs flagged by decontam (P < 0.1) and whether they were removed
# =============================================================================
s <- fread("result/decontam/decontam_scores_all_ASVs.tsv")
n_ctrl <- c(TP = 4, PT = 3, TOP = 5)
n_sink <- max(s$n_true_present / s$prev_true, na.rm = TRUE)
tot <- sum(s$reads_true, na.rm = TRUE)
tableS1 <- s[contaminant_t0.1 == TRUE][order(-remove, -reads_true), .(
  `ASV ID` = ASV,
  `Taxonomy (GTDB r220)` = gsub("^d__Bacteria; |[a-z]__", "", taxonomy),
  `decontam score` = signif(p.prev, 2),
  `Negative controls (of 12)` = n_neg_present,
  `Tail piece controls (of 4)` = round(prev_neg_TP * n_ctrl["TP"]),
  `P-trap controls (of 3)` = round(prev_neg_PT * n_ctrl["PT"]),
  `Countertop controls (of 5)` = round(prev_neg_TOP * n_ctrl["TOP"]),
  `Sink samples detected` = n_true_present,
  `Reads in negative controls` = reads_neg,
  `Reads in sink samples` = reads_true,
  `% of sink sample reads` = sprintf("%.3f", 100 * reads_true / tot),
  Removed = ifelse(remove, "Yes", "No"))]
fwrite(tableS1, "tables/output/Table_S1.csv")
cat("✓ Table S1:", nrow(tableS1), "flagged ASVs,", sum(tableS1$Removed == "Yes"), "removed\n")

# =============================================================================
# Table S2: calibration of ALDEx2, DESeq2 and edgeR with room-shuffled labels
# =============================================================================
if (file.exists("analysis/da_calibration_shuffles.csv")) {
  true_labels <- read.csv("analysis/da_calibration_true_labels.csv")
  shuffles <- read.csv("analysis/da_calibration_shuffles.csv")
  tableS2 <- shuffles %>%
    group_by(Method, Sample_Type) %>%
    summarise(med = median(n_sig), lo = min(n_sig), hi = max(n_sig),
              each = paste(n_sig, collapse = ", "), .groups = "drop") %>%
    left_join(true_labels, by = c("Method", "Sample_Type")) %>%
    mutate(Sample_Type = factor(Sample_Type, levels = c("TP", "PT", "TOP"))) %>%
    arrange(Method, Sample_Type) %>%
    transmute(Method,
              Location = c(TP = "Tail piece", PT = "P-trap", TOP = "Countertop")[as.character(Sample_Type)],
              `Genera tested` = n_tested,
              `Significant genera, true labels` = n_sig,
              `Significant genera, shuffled labels: median (range)` = sprintf("%g (%d-%d)", med, lo, hi),
              `Significant genera in each of 5 shuffles` = each)
  fwrite(tableS2, "tables/output/Table_S2.csv")
  cat("✓ Table S2\n")
} else {
  cat("⚠ Table S2 skipped: run analysis/04b_da_calibration.R first\n")
}

# =============================================================================
# Table S3: weekly PERMANOVA P values (Bray-Curtis, treatment as predictor)
# =============================================================================
weekly <- read.csv("analysis/weekly_permanova_pvalues.csv")
stars <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))
tableS3 <- weekly %>%
  mutate(P = paste0(sprintf("%.3f", P_value), stars(P_value)),
         # column order as in the published table
         Sample_Type = factor(Sample_Type, levels = c("P-Trap", "Tail Piece", "Countertop"))) %>%
  select(Week, Sample_Type, P) %>%
  arrange(Sample_Type) %>%
  pivot_wider(names_from = Sample_Type, values_from = P, values_fill = "-") %>%
  arrange(Week)
fwrite(tableS3, "tables/output/Table_S3.csv")
cat("✓ Table S3:", nrow(tableS3), "weeks\n")

cat("\n✓ Tables saved to tables/output/\n")
