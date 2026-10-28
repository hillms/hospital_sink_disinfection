# =============================================================================
# Hospital Sink Disinfection Study - Contaminant Removal (decontam)
#
# Authors: Megan Stroup Hill, Claude
# Description: Identify and remove reagent contaminant ASVs with decontam
#              (prevalence method) against the 12 negative-control libraries
# Input: result/biom-with-taxonomy.txt (raw ASV counts from 01_sequence_processing.sh)
#        data/metadata_run1_10.csv
# Output: result/decontam/
#   - biom-with-taxonomy_decontam.txt   sink samples only, contaminant ASVs removed
#   - decontam_scores_all_ASVs.tsv      every ASV: score, control prevalence, reads
#   - removed_ASVs.tsv                  ASVs meeting the removal rule
#   - per_sample_reads_removed.tsv      reads lost per sample
#   - decontam_prevalence_plot.pdf/.png
#
# Removal rule: decontam prevalence score P < 0.1 AND present in >= 50% of the
# negative controls of every control site type (TP 4, PT 3, TOP 5). This keeps
# ASVs seen in only one or two blanks, which is more likely carryover from sink
# samples than reagent contamination. Two ASVs are removed (Bradyrhizobium
# 87bd8926, Qipengyuania 94df785b; 1.4% of sink-sample reads).
#
# decontam's frequency/combined modes were tested and not used: treatment lowers
# DNA yield, so the frequency test flags taxa that are genuinely enriched in
# low-biomass treated samples (see manuscript Methods).
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(decontam)
  library(ggplot2)
})

cat("Starting decontam contaminant identification...\n")

in_fp   <- "result/biom-with-taxonomy.txt"
meta_fp <- "data/metadata_run1_10.csv"
out_dir <- "result/decontam"
threshold <- 0.1
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# Step 1: Load feature table and metadata
# =============================================================================
tab <- fread(in_fp, skip = 1, sep = "\t", header = TRUE)
setnames(tab, 1, "ASV")
taxonomy <- tab[, .(ASV, taxonomy)]
counts <- as.matrix(tab[, !c("ASV", "taxonomy")])
rownames(counts) <- tab$ASV
rm(tab); invisible(gc())

# Batch 9 columns carry "_corrected" in the feature table but not in the
# metadata; without this, all 229 Batch 9 sink samples are silently dropped.
colnames(counts) <- sub("_FastQ_corrected_", "_FastQ_", colnames(counts))

meta <- fread(meta_fp, na.strings = c("NA", ""))
meta[, is_neg := Room %in% c("C-", "C", "NegCon")]
meta[, is_true := grepl("^[0-9]+$", Room) & Sample_Type %in% c("TP", "PT", "TOP")]
meta <- meta[SeqID %in% colnames(counts) & (is_neg | is_true)]

counts <- counts[, meta$SeqID]
counts <- counts[rowSums(counts) > 0, ]
cat("✓ Samples:", sum(meta$is_true), "sink samples,", sum(meta$is_neg), "negative controls\n")
cat("✓ ASVs with reads in these samples:", nrow(counts), "\n")

# =============================================================================
# Step 2: decontam prevalence scores
# =============================================================================
res <- isContaminant(t(counts), method = "prevalence", neg = meta$is_neg,
                     threshold = threshold)
res <- as.data.table(res, keep.rownames = "ASV")

neg_pa  <- counts[, meta$is_neg]  > 0
true_pa <- counts[, meta$is_true] > 0
res[, `:=`(
  n_neg_present  = rowSums(neg_pa)[ASV],
  n_true_present = rowSums(true_pa)[ASV],
  prev_neg  = round(rowSums(neg_pa)[ASV]  / sum(meta$is_neg), 4),
  prev_true = round(rowSums(true_pa)[ASV] / sum(meta$is_true), 4),
  reads_neg  = rowSums(counts[, meta$is_neg])[ASV],
  reads_true = rowSums(counts[, meta$is_true])[ASV]
)]
res[, contaminant_t0.1 := !is.na(p) & p < threshold]

# Control prevalence by the site type each control was collected with
ctrl_types <- sort(unique(meta[is_neg == TRUE, Sample_Type]))
for (st in ctrl_types) {
  ids <- meta[is_neg == TRUE & Sample_Type == st, SeqID]
  res[, (paste0("prev_neg_", st)) :=
        round(rowSums(counts[ASV, ids, drop = FALSE] > 0) / length(ids), 4)]
}
prev_cols <- paste0("prev_neg_", ctrl_types)
res[, remove := contaminant_t0.1 & Reduce(`&`, lapply(.SD, `>=`, 0.5)), .SDcols = prev_cols]
res <- merge(res, taxonomy, by = "ASV", all.x = TRUE)
setorder(res, p, na.last = TRUE)
fwrite(res, file.path(out_dir, "decontam_scores_all_ASVs.tsv"), sep = "\t")
fwrite(res[remove == TRUE], file.path(out_dir, "removed_ASVs.tsv"), sep = "\t")

cat(sprintf("✓ Flagged at P < %.1f: %d ASVs; removed under the control-prevalence rule: %d ASVs (%.2f%% of sink-sample reads)\n",
            threshold, sum(res$contaminant_t0.1), sum(res$remove),
            100 * sum(res[remove == TRUE, reads_true]) / sum(res$reads_true)))

# =============================================================================
# Step 3: Per-sample impact and decontaminated table
# =============================================================================
contam_ids <- res[remove == TRUE, ASV]
lib  <- colSums(counts)
lost <- colSums(counts[contam_ids, , drop = FALSE])
impact <- data.table(meta[, .(SeqID, Sample, Batch, Room, Week, Sample_Type,
                              Virasept_Treated, is_neg)],
                     reads_before = lib, contaminant_reads = lost,
                     pct_removed = round(100 * lost / pmax(lib, 1), 2))
fwrite(impact, file.path(out_dir, "per_sample_reads_removed.tsv"), sep = "\t")

# Sink samples only, contaminant ASVs removed (same layout as biom-with-taxonomy.txt)
clean <- counts[setdiff(rownames(counts), contam_ids), meta$is_true]
clean <- clean[rowSums(clean) > 0, ]
out <- data.table(`#OTU ID` = rownames(clean), clean)
out[, taxonomy := taxonomy$taxonomy[match(`#OTU ID`, taxonomy$ASV)]]
clean_fp <- file.path(out_dir, "biom-with-taxonomy_decontam.txt")
writeLines("# Constructed from biom file (decontam prevalence P < 0.1 and >= 50% prevalence in every control site type)", clean_fp)
fwrite(out, clean_fp, sep = "\t", append = TRUE, col.names = TRUE)
cat("✓ Decontaminated table:", nrow(clean), "ASVs ×", ncol(clean), "sink samples\n")

# =============================================================================
# Step 4: Prevalence plot
# =============================================================================
pd <- res[n_neg_present > 0 | n_true_present > 0]
pd[, status := fifelse(remove, "Removed (P < 0.1, >= 50% of each control type)",
                fifelse(contaminant_t0.1, "Flagged P < 0.1, kept", "Not flagged"))]
p <- ggplot(pd, aes(prev_true, prev_neg, colour = status)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_colour_manual(values = c("Not flagged" = "grey60", "Flagged P < 0.1, kept" = "#e1a43a",
                                 "Removed (P < 0.1, >= 50% of each control type)" = "#c0392b"), name = NULL) +
  guides(colour = guide_legend(nrow = 2)) +
  labs(x = "Prevalence in sink samples", y = "Prevalence in negative controls") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(out_dir, "decontam_prevalence_plot.pdf"), p, width = 6.5, height = 5.5)
ggsave(file.path(out_dir, "decontam_prevalence_plot.png"), p, width = 6.5, height = 5.5, dpi = 300)

cat("\nRemoved ASVs:\n")
print(res[remove == TRUE, c("ASV", "p", prev_cols, "reads_neg", "reads_true"), with = FALSE])
cat("\n✓ decontam completed\n")
