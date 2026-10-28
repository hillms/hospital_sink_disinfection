# =============================================================================
# Hospital Sink Disinfection Study - Differential Abundance Calibration
#
# Authors: Megan Stroup Hill, Claude
# Description: Check how many genera each differential abundance method calls
#              when there is no true treatment effect (Table S2), and whether
#              the tail piece / P-trap agreement exceeds what room structure
#              alone produces.
# Input: processed_data.RData, differential_abundance.RData
# Output: analysis/da_calibration_true_labels.csv
#         analysis/da_calibration_shuffles.csv
#         analysis/da_concordance_shuffles.csv
#
# Treatment was assigned by room and each room was sampled ~35 times, so
# samples are not independent. Shuffling the treatment labels among the 31
# rooms (a room's TP, PT and TOP samples keep the same label) removes any real
# treatment effect while keeping room-level structure; a well-calibrated method
# should then call ~0 genera.
#   - ALDEx2 (04_differential_abundance.R settings)
#   - DESeq2: negative binomial Wald test, poscounts size factors
#   - edgeR: quasi-likelihood F-test, TMMwsp normalization, robust dispersion
# Same samples, genera and >= 5% prevalence filter for every method; shuffles
# use seeds 101-120 (1-5 for all methods and locations, 1-20 for ALDEx2 TP/PT).
# Result: only ALDEx2 is close to calibrated (0-6 genera with shuffled labels),
# so the manuscript reports ALDEx2 only. Runtime is roughly 10-15 minutes.
# =============================================================================

source("analysis/00_load_libraries.R")
suppressPackageStartupMessages({
  library(DESeq2)
  library(edgeR)
})

cat("Starting differential abundance calibration...\n")

load("analysis/processed_data.RData")
load("analysis/differential_abundance.RData")

site_counts <- function(st) {
  idx <- meta_clean$Sample_Type == SAMPLE_TYPE_LABELS[SAMPLE_TYPE_LEVELS == st]
  sc <- t(comm_raw[idx, ])
  list(idx = idx, counts = sc[rowSums(sc > 0) / ncol(sc) >= 0.05, ])
}

# Number of significant genera (BH/FDR < 0.05) for one method, location and labelling
n_sig <- function(method, st, treatment) {
  s <- site_counts(st)
  trt <- factor(as.character(treatment[s$idx]), levels = c("Treated", "Untreated"))
  if (method == "ALDEx2") {
    r <- run_aldex(st, treatment)
    return(sum(r$wi.eBH < 0.05 | r$we.eBH < 0.05))
  }
  if (method == "DESeq2") {
    dds <- DESeqDataSetFromMatrix(round(s$counts),
                                  data.frame(Treatment = trt, row.names = colnames(s$counts)),
                                  ~ Treatment)
    dds <- estimateSizeFactors(dds, type = "poscounts")
    dds <- tryCatch(DESeq(dds, quiet = TRUE),
                    error = function(e) DESeq(dds, fitType = "local", quiet = TRUE))
    padj <- results(dds, contrast = c("Treatment", "Untreated", "Treated"), alpha = 0.05)$padj
    return(sum(padj < 0.05, na.rm = TRUE))
  }
  y <- calcNormFactors(DGEList(s$counts), method = "TMMwsp")
  design <- model.matrix(~ trt)
  y <- estimateDisp(y, design)
  fdr <- topTags(glmQLFTest(glmQLFit(y, design, robust = TRUE), coef = 2), n = Inf)$table$FDR
  sum(fdr < 0.05)
}

# Room-level shuffle i: treatment labels permuted among rooms
rooms <- meta_clean %>% distinct(Room, Treatment) %>% filter(!duplicated(Room))
shuffle <- function(i) {
  set.seed(100 + i)
  perm <- setNames(as.character(sample(rooms$Treatment)), rooms$Room)
  perm[meta_clean$Room]
}

methods <- c("ALDEx2", "DESeq2", "edgeR")

# =============================================================================
# Step 1: True labels
# =============================================================================
cat("Step 1: Significant genera with the true treatment labels...\n")
true_labels <- expand.grid(Method = methods, Sample_Type = SAMPLE_TYPE_LEVELS,
                           stringsAsFactors = FALSE)
true_labels$n_tested <- sapply(true_labels$Sample_Type, function(st) nrow(site_counts(st)$counts))
true_labels$n_sig <- mapply(function(m, st) {
  if (m == "ALDEx2") return(sum(aldex_results$significant[aldex_results$Sample_Type == st]))
  n_sig(m, st, meta_clean$Treatment)
}, true_labels$Method, true_labels$Sample_Type)
print(true_labels)

# =============================================================================
# Step 2: Five room-level shuffles, all methods and locations (Table S2)
# =============================================================================
cat("\nStep 2: Room-level label shuffles (5 per method and location)...\n")
shuffles <- list()
for (i in 1:5) {
  trt <- shuffle(i)
  for (m in methods) for (st in SAMPLE_TYPE_LEVELS) {
    shuffles[[length(shuffles) + 1]] <- data.frame(shuffle = i, Method = m, Sample_Type = st,
                                                    n_sig = n_sig(m, st, trt))
    cat(sprintf("  shuffle %d %s %s: %d\n", i, m, st, tail(shuffles, 1)[[1]]$n_sig))
  }
}
shuffles <- bind_rows(shuffles)

# =============================================================================
# Step 3: Tail piece / P-trap agreement, ALDEx2, 20 shuffles
# =============================================================================
# Genera significant at both plumbing sites with the same direction of effect.
# (The TP-PT effect-size correlation is not used: shuffled labels give r values
# as high as the true labels because both sites share room-level structure.)
cat("\nStep 3: Shared tail piece / P-trap genera with shuffled labels (20 shuffles)...\n")
both_same_dir <- function(trt) {
  r <- bind_rows(run_aldex("TP", trt), run_aldex("PT", trt)) %>%
    mutate(sig = wi.eBH < 0.05 | we.eBH < 0.05) %>%
    select(genus, Sample_Type, effect, sig) %>%
    pivot_wider(names_from = Sample_Type, values_from = c(effect, sig))
  c(sig_TP = sum(r$sig_TP, na.rm = TRUE), sig_PT = sum(r$sig_PT, na.rm = TRUE),
    shared = sum(r$sig_TP %in% TRUE & r$sig_PT %in% TRUE & sign(r$effect_TP) == sign(r$effect_PT)))
}
concordance <- bind_rows(
  data.frame(run = "true labels", t(both_same_dir(meta_clean$Treatment))),
  bind_rows(lapply(1:20, function(i) {
    x <- data.frame(run = paste("shuffle", i), t(both_same_dir(shuffle(i))))
    print(x)
    x
  })))

# =============================================================================
# Save Results
# =============================================================================
write.csv(true_labels, "analysis/da_calibration_true_labels.csv", row.names = FALSE)
write.csv(shuffles, "analysis/da_calibration_shuffles.csv", row.names = FALSE)
write.csv(concordance, "analysis/da_concordance_shuffles.csv", row.names = FALSE)

cat("\n=== CALIBRATION SUMMARY (true labels / shuffled range) ===\n")
summ <- shuffles %>%
  group_by(Method, Sample_Type) %>%
  summarise(shuffled = sprintf("%d-%d", min(n_sig), max(n_sig)), .groups = "drop") %>%
  left_join(true_labels, by = c("Method", "Sample_Type"))
print(summ)
sh <- concordance[concordance$run != "true labels", ]
cat(sprintf("\nShared TP/PT genera: true %d vs shuffled median %g (range %d-%d)\n",
            concordance$shared[1], median(sh$shared), min(sh$shared), max(sh$shared)))
cat(sprintf("ALDEx2 significant genera with shuffled labels: TP %d-%d, PT %d-%d\n",
            min(sh$sig_TP), max(sh$sig_TP), min(sh$sig_PT), max(sh$sig_PT)))

cat("\n✓ Calibration analysis completed!\n")
