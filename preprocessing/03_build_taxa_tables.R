# =============================================================================
# Hospital Sink Disinfection Study - Build Taxa Tables
#
# Authors: Shan Sun (normalization), Megan Stroup Hill, Claude
# Description: Build the two tables used downstream from the decontaminated
#              feature table:
#   1. data/parsed_taxa_table.csv - depth-normalized table in the parsed
#      taxonomy format (volatility, instability, beta diversity, top taxa)
#   2. data/genus_counts_raw.csv  - raw read counts per GTDB genus
#      (rarefied alpha diversity and ALDEx2 differential abundance)
# Input: result/decontam/biom-with-taxonomy_decontam.txt (02_decontam.R)
# =============================================================================

suppressPackageStartupMessages(library(data.table))

cat("Starting taxa table construction...\n")

in_fp <- "result/decontam/biom-with-taxonomy_decontam.txt"
if (!file.exists(in_fp)) stop("Run preprocessing/02_decontam.R first")

tab <- fread(in_fp, skip = 1, sep = "\t", header = TRUE)
setnames(tab, 1, "ASV")
tax <- tab$taxonomy
counts <- as.matrix(tab[, !c("ASV", "taxonomy")])
asv <- tab$ASV
rm(tab); invisible(gc())
cat("✓ Loaded", nrow(counts), "ASVs ×", ncol(counts), "samples\n")

# =============================================================================
# Step 1: Raw genus-level counts
# =============================================================================
# Genus = GTDB g__ name (spaces -> "_"); ASVs without a genus are excluded
genus <- ifelse(grepl("g__[A-Za-z]", tax), gsub(" ", "_", sub(".*g__([^;]+).*", "\\1", tax)), NA)
keep <- !is.na(genus)
genus_raw <- rowsum(counts[keep, ], genus[keep])
fwrite(data.table(genus = rownames(genus_raw), genus_raw), "data/genus_counts_raw.csv")
cat("✓ data/genus_counts_raw.csv:", nrow(genus_raw), "genera\n")
rm(genus_raw); invisible(gc())

# =============================================================================
# Step 2: Depth-normalized parsed taxa table
# =============================================================================
# Shan's normalization: scale every sample to the mean library size
# tab_all = t(t(tab_all)/colSums(tab_all))*mean(colSums(tab_all))
counts <- t(t(counts) / colSums(counts)) * mean(colSums(counts))
counts[!is.finite(counts)] <- 0

# Taxonomy padded to 7 ranks with "__", spaces -> "_", ASV hash appended; one
# row per taxonomy prefix (phylum ... genus--species--hash), joined with "--".
# analysis/01_data_preprocessing.R keeps the genus-level rows.
ranks <- lapply(strsplit(tax, ";\\s*"), function(x) {
  x <- gsub(" ", "_", trimws(x))
  c(x, rep("__", 7 - length(x)))[1:7]
})
lvl <- cbind(do.call(rbind, ranks), asv)

out <- vector("list", 8)
for (k in 2:8) {
  key <- do.call(paste, c(as.data.frame(lvl[, 1:k]), sep = "--"))
  out[[k]] <- rowsum(counts, key, reorder = FALSE)
}
res <- do.call(rbind, out)
rm(out, counts); invisible(gc())

dt <- data.table(V1 = rownames(res), res)
setnames(dt, "V1", "")
fwrite(dt, "data/parsed_taxa_table.csv")
cat("✓ data/parsed_taxa_table.csv:", nrow(res), "rows ×", ncol(res), "samples\n")
cat("\n✓ Taxa tables built\n")
