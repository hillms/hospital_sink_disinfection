# =============================================================================
# Hospital Sink Disinfection Study - Data Preprocessing
#
# Authors: Megan Stroup Hill, Claude Sonnet 4
# Description: Load and preprocess data for downstream analyses
# Input: parsed_taxa_table.csv, genus_counts_raw.csv, metadata_run1_10.csv
#        (built from the decontaminated table by preprocessing/02-03)
# Output: Cleaned community matrices and metadata
#   comm      depth-normalized genus abundances (samples × genera)
#   comm_raw  raw genus read counts, same samples (rarefied alpha diversity, ALDEx2)
# =============================================================================

source("analysis/00_load_libraries.R")

cat("Starting data preprocessing...\n")
cat("Date:", as.character(Sys.Date()), "\n\n")

# =============================================================================
# Step 1: Load Data Files
# =============================================================================
cat("Step 1: Loading data files...\n")

# Load genus-level abundance table
if (!file.exists("data/parsed_taxa_table.csv")) {
  stop("Error: parsed_taxa_table.csv not found in data/ directory")
}

taxa <- read.csv("data/parsed_taxa_table.csv", header = TRUE, row.names = 1, check.names = FALSE)
cat("✓ Loaded taxa table:", nrow(taxa), "genera ×", ncol(taxa), "samples\n")

# Load metadata
if (!file.exists("data/metadata_run1_10.csv")) {
  stop("Error: metadata_run1_10.csv not found in data/ directory")
}

meta <- read.csv("data/metadata_run1_10.csv", header = TRUE, check.names = FALSE)
cat("✓ Loaded metadata:", nrow(meta), "samples with", ncol(meta), "variables\n")

# =============================================================================
# Step 2: Standard Genus-Level Filtering
# =============================================================================
cat("\nStep 2: Applying genus-level filtering...\n")

# Apply standard filtering for genus-level taxa
all_rows <- rownames(taxa)
genus_rows <- grepl("g__[A-Za-z]", all_rows) &
  !grepl("s__", all_rows) &
  !grepl("g__[^-]+--__--[0-9a-f]", all_rows) &   # exclude ASV hashes
  !grepl("g__[^-]+--__$", all_rows)              # exclude genus--__ entries

genus_table <- taxa[genus_rows, ]

# Clean genus names
genus_table$genus <- gsub(".*g__", "", rownames(genus_table))
genus_table$genus <- gsub("--.*", "", genus_table$genus)

# Aggregate duplicates (same genus from different ASVs)
genus_table <- genus_table %>%
  group_by(genus) %>%
  summarise(across(everything(), sum)) %>%
  as.data.frame()

rownames(genus_table) <- genus_table$genus
genus_table$genus <- NULL

cat("✓ Genus filtering applied:", nrow(genus_table), "unique genera retained\n")

# =============================================================================
# Step 3: Metadata Cleaning and Filtering
# =============================================================================
cat("\nStep 3: Cleaning metadata...\n")

# Clean metadata following established pipeline
meta_clean <- meta %>%
  # Remove control samples
  filter(!Room %in% c("C-", "C", "NegCon")) %>%
  # Remove samples with missing treatment information
  filter(!is.na(Virasept_Treated)) %>%
  # Convert week to numeric and filter invalid entries
  mutate(Week_num = suppressWarnings(as.numeric(Week))) %>%
  filter(!is.na(Week_num)) %>%
  # Keep only relevant sample types
  filter(Sample_Type %in% c("TP", "PT", "TOP"))

cat("✓ Metadata filtering applied\n")
cat("   Samples removed: Controls, missing treatment, invalid weeks\n")
cat("   Remaining samples:", nrow(meta_clean), "\n")

# =============================================================================
# Step 4: Match Samples Between Taxa and Metadata
# =============================================================================
cat("\nStep 4: Matching samples between data files...\n")

# Find shared samples
shared_samples <- intersect(meta_clean$SeqID, colnames(genus_table))
meta_clean <- meta_clean %>% filter(SeqID %in% shared_samples)

cat("✓ Shared samples identified:", length(shared_samples), "\n")

# Create standardized treatment labels
meta_clean$Treatment <- ifelse(meta_clean$Virasept_Treated == "Yes", "Treated", "Untreated")

# =============================================================================
# Step 5: Create Community Matrix
# =============================================================================
cat("\nStep 5: Creating community matrix...\n")

# Transpose and subset genus table to match metadata
comm <- t(genus_table[, meta_clean$SeqID])

# Remove zero-sum samples (if any)
keep <- rowSums(comm) > 0
comm <- comm[keep, ]
meta_clean <- meta_clean[keep, ]

cat("✓ Community matrix created:", nrow(comm), "samples ×", ncol(comm), "genera\n")
cat("   Zero-sum samples removed:", sum(!keep), "\n")

# Raw genus read counts for the same samples (all genera, alphabetical)
if (!file.exists("data/genus_counts_raw.csv")) {
  stop("Error: genus_counts_raw.csv not found in data/ directory (run preprocessing/03_build_taxa_tables.R)")
}
genus_raw <- read.csv("data/genus_counts_raw.csv", row.names = 1, check.names = FALSE)
comm_raw <- t(as.matrix(genus_raw[, meta_clean$SeqID]))
rm(genus_raw)
cat("✓ Raw count matrix created:", nrow(comm_raw), "samples ×", ncol(comm_raw), "genera\n")
cat("   Median reads per sample:", median(rowSums(comm_raw)), "\n")

# =============================================================================
# Step 6: Calculate Relative Abundances
# =============================================================================
cat("\nStep 6: Calculating relative abundances...\n")

# Calculate relative abundances
comm_rel <- comm / rowSums(comm)

cat("✓ Relative abundances calculated\n")

# =============================================================================
# Step 7: Apply Final Factor Ordering
# =============================================================================
cat("\nStep 7: Applying factor ordering...\n")

# Apply consistent factor ordering (Tail Piece → P-Trap → Countertop)
meta_clean$Sample_Type <- factor(meta_clean$Sample_Type,
                                levels = SAMPLE_TYPE_LEVELS,
                                labels = SAMPLE_TYPE_LABELS)

meta_clean$Treatment <- factor(meta_clean$Treatment,
                              levels = c("Treated", "Untreated"))

cat("✓ Factor ordering applied\n")
cat("   Sample type levels:", paste(levels(meta_clean$Sample_Type), collapse = " → "), "\n")
cat("   Treatment levels:", paste(levels(meta_clean$Treatment), collapse = ", "), "\n")

# =============================================================================
# Step 8: Generate Summary Statistics
# =============================================================================
cat("\nStep 8: Generating summary statistics...\n")

# Sample distribution summary
sample_summary <- meta_clean %>%
  group_by(Sample_Type, Treatment) %>%
  summarise(n_samples = n(), .groups = "drop") %>%
  pivot_wider(names_from = Treatment, values_from = n_samples, values_fill = 0)

print(sample_summary)

# Week distribution
week_summary <- meta_clean %>%
  summarise(
    total_weeks = length(unique(Week_num)),
    week_range = paste(range(Week_num), collapse = "-"),
    missing_weeks = paste(setdiff(1:35, unique(Week_num)), collapse = ", ")
  )

cat("\n✓ Data covers", week_summary$total_weeks, "weeks (", week_summary$week_range, ")\n")
if (week_summary$missing_weeks != "") {
  cat("   Missing weeks:", week_summary$missing_weeks, "\n")
}

# =============================================================================
# Step 9: Save Processed Data
# =============================================================================
cat("\nStep 9: Saving processed data...\n")

# Save main objects for downstream analysis
save(comm, comm_rel, comm_raw, meta_clean, genus_table, sample_summary,
     file = "analysis/processed_data.RData")

# Also save as individual CSV files
write.csv(comm_rel, "analysis/community_relative_abundances.csv")
write.csv(meta_clean, "analysis/metadata_clean.csv", row.names = FALSE)

cat("✓ Processed data saved:\n")
cat("   - analysis/processed_data.RData (all objects)\n")
cat("   - analysis/community_relative_abundances.csv\n")
cat("   - analysis/metadata_clean.csv\n")

# =============================================================================
# Step 10: Data Quality Checks
# =============================================================================
cat("\nStep 10: Final quality checks...\n")

# Check for any remaining issues
issues <- c()

# Check for zero-variance genera
zero_var_genera <- colnames(comm_rel)[apply(comm_rel, 2, var) == 0]
if (length(zero_var_genera) > 0) {
  issues <- c(issues, paste("Zero-variance genera:", length(zero_var_genera)))
}

# Check sample size balance
min_group_size <- min(sample_summary[, -1])
if (min_group_size < 10) {
  issues <- c(issues, paste("Small group size detected:", min_group_size))
}

# Check for extreme library sizes
lib_sizes <- rowSums(comm)
cv_lib_size <- sd(lib_sizes) / mean(lib_sizes)
if (cv_lib_size > 2) {
  issues <- c(issues, paste("High library size variation: CV =", round(cv_lib_size, 2)))
}

if (length(issues) == 0) {
  cat("✓ All quality checks passed\n")
} else {
  cat("⚠ Quality issues detected:\n")
  for (issue in issues) {
    cat("  -", issue, "\n")
  }
}

cat("\n=== DATA PREPROCESSING COMPLETED ===\n")
cat("Final dataset:\n")
cat("  Samples:", nrow(comm_rel), "\n")
cat("  Genera:", ncol(comm_rel), "\n")
cat("  Sample types:", length(levels(meta_clean$Sample_Type)), "\n")
cat("  Treatments:", length(levels(meta_clean$Treatment)), "\n")
cat("  Study duration:", week_summary$total_weeks, "weeks\n")
cat("\nReady for downstream analyses!\n")