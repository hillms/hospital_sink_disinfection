# =============================================================================
# Hospital Sink Disinfection Study - Data Normalization
#
# Author: Shan Sun
# Description: Normalize feature table using proportional scaling
# Input: biom-with-taxonomy.txt from QIIME2 processing
# Output: parsed_taxa_table.csv (normalized genus-level table)
# =============================================================================

# Load required libraries
library(dplyr)
library(readr)

cat("Starting data normalization pipeline...\n")
cat("Author: Shan Sun\n")
cat("Date:", as.character(Sys.Date()), "\n\n")

# =============================================================================
# Step 1: Load and parse BIOM output
# =============================================================================
cat("Step 1: Loading BIOM table...\n")

# Read the BIOM table (skip first line which is a comment)
biom_data <- read_tsv("result/biom-with-taxonomy.txt", skip = 1)

# Extract feature table (remove taxonomy column)
feature_table <- biom_data %>%
  select(-taxonomy) %>%
  column_to_rownames("# Constructed from biom file")

cat("✓ Loaded", nrow(feature_table), "features across", ncol(feature_table), "samples\n")

# =============================================================================
# Step 2: Apply proportional normalization
# =============================================================================
cat("Step 2: Applying proportional normalization...\n")

# Shan's normalization: scale by mean library size
# Formula: tab_all = t(t(tab_all)/colSums(tab_all))*mean(colSums(tab_all))

# Convert to matrix for calculations
tab_all <- as.matrix(feature_table)

# Calculate library sizes (column sums)
lib_sizes <- colSums(tab_all)
mean_lib_size <- mean(lib_sizes)

# Apply proportional scaling
tab_normalized <- t(t(tab_all) / lib_sizes) * mean_lib_size

cat("✓ Normalization applied\n")
cat("   Mean library size:", round(mean_lib_size, 0), "\n")
cat("   Library size range:", round(min(lib_sizes), 0), "-", round(max(lib_sizes), 0), "\n")

# =============================================================================
# Step 3: Parse taxonomy and aggregate to genus level
# =============================================================================
cat("Step 3: Parsing taxonomy and aggregating to genus level...\n")

# Extract taxonomy information
taxonomy_info <- biom_data %>%
  select(`# Constructed from biom file`, taxonomy)

# Parse taxonomy string and extract genus
taxonomy_parsed <- taxonomy_info %>%
  mutate(
    # Split taxonomy string and extract components
    tax_split = strsplit(taxonomy, ";"),
    # Extract genus (format: "g__GenusName")
    genus_raw = sapply(tax_split, function(x) {
      genus_level <- grep("g__", x, value = TRUE)
      if (length(genus_level) > 0) return(genus_level[1]) else return(NA)
    })
  ) %>%
  filter(!is.na(genus_raw)) %>%
  mutate(
    # Clean genus names
    genus_clean = gsub("g__", "", genus_raw),
    genus_clean = gsub("--.*", "", genus_clean)  # Remove any suffix after --
  ) %>%
  filter(
    # Apply standard genus filtering (from improved_plot.R)
    grepl("^[A-Za-z]", genus_clean),  # Must start with letter
    !grepl("__", genus_raw),          # No double underscores
    genus_clean != ""                 # Not empty
  )

# Create genus-level abundance table
genus_table <- tab_normalized[taxonomy_parsed$`# Constructed from biom file`, ] %>%
  as.data.frame() %>%
  mutate(genus = taxonomy_parsed$genus_clean) %>%
  group_by(genus) %>%
  summarise(across(everything(), sum), .groups = "drop") %>%
  column_to_rownames("genus")

cat("✓ Aggregated to", nrow(genus_table), "genera\n")

# =============================================================================
# Step 4: Final filtering and quality control
# =============================================================================
cat("Step 4: Quality control filtering...\n")

# Remove genera with very low abundance across all samples
min_total_abundance <- 0.001 * ncol(genus_table)  # At least 0.001 average across samples
genus_table_filtered <- genus_table[rowSums(genus_table) > min_total_abundance, ]

# Remove samples with zero reads (if any)
sample_totals <- colSums(genus_table_filtered)
genus_table_final <- genus_table_filtered[, sample_totals > 0]

cat("✓ Final table:", nrow(genus_table_final), "genera ×", ncol(genus_table_final), "samples\n")

# =============================================================================
# Step 5: Save processed data
# =============================================================================
cat("Step 5: Saving processed data...\n")

# Save genus-level table
write.csv(genus_table_final, "parsed_taxa_table.csv", row.names = TRUE)

# Save processing metadata
processing_info <- data.frame(
  step = c("raw_features", "after_normalization", "after_genus_aggregation", "final_filtered"),
  n_features = c(nrow(feature_table), nrow(tab_normalized), nrow(genus_table), nrow(genus_table_final)),
  n_samples = c(ncol(feature_table), ncol(tab_normalized), ncol(genus_table), ncol(genus_table_final)),
  timestamp = Sys.time()
)

write.csv(processing_info, "normalization_summary.csv", row.names = FALSE)

cat("✓ Files saved:\n")
cat("   - parsed_taxa_table.csv (genus-level normalized abundances)\n")
cat("   - normalization_summary.csv (processing metadata)\n")

# =============================================================================
# Step 6: Generate summary statistics
# =============================================================================
cat("\nStep 6: Summary statistics...\n")

cat("=== NORMALIZATION SUMMARY ===\n")
cat("Input features:", nrow(feature_table), "\n")
cat("Final genera:", nrow(genus_table_final), "\n")
cat("Final samples:", ncol(genus_table_final), "\n")
cat("Mean library size:", round(mean_lib_size, 0), "\n")
cat("Normalization method: Proportional scaling to mean library size\n")
cat("Taxonomic level: Genus (GTDB)\n")
cat("Date processed:", as.character(Sys.Date()), "\n")

cat("\n✓ Normalization pipeline completed successfully!\n")
cat("Next step: Load parsed_taxa_table.csv for downstream analyses\n")