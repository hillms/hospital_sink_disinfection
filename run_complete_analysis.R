# =============================================================================
# Hospital Sink Disinfection Study - Complete Analysis Pipeline
#
# Authors: Megan Stroup Hill, Claude Sonnet 4, Claude
# Description: Master script to run the complete analysis pipeline
# Requires data/ tables built by preprocessing/02_decontam.R and
# preprocessing/03_build_taxa_tables.R.
# Set RUN_CALIBRATION <- FALSE to skip 04b_da_calibration.R (~1-2 hours).
# =============================================================================

if (!exists("RUN_CALIBRATION")) RUN_CALIBRATION <- TRUE

# Set working directory to repository root
if (!file.exists("run_complete_analysis.R")) {
  stop("Please run this script from the repository root directory")
}

cat("=== HOSPITAL SINK DISINFECTION ANALYSIS PIPELINE ===\n")
cat("Starting complete analysis...\n")
cat("Timestamp:", as.character(Sys.time()), "\n\n")

# =============================================================================
# Step 0: Check Prerequisites
# =============================================================================
cat("Step 0: Checking prerequisites...\n")

# Check for required data files
required_files <- c(
  "data/parsed_taxa_table.csv",
  "data/genus_counts_raw.csv",
  "data/metadata_run1_10.csv"
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required data files:\n", paste(missing_files, collapse = "\n"))
}

# Test library loading
tryCatch({
  source("analysis/00_load_libraries.R")
  cat("✓ All required libraries available\n")
}, error = function(e) {
  stop("Error loading required libraries. Please install packages from requirements.txt\n",
       "Error: ", e$message)
})

# =============================================================================
# Step 1: Data Preprocessing
# =============================================================================
cat("\nStep 1: Running data preprocessing...\n")

tryCatch({
  source("analysis/01_data_preprocessing.R")
  cat("✓ Data preprocessing completed\n")
}, error = function(e) {
  stop("Error in data preprocessing: ", e$message)
})

# =============================================================================
# Step 2: Alpha Diversity Analysis
# =============================================================================
cat("\nStep 2: Running alpha diversity analysis...\n")

tryCatch({
  source("analysis/02_alpha_diversity.R")
  cat("✓ Alpha diversity analysis completed\n")
}, error = function(e) {
  stop("Error in alpha diversity analysis: ", e$message)
})

# =============================================================================
# Step 3: Beta Diversity Analysis (if available)
# =============================================================================
cat("\nStep 3: Running beta diversity analysis...\n")

if (file.exists("analysis/03_beta_diversity.R")) {
  tryCatch({
    source("analysis/03_beta_diversity.R")
    cat("✓ Beta diversity analysis completed\n")
  }, error = function(e) {
    cat("⚠ Warning: Beta diversity analysis failed: ", e$message, "\n")
  })
} else {
  cat("ℹ Beta diversity script not found, skipping...\n")
}

# =============================================================================
# Step 4: Differential Abundance Analysis
# =============================================================================
cat("\nStep 4: Running differential abundance analysis (ALDEx2, raw counts)...\n")

tryCatch({
  source("analysis/04_differential_abundance.R")
  cat("✓ Differential abundance analysis completed\n")
}, error = function(e) {
  stop("Error in differential abundance analysis: ", e$message)
})

if (RUN_CALIBRATION) {
  cat("\nStep 4b: Running differential abundance calibration (room-shuffled labels)...\n")
  tryCatch({
    source("analysis/04b_da_calibration.R")
    cat("✓ Calibration analysis completed\n")
  }, error = function(e) {
    cat("⚠ Warning: Calibration analysis failed: ", e$message, "\n")
  })
} else {
  cat("\nStep 4b: Calibration skipped (RUN_CALIBRATION = FALSE)\n")
}

# =============================================================================
# Step 5: Community Volatility Analysis
# =============================================================================
cat("\nStep 5: Running community volatility analysis...\n")

tryCatch({
  source("analysis/05_community_volatility.R")
  cat("✓ Community volatility analysis completed\n")
}, error = function(e) {
  stop("Error in community volatility analysis: ", e$message)
})

# =============================================================================
# Step 6: Community Instability Analysis
# =============================================================================
cat("\nStep 6: Running community instability analysis...\n")

tryCatch({
  source("analysis/06_community_instability.R")
  cat("✓ Community instability analysis completed\n")
}, error = function(e) {
  stop("Error in community instability analysis: ", e$message)
})

# =============================================================================
# Step 6b: Weekly Treatment PERMANOVA (Table S3)
# =============================================================================
cat("\nStep 6b: Running weekly PERMANOVA...\n")

tryCatch({
  source("analysis/07_weekly_permanova.R")
  cat("✓ Weekly PERMANOVA completed\n")
}, error = function(e) {
  stop("Error in weekly PERMANOVA: ", e$message)
})

# =============================================================================
# Step 7: Generate Figures
# =============================================================================
cat("\nStep 7: Generating manuscript figures...\n")

figure_scripts <- list.files("figures", pattern = "figure_.*\\.R$", full.names = TRUE)

for (script in figure_scripts) {
  figure_name <- basename(script)
  cat("  Generating", figure_name, "...\n")

  tryCatch({
    source(script)
    cat("  ✓", figure_name, "completed\n")
  }, error = function(e) {
    cat("  ⚠ Warning:", figure_name, "failed:", e$message, "\n")
  })
}

# Tables 1 and S1-S3 (Table S2 needs the calibration outputs)
cat("\nGenerating manuscript tables...\n")
tryCatch({
  source("tables/make_tables.R")
}, error = function(e) {
  cat("⚠ Warning: table generation failed:", e$message, "\n")
})

# =============================================================================
# Step 8: Generate Analysis Summary
# =============================================================================
cat("\nStep 8: Generating analysis summary...\n")

# Collect key results from various analyses
summary_data <- list()

# Alpha diversity results
if (file.exists("analysis/alpha_diversity_statistics.csv")) {
  summary_data$alpha <- read.csv("analysis/alpha_diversity_statistics.csv")
}

# Volatility results
if (file.exists("analysis/volatility_statistics.csv")) {
  summary_data$volatility <- read.csv("analysis/volatility_statistics.csv")
}

# Instability results
if (file.exists("analysis/instability_statistics.csv")) {
  summary_data$instability <- read.csv("analysis/instability_statistics.csv")
}

# Generate summary report
summary_file <- "analysis/complete_analysis_summary.txt"
sink(summary_file)

cat("=== HOSPITAL SINK DISINFECTION ANALYSIS SUMMARY ===\n")
cat("Generated:", as.character(Sys.time()), "\n\n")

cat("STUDY DESIGN:\n")
cat("- Longitudinal study: 35 weeks\n")
cat("- Hospital rooms: ~31 rooms\n")
cat("- Sample locations: Tail Piece, P-Trap, Countertop\n")
cat("- Treatment: Virasept (H2O2/PAA) vs. Untreated control\n")
cat("- Sequencing: Full-length 16S rRNA (PacBio)\n")
cat("- Taxonomy: GTDB classification\n\n")

if (exists("meta_clean")) {
  cat("FINAL DATASET:\n")
  cat("- Total samples:", nrow(meta_clean), "\n")
  cat("- Genera analyzed:", ncol(comm_rel), "\n")

  sample_counts <- meta_clean %>%
    group_by(Sample_Type, Treatment) %>%
    summarise(n = n(), .groups = "drop")

  for (i in 1:nrow(sample_counts)) {
    cat(sprintf("- %s %s: %d samples\n",
                sample_counts$Sample_Type[i],
                sample_counts$Treatment[i],
                sample_counts$n[i]))
  }
  cat("\n")
}

cat("KEY FINDINGS:\n")

if (!is.null(summary_data$alpha)) {
  sig_alpha <- summary_data$alpha[summary_data$alpha$p_adjusted < 0.05, ]
  if (nrow(sig_alpha) > 0) {
    cat("Alpha Diversity:\n")
    for (i in 1:nrow(sig_alpha)) {
      result <- sig_alpha[i, ]
      change <- ((result$treated_mean - result$untreated_mean) / result$untreated_mean) * 100
      cat(sprintf("  - %s %s: %.1f%% change (p = %.3f)\n",
                  result$Sample_Type, result$Metric, change, result$p_adjusted))
    }
    cat("\n")
  }
}

if (!is.null(summary_data$volatility)) {
  cat("Community Volatility:\n")
  for (i in 1:nrow(summary_data$volatility)) {
    result <- summary_data$volatility[i, ]
    direction <- ifelse(result$treated_mean > result$untreated_mean, "higher", "lower")
    cat(sprintf("  - %s: %s volatility in treated (p = %.3f)\n",
                result$Sample_Type, direction, result$p_value))
  }
  cat("\n")
}

if (!is.null(summary_data$instability)) {
  cat("Community Instability:\n")
  for (i in 1:nrow(summary_data$instability)) {
    result <- summary_data$instability[i, ]
    direction <- ifelse(result$treated_mean < result$untreated_mean, "reduced", "increased")
    metric_name <- ifelse(result$metric == "lag1_acf", "predictability", "directional drift")
    cat(sprintf("  - %s: %s %s in treated (p = %.3f)\n",
                result$Sample_Type, direction, metric_name, result$p_value))
  }
}

cat("\nFILES GENERATED:\n")
analysis_files <- list.files("analysis", pattern = "\\.(csv|RData)$", full.names = FALSE)
figure_files <- list.files("figures/output", pattern = "\\.(pdf|png)$", full.names = FALSE)

cat("Analysis outputs:\n")
for (file in analysis_files) {
  cat("  -", file, "\n")
}

if (length(figure_files) > 0) {
  cat("Figure outputs:\n")
  for (file in figure_files) {
    cat("  -", file, "\n")
  }
}

sink()

cat("✓ Analysis summary saved to", summary_file, "\n")

# =============================================================================
# Pipeline Completion
# =============================================================================
cat("\n=== ANALYSIS PIPELINE COMPLETED SUCCESSFULLY ===\n")
cat("Completion time:", as.character(Sys.time()), "\n")
cat("\nNext steps:\n")
cat("1. Review analysis outputs in analysis/ directory\n")
cat("2. Check generated figures in figures/output/ directory\n")
cat("3. Read analysis summary:", summary_file, "\n")
cat("4. Adapt figure scripts as needed for publication\n")

cat("\n✓ All analyses completed!\n")