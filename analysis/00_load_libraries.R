# =============================================================================
# Hospital Sink Disinfection Study - Load Required Libraries
#
# Authors: Megan Stroup Hill, Claude Sonnet 4
# Description: Load all required R packages for analysis pipeline
# =============================================================================

cat("Loading required libraries for hospital sink disinfection analysis...\n")

# =============================================================================
# Core Data Manipulation
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr)         # Data manipulation
  library(tidyr)         # Data reshaping
  library(readr)         # File reading
  library(tibble)        # Modern data frames
})

# =============================================================================
# Statistical Analysis
# =============================================================================
suppressPackageStartupMessages({
  library(vegan)         # Ecological statistics and ordination
  library(ALDEx2)        # Differential abundance analysis
  library(lme4)          # Linear mixed-effects models
  library(broom)         # Tidy statistical outputs
})

# =============================================================================
# Data Visualization
# =============================================================================
suppressPackageStartupMessages({
  library(ggplot2)       # Grammar of graphics
  library(ggrepel)       # Text repelling for plots
  library(patchwork)     # Combining plots
  library(cowplot)       # Plot themes and arrangements
  library(gridExtra)     # Grid-based plot arrangements
  library(scales)        # Scale functions for ggplot2
})

# =============================================================================
# Color Palettes
# =============================================================================
suppressPackageStartupMessages({
  library(RColorBrewer)  # Color palettes
  library(viridis)       # Perceptually uniform color scales
})

# =============================================================================
# Specialized Analysis
# =============================================================================
suppressPackageStartupMessages({
  library(indicspecies)  # Indicator species analysis
  library(igraph)        # Network analysis
  library(corrplot)      # Correlation plots
})

# Resolve masked functions: MASS (attached via ALDEx2's dependencies) masks
# dplyr::select, and igraph masks vegan::diversity
select <- dplyr::select
diversity <- vegan::diversity

# =============================================================================
# Define Global Variables and Functions
# =============================================================================

# Treatment colors (consistent across all figures)
TREATMENT_COLORS <- c("Treated" = "#457B9D", "Untreated" = "#D62828")

# Sample type ordering (Tail Piece → P-Trap → Countertop)
SAMPLE_TYPE_LEVELS <- c("TP", "PT", "TOP")
SAMPLE_TYPE_LABELS <- c("Tail Piece", "P-Trap", "Countertop")

# ESKAPEE pathogen genera (including GTDB variants)
ESKAPEE_GENERA <- c(
  # Enterococcus variants
  "Enterococcus", "Enterococcus_A", "Enterococcus_B", "Enterococcus_C",
  "Enterococcus_D", "Enterococcus_E", "Enterococcus_G", "Enterococcus_I",
  # Staphylococcus
  "Staphylococcus",
  # Klebsiella
  "Klebsiella",
  # Acinetobacter
  "Acinetobacter",
  # Pseudomonas variants
  "Pseudomonas", "Pseudomonas_B", "Pseudomonas_E", "Pseudomonas_O",
  # Enterobacter
  "Enterobacter",
  # Escherichia
  "Escherichia"
)

# =============================================================================
# Utility Functions
# =============================================================================

# Standard genus-level filtering function
filter_genus_level <- function(taxa_table) {
  all_rows <- rownames(taxa_table)
  genus_rows <- grepl("g__[A-Za-z]", all_rows) &
    !grepl("s__", all_rows) &
    !grepl("g__[^-]+--__--[0-9a-f]", all_rows) &   # exclude ASV hashes
    !grepl("g__[^-]+--__$", all_rows)              # exclude genus--__ entries

  return(taxa_table[genus_rows, ])
}

# Clean genus names function
clean_genus_names <- function(genus_names) {
  # Extract genus from taxonomic string
  clean_names <- gsub(".*g__", "", genus_names)
  clean_names <- gsub("--.*", "", clean_names)
  return(clean_names)
}

# Format p-values for significance
format_p_value <- function(p) {
  case_when(
    p < 0.0001 ~ "****",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

# CLR transformation function
clr_transform <- function(x) {
  # Add pseudocount for zeros
  pseudocount <- min(x[x > 0]) * 0.5
  x_pseudo <- x + pseudocount
  # CLR transform
  log_x <- log(x_pseudo)
  clr_x <- log_x - rowMeans(log_x)
  return(clr_x)
}

# =============================================================================
# Set ggplot2 theme defaults
# =============================================================================
theme_set(theme_bw(base_size = 11))

# Update geom defaults for consistent styling
update_geom_defaults("point", list(size = 2, alpha = 0.7))
update_geom_defaults("line", list(linewidth = 0.8))

cat("✓ All libraries loaded successfully!\n")
cat("✓ Global variables and functions defined\n")
cat("✓ Plot themes configured\n")
cat("\nReady for analysis pipeline.\n")