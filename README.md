# Hospital Sink Disinfection Microbiome Study

This repository contains the complete analysis pipeline for the longitudinal study of hospital sink microbiomes under hydrogen peroxide/peracetic acid (Virasept) treatment.

## Study Overview

**Title**: Ecological replacement rather than sterilization drives hospital sink disinfection efficacy

**Study Design**: 35-week longitudinal study across 31 hospital rooms comparing treated vs. untreated sink microbiomes

**Sampling**: 3 locations per sink (tail piece, P-trap, countertop), 2,186 sink samples analysed

**Sequencing**: Full-length 16S rRNA gene sequencing (PacBio) with GTDB taxonomy

## Repository Structure

```
├── preprocessing/          # Sequence processing, decontam, taxa tables
├── analysis/              # Statistical analyses and data processing
├── figures/               # Figure generation scripts
├── tables/                # Table 1 and Tables S1-S3
├── data/                  # Input tables and metadata (not tracked)
├── docs/                  # Additional documentation
├── README.md              # This file
├── requirements.txt       # R package dependencies
└── environment.yml        # Conda environment specification
```

## Key Findings

1. **Spatial specificity**: Treatment shifts tail-piece and P-trap communities; countertops are unaffected
2. **Pathogen control**: Genera containing ESKAPEE pathogens (Enterobacter, Klebsiella, Acinetobacter) are relatively enriched in untreated sinks
3. **Compositional shift**: ALDEx2 finds 51 (tail piece) and 75 (P-trap) differentially abundant genera; 31 shift in the same direction at both sites, against a median of 0 with room-shuffled treatment labels
4. **Richness**: Observed genera about 20% lower in treated plumbing (rarefied to 1,000 reads)
5. **Temporal dynamics**: Treated tail pieces show lower week-to-week predictability (lag-1 autocorrelation) and less directional drift

## Usage

### Prerequisites

- R (>= 4.0.0)
- QIIME2 (>= 2023.5) for preprocessing
- Required R packages (see `requirements.txt`)

### Quick Start

1. **Clone the repository**:
```bash
git clone https://github.com/hillms/hospital_sink_disinfection.git
cd hospital_sink_disinfection
```

2. **Install dependencies**:
```r
# Install required R packages
source("requirements.txt")
```

3. **Run preprocessing** (if starting from raw sequences), from the repository root:
```bash
bash preprocessing/01_sequence_processing.sh   # QIIME2: result/biom-with-taxonomy.txt
Rscript preprocessing/02_decontam.R            # remove contaminant ASVs (decontam)
Rscript preprocessing/03_build_taxa_tables.R   # data/parsed_taxa_table.csv, data/genus_counts_raw.csv
```

4. **Run main analyses**:
```r
# Load required libraries
source("analysis/00_load_libraries.R")

# Data preprocessing
source("analysis/01_data_preprocessing.R")

# Statistical analyses
source("analysis/02_alpha_diversity.R")          # rarefied observed genera + Shannon
source("analysis/04_differential_abundance.R")   # ALDEx2 on raw genus counts
source("analysis/04b_da_calibration.R")          # room-shuffle calibration (~10-15 min)
source("analysis/05_community_volatility.R")
source("analysis/06_community_instability.R")
source("analysis/07_weekly_permanova.R")         # Table S3
```

Or run everything with `source("run_complete_analysis.R")` (set `RUN_CALIBRATION <- FALSE` first to skip the calibration step).

5. **Generate figures**:
```r
# Main manuscript figures
source("figures/figure_1_volatility.R")
source("figures/figure_2_instability.R")
source("figures/figure_3_alpha_diversity.R")
source("figures/figure_4_top_taxa.R")
source("figures/figure_5_differential_abundance.R")

# Table 1 and Tables S1-S3
source("tables/make_tables.R")
```

## Methods Notes

- **Contaminant removal**: decontam prevalence method against 12 negative controls. ASVs are removed if they score P < 0.1 and occur in at least 50% of the controls of every site type. Two ASVs are removed (Bradyrhizobium and Qipengyuania; 1.4% of reads). All flagged ASVs are listed in Table S1.
- **Read depth**: Treated samples have far fewer reads than untreated ones, so alpha diversity uses genus counts rarefied to 1,000 reads, and differential abundance uses ALDEx2 on raw counts. Other analyses use the depth-normalized table.
- **Differential abundance calibration**: Treatment was assigned by room, so treatment labels were shuffled among the 31 rooms to check how many genera each method calls with no true effect. ALDEx2 called 0-6 genera; DESeq2 and edgeR called many more (Table S2), so only ALDEx2 is reported.

## Revision History

**October 2026 (manuscript revision)**
- Added contaminant removal with decontam (`preprocessing/02_decontam.R`).
- Fixed sample matching for Batch 9: 229 sink samples had been dropped because their feature-table IDs carry `_corrected`. The analysis now uses 2,186 samples (previously 1,957).
- `preprocessing/02_normalization.R` was replaced by `03_build_taxa_tables.R`. The new script reproduces the parsed table format the analysis scripts expect, and also writes raw genus counts.
- Alpha diversity now uses observed genera and Shannon on counts rarefied to 1,000 reads. Unrarefied Chao1 was dropped; it largely tracked read depth.
- Differential abundance now uses ALDEx2 on raw counts (previously on depth-normalized counts), with a room-shuffle calibration check (`analysis/04b_da_calibration.R`).
- Added the weekly PERMANOVA (`analysis/07_weekly_permanova.R`) and table scripts (`tables/make_tables.R`).
- Removed the "85% soil/water-associated" claim and the Figure 5 subtitle; the key is now in the figure legend.

## Data Files

**⚠️ IMPORTANT: Large data files are not included in this repository due to GitHub size limits.**

The following data files must be downloaded separately and placed in the `data/` directory:

- `metadata_run1_10.csv`: Sample metadata with treatment, location, and temporal information (~127 KB)
- `parsed_taxa_table.csv`: Depth-normalized abundance table (~500 MB), built by `preprocessing/03_build_taxa_tables.R`
- `genus_counts_raw.csv`: Raw read counts per genus (1,412 genera), built by `preprocessing/03_build_taxa_tables.R`

**To obtain data files:**
1. Download from NCBI SRA, BioProject ID: PRJNA1433472
2. Or contact: Megan Hill - msthoemmes@gmail.com

**After downloading, place files in:**
```
data/
├── metadata_run1_10.csv     # ← Download and place here
├── parsed_taxa_table.csv    # ← Built by preprocessing/03_build_taxa_tables.R
├── genus_counts_raw.csv     # ← Built by preprocessing/03_build_taxa_tables.R
└── README.md                # Already included
```

## Credits

**Sequence Processing & Normalization**: Shan Sun
**Statistical Analysis & Figure Generation**: Megan Stroup Hill, Claude Sonnet 4, Claude
**Study Design & Interpretation**: [Principal Investigators]

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contact

For questions about this analysis, please contact:
- Megan Hill - msthoemmes@gmail.com

## Acknowledgments

- Shan Sun for sequence processing and initial data preparation
- Claude Sonnet 4 and Claude for analysis pipeline development and optimization
- This work was supported primarily by the Engineering Research Centers Program of the National Science Foundation under NSF Cooperative Agreement No. EEC-2133504