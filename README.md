# Hospital Sink Disinfection Microbiome Study

This repository contains the complete analysis pipeline for the longitudinal study of hospital sink microbiomes under hydrogen peroxide/peracetic acid (Virasept) treatment.

## Study Overview

**Title**: Ecological replacement rather than sterilization drives hospital sink disinfection efficacy

**Study Design**: 35-week longitudinal study across 31 hospital rooms comparing treated vs. untreated sink microbiomes

**Sampling**: 3 locations per sink (tail piece, P-trap, countertop), 1,957 total samples

**Sequencing**: Full-length 16S rRNA gene sequencing (PacBio) with GTDB taxonomy

## Repository Structure

```
├── preprocessing/          # Sequence processing and normalization
├── analysis/              # Statistical analyses and data processing
├── figures/               # Figure generation scripts
├── data/                  # Example data files and metadata
├── docs/                  # Additional documentation
├── README.md              # This file
├── requirements.txt       # R package dependencies
└── environment.yml        # Conda environment specification
```

## Key Findings

1. **Spatial specificity**: Treatment effects strongest in tail pieces > P-traps > countertops
2. **Pathogen control**: ESKAPEE pathogens significantly enriched in untreated sinks
3. **Ecological replacement**: 85% of treatment-enriched genera were soil/water-associated
4. **Controlled instability**: Treatment disrupts biofilm succession without creating chaos
5. **Sustained effects**: Treatment efficacy maintained throughout 35-week study

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

3. **Run preprocessing** (if starting from raw sequences):
```bash
cd preprocessing/
bash 01_sequence_processing.sh
```

4. **Run main analyses**:
```r
# Load required libraries
source("analysis/00_load_libraries.R")

# Data preprocessing
source("analysis/01_data_preprocessing.R")

# Statistical analyses
source("analysis/02_alpha_diversity.R")
source("analysis/03_beta_diversity.R")
source("analysis/04_differential_abundance.R")
source("analysis/05_community_volatility.R")
source("analysis/06_community_instability.R")
source("analysis/07_temporal_dynamics.R")
```

5. **Generate figures**:
```r
# Main manuscript figures
source("figures/figure_1_volatility.R")
source("figures/figure_2_instability.R")
source("figures/figure_3_alpha_diversity.R")
source("figures/figure_4_top_taxa.R")
source("figures/figure_5_differential_abundance.R")
```

## Data Files

**⚠️ IMPORTANT: Large data files are not included in this repository due to GitHub size limits.**

The following data files must be downloaded separately and placed in the `data/` directory:

- `parsed_taxa_table.csv`: Genus-level abundance table (1,456 genera × 1,957 samples, ~498 MB)
- `metadata_run1_10.csv`: Sample metadata with treatment, location, and temporal information (~127 KB)

**To obtain data files:**
1. Download from NCBI SRA, BioProject ID: PRJNA1433472
2. Or contact: Megan Hill - msthoemmes@gmail.com

**After downloading, place files in:**
```
data/
├── parsed_taxa_table.csv    # ← Download and place here
├── metadata_run1_10.csv     # ← Download and place here
└── README.md                # Already included
```

## Credits

**Sequence Processing & Normalization**: Shan Sun
**Statistical Analysis & Figure Generation**: Megan Stroup Hill, Claude Sonnet 4
**Study Design & Interpretation**: [Principal Investigators]

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contact

For questions about this analysis, please contact:
- Megan Hill - msthoemmes@gmail.com

## Acknowledgments

- Shan Sun for sequence processing and initial data preparation
- Claude Sonnet 4 for analysis pipeline development and optimization
- This work was supported primarily by the Engineering Research Centers Program of the National Science Foundation under NSF Cooperative Agreement No. EEC-2133504