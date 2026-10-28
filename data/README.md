# Data Directory

This directory should contain the input data files for the analysis pipeline.

## Required Files

### `metadata_run1_10.csv`
- **Description**: Sample metadata with treatment and temporal information
- **Format**: CSV with samples as rows
- **Required columns**:
  - `SeqID`: Sample identifier
  - `Sample`: Sample name
  - `Batch`: Sequencing batch (batch1-batch10)
  - `Room`: Hospital room identifier (`C-`, `C` and `NegCon` are negative controls)
  - `Week`: Week number (1-35, some contain "P" and "Y" - filtered out)
  - `Sample_Type`: Sample location (TP=tail piece, PT=P-trap, TOP=countertop)
  - `Virasept_Treated`: Treatment status ("Yes"/"No")

### `parsed_taxa_table.csv`
- **Description**: Depth-normalized abundance table (proportional scaling to the mean library size) in parsed-taxonomy format; `analysis/01_data_preprocessing.R` keeps the genus-level rows
- **Source**: `preprocessing/03_build_taxa_tables.R`, from the decontaminated feature table
- **Used for**: volatility, instability, beta diversity, top taxa, weekly PERMANOVA

### `genus_counts_raw.csv`
- **Description**: Raw read counts per GTDB genus (genera as rows, samples as columns)
- **Source**: `preprocessing/03_build_taxa_tables.R`, from the decontaminated feature table
- **Used for**: alpha diversity (rarefied to 1,000 reads) and ALDEx2 differential abundance

## Building the Tables

From the repository root, after `preprocessing/01_sequence_processing.sh` has written `result/biom-with-taxonomy.txt`:

```bash
Rscript preprocessing/02_decontam.R          # result/decontam/ (contaminant ASVs removed)
Rscript preprocessing/03_build_taxa_tables.R # data/parsed_taxa_table.csv, data/genus_counts_raw.csv
```

The Batch 9 columns of the feature table carry `_FastQ_corrected_` in their names; `02_decontam.R` renames them to match the metadata. Without this, the 229 Batch 9 sink samples are dropped.

## Data Access

Due to file size limitations, these data files are not included in the repository. To obtain the data files:

1. Contact the corresponding authors
2. Download files from the NCBI SRA, BioProject ID: PRJNA1433472

## Data Processing Notes

- **Contaminants**: 2 ASVs removed with decontam (prevalence method; Bradyrhizobium and Qipengyuania, 1.4% of sink-sample reads)
- **Taxonomy**: GTDB classification (release 220)
- **Controls**: Negative controls and samples with missing treatment or non-numeric week removed
- **Missing weeks**: Weeks 19, 23, and 25 have no data

## Final Dataset

- Sink samples analysed: 2,186
- Genera with reads in the analysed samples: 1,412
- Study duration: 35 weeks
