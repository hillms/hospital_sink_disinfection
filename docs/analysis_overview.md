# Hospital Sink Disinfection Study - Analysis Overview

## Study Background

Hospital-acquired infections represent a major public health challenge, with hospital sink drainage systems serving as reservoirs for opportunistic pathogens. This study investigates the ecological mechanisms underlying effective antimicrobial control using a hydrogen peroxide/peracetic acid treatment (Virasept).

## Novel Analytical Approaches

### 1. Community Volatility Analysis
- **Concept**: Week-to-week compositional turnover
- **Method**: Aitchison distances between consecutive timepoints
- **Innovation**: Reveals episodic treatment effects vs. baseline differences
- **Key finding**: Similar overall volatility but different weekly response patterns

### 2. Community Instability Framework
- **Lag-1 Autocorrelation**: Measures short-term predictability
- **PERMANOVA R²**: Quantifies directional drift over time
- **Joint analysis**: Creates instability landscape with four quadrants
- **Controlled instability**: Treatment signature of low ACF + low R²

### 3. Compositional Shift Under Treatment
- **Traditional view**: Sterilization eliminates all microbes
- **Our finding**: Treated plumbing communities are relatively enriched in environmental taxa and depleted in ESKAPEE genera
- **Open question**: Relative abundance data cannot distinguish competitive exclusion from differential survival of oxidant-tolerant bacteria

## Key Methodological Choices

### Contaminant Removal
- **decontam (prevalence method)** against 12 negative-control libraries
- **Removal rule**: score P < 0.1 and present in >= 50% of the controls of every site type; 2 ASVs removed (1.4% of reads)
- **Frequency/combined modes not used**: treatment lowers DNA yield, so the frequency test flags taxa genuinely enriched in treated sinks

### Read Depth and Rarefaction
- Read depth is confounded with treatment (median tail-piece reads 1,459 treated vs 67,339 untreated)
- **Alpha diversity**: raw genus counts rarefied to 1,000 reads
- **Differential abundance**: ALDEx2 on raw (unrarefied) genus counts

### Differential Abundance Calibration
- Treatment was assigned by room and each room sampled ~35 times, so samples are not independent
- Treatment labels shuffled among the 31 rooms remove any true effect; a well-calibrated method should then call ~0 genera
- ALDEx2 called 0-6 genera with shuffled labels; DESeq2 (6-27) and edgeR (30-74) were poorly calibrated, so only ALDEx2 is reported
- Tail piece / P-trap agreement is measured as genera significant at both sites in the same direction (31 with true labels vs median 0 with shuffled labels), not by effect-size correlation, which shuffled labels also produce

### Compositional Data Analysis
- **CLR transformation**: Properly handles compositional data constraints
- **Aitchison distances**: Appropriate metric for microbiome data
- **Pseudocount strategy**: Handles zero values in sparse datasets

### Temporal Dynamics
- **Per-room analysis**: Accounts for spatial heterogeneity
- **Autocorrelation functions**: Reveals temporal structure
- **Weekly PERMANOVA**: Identifies episodic effects
- **Temporal clustering**: Detects facility-wide perturbations

### Statistical Rigor
- **Multiple testing correction**: Benjamini-Hochberg procedure
- **Robust methods**: Wilcoxon tests for non-parametric data
- **Spatial controls**: Countertop samples as negative controls

## Ecological Insights

### Spatial Gradient
Treatment effects follow a clear spatial pattern:
- **Tail Piece** > **P-Trap** > **Countertop**
- Direct application site shows strongest effects
- Spatial attenuation due to dilution and physical barriers

### Temporal Patterns
- **Immediate onset**: Effects detectable from week 1
- **Sustained efficacy**: Maintained throughout 35-week study
- **Episodic responses**: Weekly variations suggest perturbation-response dynamics

### Pathogen Control
- **ESKAPEE enrichment**: Enterobacter, Klebsiella and Acinetobacter relatively enriched in untreated sinks
- **Environmental taxa**: Treated plumbing relatively enriched in environmental genera such as Afipia, Bradyrhizobium and Pelomonas
- **Caveat**: Relative abundances; treated samples yield less DNA, so enrichment does not imply higher absolute abundance

## Clinical Relevance

### Healthcare Applications
- **Environmental management**: Supports evaluating cleaning by its effect on community composition, not only by bioburden reduction
- **Infection prevention**: Lower relative abundance of ESKAPEE genera in treated plumbing

### Broader Implications
- **Antimicrobial strategy**: Shifting community composition as a complement to broad-spectrum elimination

## Future Directions

### Mechanistic Studies
- **Metabolic profiling**: Understand biochemical basis of ecological replacement
- **Strain-level analysis**: Identify specific beneficial organisms
- **Interaction networks**: Map ecological relationships

### Applied Research
- **Optimization studies**: Refine treatment protocols for maximum efficacy
- **Broader applications**: Test in other healthcare environments
- **Long-term monitoring**: Assess sustainability over extended periods

### Technological Development
- **Biomarker identification**: Early warning systems for pathogen colonization
- **Ecological engineering**: Design treatments to promote beneficial communities
- **Personalized approaches**: Tailor treatments to specific environmental conditions

## Statistical Framework Summary

The analysis pipeline implements a comprehensive statistical framework:

1. **Contaminant removal**: decontam prevalence method (`preprocessing/02_decontam.R`)
2. **Data preprocessing**: Normalized and raw genus tables, sample filtering
3. **Alpha diversity**: Observed genera and Shannon, rarefied to 1,000 reads
4. **Beta diversity**: Aitchison PCA (Figure 3B)
5. **Differential abundance**: ALDEx2 on raw counts, with a room-shuffle calibration check
6. **Community volatility**: Temporal compositional turnover
7. **Community instability**: Temporal autocorrelation and drift
8. **Weekly PERMANOVA**: Treatment effect by week (Table S3)
9. **Visualization and tables**: Figures 1-5, Table 1, Tables S1-S3

This framework provides a template for similar microbiome studies in healthcare and other built environments.