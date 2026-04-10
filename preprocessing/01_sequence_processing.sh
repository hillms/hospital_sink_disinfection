#!/bin/bash

# =============================================================================
# Hospital Sink Disinfection Study - Sequence Processing Pipeline
#
# Author: Shan Sun
# Description: QIIME2-based processing of PacBio CCS sequences
# Input: Raw fastq files with manifest
# Output: Processed feature table and taxonomy assignments
# =============================================================================

set -e  # Exit on any error

# Create output directories
mkdir -p result/{rep-seqs,stats,table,taxonomy}

echo "Starting sequence processing pipeline..."
echo "Timestamp: $(date)"

# =============================================================================
# Step 1: Import sequences into QIIME2
# =============================================================================
echo "Step 1: Importing sequences..."

qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path manifest.txt \
  --output-path single-end-demux.qza \
  --input-format SingleEndFastqManifestPhred33V2

echo "✓ Sequences imported successfully"

# =============================================================================
# Step 2: Denoise with DADA2 (optimized for PacBio CCS)
# =============================================================================
echo "Step 2: Denoising sequences with DADA2..."

qiime dada2 denoise-ccs \
	--i-demultiplexed-seqs single-end-demux.qza \
	--p-front AGRGTTYGATYMTGGCTCAG \
	--p-adapter RGYTACCTTGTTACGACTT \
	--p-min-len 1000 \
	--p-max-len 1600 \
	--p-n-threads 32 \
	--o-representative-sequences result/rep-seqs.qza \
	--o-table result/table.qza \
	--o-denoising-stats result/stats.qza

echo "✓ Denoising completed"

# =============================================================================
# Step 3: Taxonomic classification using GTDB
# =============================================================================
echo "Step 3: Assigning taxonomy with GTDB classifier..."

qiime feature-classifier classify-sklearn \
      --i-classifier /projects/afodor_premier1/ssun5/database/gtdb_classifier_r220.qza \
      --i-reads result/rep-seqs.qza \
      --o-classification result/taxonomy.qza

echo "✓ Taxonomy assigned"

# =============================================================================
# Step 4: Export QIIME2 artifacts
# =============================================================================
echo "Step 4: Exporting QIIME2 artifacts..."

# Export representative sequences
qiime tools export \
  --input-path result/rep-seqs.qza \
  --output-path result/rep-seqs

# Export denoising statistics
qiime tools export \
  --input-path result/stats.qza \
  --output-path result/stats

# Export feature table
qiime tools export \
  --input-path result/table.qza \
  --output-path result/table

# Export taxonomy
qiime tools export \
  --input-path result/taxonomy.qza \
  --output-path result/taxonomy

echo "✓ All artifacts exported"

# =============================================================================
# Step 5: Create BIOM table with taxonomy
# =============================================================================
echo "Step 5: Creating BIOM table with taxonomy..."

# Add taxonomy metadata to BIOM table
biom add-metadata \
 --input-fp result/table/feature-table.biom \
 --observation-metadata-fp result/taxonomy/taxonomy.tsv \
 --output-fp result/biom-with-taxonomy.biom  \
 --sc-separated taxonomy --observation-header OTUID,taxonomy

# Convert BIOM to TSV format
biom convert -i result/biom-with-taxonomy.biom \
             -o result/biom-with-taxonomy.txt \
             --to-tsv --header-key taxonomy

echo "✓ BIOM table with taxonomy created"

# =============================================================================
# Step 6: Generate processing summary
# =============================================================================
echo "Step 6: Generating processing summary..."

echo "=== SEQUENCE PROCESSING SUMMARY ===" > processing_summary.txt
echo "Date: $(date)" >> processing_summary.txt
echo "Pipeline: QIIME2 + DADA2 + GTDB" >> processing_summary.txt
echo "" >> processing_summary.txt

# Get sample counts and feature counts from BIOM file
echo "Output files:" >> processing_summary.txt
echo "- Feature table: result/table/feature-table.biom" >> processing_summary.txt
echo "- Taxonomy: result/taxonomy/taxonomy.tsv" >> processing_summary.txt
echo "- Combined: result/biom-with-taxonomy.txt" >> processing_summary.txt
echo "- Representative sequences: result/rep-seqs/dna-sequences.fasta" >> processing_summary.txt

echo "✓ Processing summary generated"
echo ""
echo "=== PIPELINE COMPLETED SUCCESSFULLY ==="
echo "Timestamp: $(date)"
echo ""
echo "Next steps:"
echo "1. Transfer result/biom-with-taxonomy.txt to R analysis pipeline"
echo "2. Apply normalization in R (see 02_normalization.R)"
echo "3. Run downstream statistical analyses"