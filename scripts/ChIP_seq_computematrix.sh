#!/bin/bash

# -----------------------------------------------
# Script for ChIP-seq computeMatrix Analysis at TSS, TES and scaled regions
# Author: [Jian Yi Kok]
# Date: [2025-07-17]
# Description: This script calculates signal matrices for ChIP-seq data
#              at TSS, TES and scaled regions (TSS to TES) using deepTools'
#              computeMatrix and generates heatmaps using plotHeatmap.
# -----------------------------------------------

# Load required modules
ml deeptools/3.3.1-foss-2018b-python-3.6.6

# Base directories
BASE_DIR="/path/to/project"  # Replace with the base directory of your project
BIGWIG_DIR="$BASE_DIR/bigwigs"  # Directory containing bigWig files
REGIONS_FILE="/path/to/annotations/pcg.bed"  # Replace with the path to your BED file
OUTPUT_DIR="$BASE_DIR/deepTools/computeMatrix"  # Output directory for computeMatrix results

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Input bigWig files (12 files)
BIGWIG_FILES=(
    "$BIGWIG_DIR/sample1.bigwig"
    "$BIGWIG_DIR/sample2.bigwig"
    "$BIGWIG_DIR/sample3.bigwig"
    "$BIGWIG_DIR/sample4.bigwig"
    "$BIGWIG_DIR/sample5.bigwig"
    "$BIGWIG_DIR/sample6.bigwig"
    "$BIGWIG_DIR/sample7.bigwig"
    "$BIGWIG_DIR/sample8.bigwig"
    "$BIGWIG_DIR/sample9.bigwig"
    "$BIGWIG_DIR/sample10.bigwig"
    "$BIGWIG_DIR/sample11.bigwig"
    "$BIGWIG_DIR/sample12.bigwig"
)

# Labels for samples (12 labels corresponding to the bigWig files)
SAMPLES_LABELS=("Sample1" "Sample2" "Sample3" "Sample4" "Sample5" "Sample6" "Sample7" "Sample8" "Sample9" "Sample10" "Sample11" "Sample12")

# ComputeMatrix parameters
BEFORE_REGION_LENGTH=1000
AFTER_REGION_LENGTH=1000
THREADS=max

# -----------------------------------------------
# Functions
# -----------------------------------------------

# Function to run computeMatrix for reference-point (TSS or TES)
run_computeMatrix_referencePoint() {
    local reference_point=$1  # TSS or TES
    local output_prefix=$2    # Prefix for output files

    computeMatrix reference-point \
        --referencePoint "$reference_point" \
        --scoreFileName "${BIGWIG_FILES[@]}" \
        --regionsFileName "$REGIONS_FILE" \
        -out "${OUTPUT_DIR}/${output_prefix}.mat.gz" \
        --outFileNameMatrix "${OUTPUT_DIR}/${output_prefix}.tab" \
        --beforeRegionStartLength "$BEFORE_REGION_LENGTH" \
        --afterRegionStartLength "$AFTER_REGION_LENGTH" \
        --missingDataAsZero \
        -p "$THREADS" \
        --samplesLabel "${SAMPLES_LABELS[@]}"
}

# Function to run computeMatrix for scale-regions (TSS to TES)
run_computeMatrix_scaleRegions() {
    local output_prefix=$1  # Prefix for output files

    computeMatrix scale-regions \
        --scoreFileName "${BIGWIG_FILES[@]}" \
        --regionsFileName "$REGIONS_FILE" \
        -out "${OUTPUT_DIR}/${output_prefix}.mat.gz" \
        --outFileNameMatrix "${OUTPUT_DIR}/${output_prefix}.tab" \
        --regionBodyLength 5000 \  # Length of the scaled region (adjust as needed)
        --beforeRegionStartLength "$BEFORE_REGION_LENGTH" \
        --afterRegionStartLength "$AFTER_REGION_LENGTH" \
        --missingDataAsZero \
        -p "$THREADS" \
        --samplesLabel "${SAMPLES_LABELS[@]}"
}

# Function to run plotHeatmap
run_plotHeatmap() {
    local matrix_file=$1  # Input matrix file
    local output_file=$2  # Output heatmap file

    plotHeatmap \
        -m "$matrix_file" \
        -out "$output_file" \
        --colorMap "viridis" \
        --whatToShow "heatmap and colorbar" \
        --zMin 0 \
        --zMax 10 \
        --heatmapHeight 10 \
        --heatmapWidth 6
}

# -----------------------------------------------
# Main Script
# -----------------------------------------------

# Run computeMatrix for TSS
echo "Running computeMatrix for TSS..."
run_computeMatrix_referencePoint "TSS" "TSS_analysis"

# Generate heatmap for TSS
echo "Generating heatmap for TSS..."
run_plotHeatmap "${OUTPUT_DIR}/TSS_analysis.mat.gz" "${OUTPUT_DIR}/TSS_analysis_heatmap.pdf"

# Run computeMatrix for TES
echo "Running computeMatrix for TES..."
run_computeMatrix_referencePoint "TES" "TES_analysis"

# Generate heatmap for TES
echo "Generating heatmap for TES..."
run_plotHeatmap "${OUTPUT_DIR}/TES_analysis.mat.gz" "${OUTPUT_DIR}/TES_analysis_heatmap.pdf"

# Run computeMatrix for scale-regions (TSS to TES)
echo "Running computeMatrix for scale-regions (TSS to TES)..."
run_computeMatrix_scaleRegions "ScaleRegions_analysis"

# Generate heatmap for scale-regions
echo "Generating heatmap for scale-regions..."
run_plotHeatmap "${OUTPUT_DIR}/ScaleRegions_analysis.mat.gz" "${OUTPUT_DIR}/ScaleRegions_analysis_heatmap.pdf"

echo "Analysis completed. Output files saved to: $OUTPUT_DIR"