#!/bin/bash

# -----------------------------------------------
# Script for MNase-seq data analysis
# Author: Jian Yi Kok
# Date: [2025-07-17]
# Description: This script generates matrices of MNase-seq signal
#              using deepTools.
# -----------------------------------------------

# Load required modules
ml deeptools/3.3.1-foss-2018b-python-3.6.6

# Figure 1B
# 1. Define input files and directories
wt_bigwig=/path/to/WT.mRp.clN.Fnor.smooth.bigWig
remodeler_bigwigs=(
    /path/to/fft1.mRp.clN.Fnor.smooth.bigWig
    /path/to/fft2.mRp.clN.Fnor.smooth.bigWig
    /path/to/fft3.mRp.clN.Fnor.smooth.bigWig
    /path/to/hrp1.mRp.clN.Fnor.smooth.bigWig
    /path/to/hrp3.mRp.clN.Fnor.smooth.bigWig
    /path/to/mit1.mRp.clN.Fnor.smooth.bigWig
    /path/to/rrp1.mRp.clN.Fnor.smooth.bigWig
    /path/to/rrp2.mRp.clN.Fnor.smooth.bigWig
    /path/to/snf22.mRp.clN.Fnor.smooth.bigWig
    /path/to/irc20.mRp.clN.Fnor.smooth.bigWig
    /path/to/swr1.mRp.clN.Fnor.smooth.bigWig
)
remodeler_labels=("fft1" "fft2" "fft3" "hrp1" "hrp3" "mit1" "rrp1" "rrp2" "snf22" "irc20" "swr1")
regions_file=/path/to/plusonenuc_filtered_final_WT.bed
output_dir=/path/to/output_directory
mkdir -p $output_dir

# 2. Define parameters
before_region_length=800
after_region_length=800

# 3. Loop through each remodeler and compare to WT
for i in "${!remodeler_bigwigs[@]}"; do
    remodeler_bigwig=${remodeler_bigwigs[$i]}
    remodeler_label=${remodeler_labels[$i]}
    
    # Define output files
    output_matrix=$output_dir/${remodeler_label}_vs_WT.mat.gz
    output_tab=$output_dir/${remodeler_label}_vs_WT.tab

    # Run computeMatrix
    computeMatrix reference-point \
        --referencePoint TSS \
        --scoreFileName $wt_bigwig $remodeler_bigwig \
        --regionsFileName $regions_file \
        -out $output_matrix \
        --outFileNameMatrix $output_tab \
        --beforeRegionStartLength $before_region_length \
        --afterRegionStartLength $after_region_length \
        --missingDataAsZero \
        -p max \
        --samplesLabel "WT" "$remodeler_label"

    echo "Matrix generation for $remodeler_label vs WT completed."
done

# Completion message
echo "All remodeler comparisons to WT completed successfully."


# Figure 1C
computeMatrix reference-point \
    --referencePoint TSS \
    --scoreFileName /path/to/WT.mRp.clN.Fnor.smooth.bigWig \
                    /path/to/hrp1.mRp.clN.Fnor.smooth.bigWig \
                    /path/to/hrp3.mRp.clN.Fnor.smooth.bigWig \
                    /path/to/hrp1hrp3.mRp.clN.Fnor.smooth.bigWig \
    --regionsFileName /path/to/plusonenuc_filtered_final_WT.bed \
    --outFileName /path/to/output/CHD1_mnase.mat.gz \
    --outFileNameMatrix /path/to/output/CHD1_mnase.tab \
    --beforeRegionStartLength 800 \
    --afterRegionStartLength 800 \
    --missingDataAsZero \
    -p max \
    --samplesLabel "WT" "hrp1" "hrp3" "hrp1hrp3"
done


# Figure 3D
computeMatrix reference-point \
    --referencePoint TSS \
    --scoreFileName /path/to/WT.mRp.clN.Fnor.smooth.bigWig \
                    /path/to/hrp1.mRp.clN.Fnor.smooth.bigWig \
                    /path/to/hrp3.mRp.clN.Fnor.smooth.bigWig \
                    /path/to/hrp1hrp3.mRp.clN.Fnor.smooth.bigWig \
    --regionsFileName /path/to/as_genes.bed \
    --outFileName /path/to/output/as_mnase.mat.gz \
    --outFileNameMatrix /path/to/output/as_mnase.tab \
    --beforeRegionStartLength 800 \
    --afterRegionStartLength 800 \
    --missingDataAsZero \
    -p max \
    --samplesLabel "WT" "hrp1" "hrp3" "hrp1hrp3"
done


# Figure 7F
computeMatrix reference-point \
    --referencePoint TSS \
    --scoreFileName /path/to/WT.mRp.clN.Fnor.smooth.bigWig \
                    /path/to/prf1.mRp.clN.Fnor.smooth.bigWig \
                    /path/to/hrp1hrp3prf1.mRp.clN.Fnor.smooth.bigWig \
    --regionsFileName /path/to/plusonenuc_filtered_final_WT.bed \
    --outFileName /path/to/output/prf1_mnase.mat.gz \
    --outFileNameMatrix /path/to/output/prf1_mnase.tab \
    --beforeRegionStartLength 800 \
    --afterRegionStartLength 800 \
    --missingDataAsZero \
    -p max \
    --samplesLabel "WT" "prf1" "hrp1hrp3prf1"
done