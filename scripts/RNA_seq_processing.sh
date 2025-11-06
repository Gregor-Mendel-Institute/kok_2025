#!/bin/bash
#
#SBATCH --job-name=RNA_analysis
#SBATCH --output=output/RNA_alignment_%A_%a.out
#SBATCH --error=output/RNA_alignment_%A_%a.err
#
#SBATCH --time=6:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --qos=medium

# -----------------------------------------------
# Script for mRNA-seq data analysis pipeline
# Author: Jian Yi Kok, Zachary Harvey
# Date: [2025-07-17]
# Description: This script processes mRNA-seq data, including quality control,
#              adaptor trimming, and transcript quantification using Kallisto.
# -----------------------------------------------

# 1. Set up directories
work_folder=/scratch/RNA_analysis
mkdir -p $work_folder

# 2. Define paths to required resources
index_kallisto=/path/to/kallisto_index
index_star=/path/to/star_index
input_dir=/path/to/demultiplexed_reads
output_dir=/path/to/output_directory
mkdir -p $output_dir

# 3. Define sample IDs (replace with actual sample IDs)
samples=(sample1 sample2 sample3 sample4 sample5 sample6 sample7 sample8 sample9 sample10)
output=${samples[$SLURM_ARRAY_TASK_ID]}

# 4. Load required modules
module load build-env/2020
module load fastqc/0.11.8-java-1.8
module load samtools/1.9-foss-2018b
module load bedtools/2.27.1-foss-2018b
module load trimmomatic/0.38-java-1.8
module load kallisto/0.46.0-foss-2018b
module load star/2.7.1a-foss-2018b

# 5. Locate input files
read1=$(ls $input_dir/${output}/* | grep '_R1_.*\.fastq.gz')
read2=$(ls $input_dir/${output}/* | grep '_R2_.*\.fastq.gz')

# 6. Quality control with FastQC
qc_dir=$output_dir/fastQC/$output
mkdir -p $qc_dir
fastqc -t 2 -o $qc_dir $read1 $read2

# 7. Adaptor trimming with Trimmomatic
trimmed_dir=$work_folder/trimmed_reads
mkdir -p $trimmed_dir
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.38.jar PE -threads 16 \
    $read1 $read2 \
    $trimmed_dir/${output}_R1_paired.fq.gz \
    $trimmed_dir/${output}_R1_unpaired.fq.gz \
    $trimmed_dir/${output}_R2_paired.fq.gz \
    $trimmed_dir/${output}_R2_unpaired.fq.gz \
    ILLUMINACLIP:/path/to/adaptor_sequences/TruSeq3-PE-2.fa:2:30:10:2:keepBothReads \
    LEADING:3 TRAILING:3 MINLEN:36 SLIDINGWINDOW:4:20

# 8. Transcript quantification with Kallisto
kallisto_dir=$output_dir/kallisto_counts/$output
mkdir -p $kallisto_dir
kallisto quant -i $index_kallisto -o $kallisto_dir -b 100 -t 16 --rf-stranded \
    $trimmed_dir/${output}_R1_paired.fq.gz \
    $trimmed_dir/${output}_R2_paired.fq.gz

# 9. Completion message
echo "Analysis for sample $output completed successfully."