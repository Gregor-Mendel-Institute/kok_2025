#!/bin/bash
#SBATCH --job-name=minimap2_trans
#SBATCH --output=minimap2_%j.out
#SBATCH --error=minimap2_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G

# Load required modules
module load build-env/f2022
module load minimap2/2.24-gcccore-11.2.0
module load samtools/1.17-gcc-12.2.0

# Define input and output directories
OUTPUT_DIR="/scratch-cbe/users/elin.axelsson/minimap2_out_fasta" # fa files on scratch 
#REFERENCE="/groups/berger/user/elin.axelsson/projects/projects_2025/pacbio/paper25/annotations/Schizosaccharomyces_pombe.ASM294v2.dna.toplevel.fa"
REFERENCE="/groups/berger/user/elin.axelsson/projects/projects_2025/ONS/Schizosaccharomyces_pombe_all_chromosomes.fa"


# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# List of full paths to FASTA files
FASTA_FILES=(
        "/scratch-cbe/users/elin.axelsson/IsoSeq_Results/library1/library1.transcripts.fasta"
        "/scratch-cbe/users/elin.axelsson/IsoSeq_Results/library2/library2.transcripts.fasta"
)

# Loop through the FASTA files
for FASTA_FILE in "${FASTA_FILES[@]}"; do
    # gunzip the FASTA file if it is gzipped
    if [[ $FASTA_FILE == *.gz ]]; then
        gunzip -c $FASTA_FILE > "${FASTA_FILE%.gz}"
        FASTA_FILE="${FASTA_FILE%.gz}"
    fi
    # Get the base name of the FASTA file without the extension
    # and create output file names
    BASENAME=$(basename $FASTA_FILE .fasta)
    echo "Basename: $BASENAME"
    SAM_FILE="$OUTPUT_DIR/${BASENAME}.sam"
    FINAL_BAM_FILE="$OUTPUT_DIR/${BASENAME}.bam"
    SORTED_BAM_FILE="$OUTPUT_DIR/${BASENAME}_sorted.bam"

    echo "Processing ${FASTA_FILE}..."
    
    # Run minimap2 directly on the FASTA file
    minimap2 -ax splice:hq -uf $REFERENCE $FASTA_FILE > $SAM_FILE
    
    # Convert SAM to BAM
    samtools view -S -b $SAM_FILE > $FINAL_BAM_FILE
    
    # Sort the BAM file
    samtools sort $FINAL_BAM_FILE -o $SORTED_BAM_FILE
    
    # Index the sorted BAM file
    samtools index $SORTED_BAM_FILE

    echo "Finished processing ${FASTA_FILE}."
done