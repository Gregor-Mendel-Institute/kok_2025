#!/bin/bash
#SBATCH --job-name=tama_collapse
#SBATCH --output=tama_collapse_%j.out
#SBATCH --error=tama_collapse_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G

module load build-env/f2022
module load samtools/1.17-gcc-12.2.0

# Define paths
CONTAINER="https://depot.galaxyproject.org/singularity/gs-tama:1.0.3--hdfd78af_0"
OUTPUT_DIR="/groups/berger/user/elin.axelsson/projects/projects_2025/pacbio/paper25/output/tama_out"
#REFERENCE="/groups/berger/user/elin.axelsson/projects/projects_2025/pacbio/paper25/annotations/Schizosaccharomyces_pombe.ASM294v2.dna.toplevel.fa"
REFERENCE="/groups/berger/user/elin.axelsson/projects/projects_2025/ONS/Schizosaccharomyces_pombe_all_chromosomes.fa"

echo $REFERENCE
# List of full paths to BAM files
BAM_FILES=(
    "/scratch-cbe/users/elin.axelsson/minimap2_out_fasta/library1.transcripts_sorted.bam"
    "/scratch-cbe/users/elin.axelsson/minimap2_out_fasta/library2.transcripts_sorted.bam"
)

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Loop through BAM files and run tama_collapse.py
for BAM_FILE in "${BAM_FILES[@]}"; do
    BASENAME=$(basename $BAM_FILE .bam)
    SORTED_SAM_FILE="$OUTPUT_DIR/${BASENAME}_sorted.sam"
    OUTPUT_PREFIX="$OUTPUT_DIR/${BASENAME}"

    # Convert BAM to sorted SAM
    samtools view $BAM_FILE > $SORTED_SAM_FILE

    echo "Converted ${BAM_FILE} to sorted SAM: ${SORTED_SAM_FILE}"


    echo "Processing ${SORTED_SAM_FILE}..."

   singularity run $CONTAINER tama_collapse.py \
        -s $SORTED_SAM_FILE \
        -f $REFERENCE \
        -p $OUTPUT_PREFIX \
        -m 100 \
        -a 100 \
        -x no_cap 

    echo "Finished processing ${SORTED_SAM_FILE}."
done

singularity run $CONTAINER  tama_merge.py \
    -f /groups/berger/user/elin.axelsson/projects/projects_2025/pacbio/paper25/meta/col_trans.txt \
    -p merged_default \
    -d merge_dup 
    
  

mkdir -p $OUTPUT_DIR/merged_default
mv merged_default* $OUTPUT_DIR/merged_default
echo "Merged results saved to $OUTPUT_DIR/merged_default"