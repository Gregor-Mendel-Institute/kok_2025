#!/bin/bash
#SBATCH --job-name=gffcompare
#SBATCH --output=gffcompare_%j.out
#SBATCH --error=gffcompare_%j.err
#SBATCH --time=00:30:00
#SBATCH --mem=1G


# Define variables
SINGULARITY_IMAGE="https://depot.galaxyproject.org/singularity/gffcompare:0.12.6--h9948957_4"
REFERENCE_GTF="/groups/berger/user/elin.axelsson/projects/projects_2025/pacbio/paper25/annotations/Schizosaccharomyces_pombe.ASM294v2.60.gff3"
INPUT_BED="/groups/berger/user/elin.axelsson/projects/projects_2025/pacbio/paper25/output/tama_out/merged_default/merged_default.bed"

mkdir -p /groups/berger/user/elin.axelsson/projects/projects_2025/pacbio/paper25/modified/
FILTERED_GTF="/groups/berger/user/elin.axelsson/projects/projects_2025/pacbio/paper25/modified/Schizosaccharomyces_pombe.ASM294v2.60.pc_filtered.gff3"
OUTPUT_PREFIX="merged_default_w_anno"
OUTPUT_DIR="/groups/berger/user/elin.axelsson/projects/projects_2025/pacbio/paper25/output/gffcompare"
mkdir -p $OUTPUT_DIR


# First we need to remove all non protein coding genes from the GFF file before running gffcompare
# This way we will find "all" antitsense transcripts (comprared to protein coding genes) in the gffcompare output
# We do not need to merge the "novel" as with the known ones
# We will miss noverl AS to non coding genes, but that is not the focus of this analysis

# Step 1: Extract names of ncRNA entries
echo "Extracting names of ncRNA entries from $REFERENCE_GTF..."
grep -w "ncRNA_gene" $REFERENCE_GTF | awk -F'\t' '{match($9, /ID=gene:([^;]+)/, arr); if (arr[1] != "") print arr[1]}' > ncRNA_names.txt
echo "Extracted names saved to ncRNA_names.txt."

# Step 2: Filter out lines containing these names in the 4th column
echo "Filtering out lines with ncRNA names from $REFERENCE_GTF..."
grep -Fwvf ncRNA_names.txt $REFERENCE_GTF > $FILTERED_GTF
echo "Filtered GFF3 file saved to $FILTERED_GTF."

# Step 3: Run gffcompare
echo "Running gffcompare..."
singularity run $SINGULARITY_IMAGE gffcompare -R -r $FILTERED_GTF $INPUT_BED -o $OUTPUT_PREFIX
echo "gffcompare has finished. Results are saved with prefix: $OUTPUT_PREFIX."

# Step 4: Clean up
echo "Cleaning up temporary files..."
rm ncRNA_names.txt

# Step 5: Mv output files to the correct folder
mv $OUTPUT_PREFIX.* $OUTPUT_DIR
echo "Output files moved to $OUTPUT_DIR."
echo "All steps completed successfully."

