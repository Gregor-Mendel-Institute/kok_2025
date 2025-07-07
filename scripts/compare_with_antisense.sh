#!/bin/bash

# Define variables
SINGULARITY_IMAGE="https://depot.galaxyproject.org/singularity/gffcompare:0.12.6--h9948957_4"
REFERENCE_GTF="/groups/berger/user/elin.axelsson/projects/projects_2025/pacbio/paper25/annotations/Schizosaccharomyces_pombe.ASM294v2.60.gff3"
INPUT_BED="/groups/berger/user/elin.axelsson/projects/projects_2025/pacbio/paper25/output/tama_out/merged_default/merged_default.bed"

FILTERED_GTF="/groups/berger/user/elin.axelsson/projects/projects_2025/pacbio/paper25/modified/Schizosaccharomyces_pombe.ASM294v2.60.asonly.gff3"
OUTPUT_PREFIX="only_antisense"


# First we need to remove all non protein coding genes from the GFF file before running gffcompare
# This way we will find "all" antitsense transcripts (comprared to protein coding genes) in the gffcompare output
# We do not need to merge the "novel" as with the known ones
# We will miss noverl AS to non coding genes, but that is not the focus of this analysis

# Step 1: Extract names of ncRNA entries
echo "Extracting names of ncRNA entries from $REFERENCE_GTF..."
grep  "antisense" $REFERENCE_GTF | awk -F'\t' '{match($9, /ID=gene:([^;]+)/, arr); if (arr[1] != "") print arr[1]}' > ncRNA_names.txt
echo "Extracted names saved to ncRNA_names.txt."

# Step 2:Keep onlu lines that contain the identifiers in the 4th column
grep -Fwf ncRNA_names.txt $REFERENCE_GTF > $FILTERED_GTF
echo "Filtered GFF3 file saved to $FILTERED_GTF."

# Step 3: Run gffcompare
echo "Running gffcompare..."
singularity run $SINGULARITY_IMAGE gffcompare -R -r $FILTERED_GTF $INPUT_BED -o $OUTPUT_PREFIX
echo "gffcompare has finished. Results are saved with prefix: $OUTPUT_PREFIX."

# Step 4: Clean up
echo "Cleaning up temporary files..."
#rm ncRNA_names.txt

# Step 5: Mv output files to the correct folder
mv $OUTPUT_PREFIX.* /groups/berger/user/elin.axelsson/projects/projects_2025/pacbio/paper25/output/gffcompare/
echo "Output files moved to ../output/gffcompare/."
echo "All steps completed successfully."

# for the file ../output/gffcompare/only_antisense.tracking, use column 4 and count the number of each unique entry 

awk -F'\t' '{print $4}' /groups/berger/user/elin.axelsson/projects/projects_2025/pacbio/paper25/output/gffcompare/only_antisense.tracking | sort | uniq -c | sort -nr > /groups/berger/user/elin.axelsson/projects/projects_2025/pacbio/paper25/output/gffcompare/only_antisense_counts.txt

# Filter ../output/gffcompare/only_antisense.tracking for rows there the 4th is not x, u,s, eller p and count unique entries in the 3rd column
awk -F'\t' '$4 != "x" && $4 != "u" && $4 != "s" && $4 != "p" {print $3}'  /groups/berger/user/elin.axelsson/projects/projects_2025/pacbio/paper25/output/gffcompare/only_antisense.tracking | sort | uniq -c | sort -nr > /groups/berger/user/elin.axelsson/projects/projects_2025/pacbio/paper25/output/gffcompare/only_antisense_filtered_counts.txt

# print summary of the counts first columns of only_antisense_filtered_counts.txt
echo "Summary of counts in only_antisense_filtered_counts.txt:"
awk '{sum += $1} END {print "Total counts:", sum}' /groups/berger/user/elin.axelsson/projects/projects_2025/pacbio/paper25/output/gffcompare/only_antisense_filtered_counts.txt
q