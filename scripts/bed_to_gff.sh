#!/bin/bash
#SBATCH --job-name=bed_to_gtf_gff
#SBATCH --output=bed_to_gtf_gff_%j.out
#SBATCH --error=bed_to_gtf_gff_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G

# Load required modules
module load build-env/.f2021
module load build-env/f2021
module load kent_tools/411-gcc-10.2.0
module load bedtools/2.30.0-gcc-10.2.0

## create bed file with only antisense transcription
## create gff and gtf files from the bed file


# Define input 
BED_FILE="/groups/berger/user/elin.axelsson/projects/projects_2025/pacbio/paper25/output/tama_out/merged_default/merged_default.bed"  # Folder containing .bed file from tama
TRACKING_FILE="/groups/berger/user/elin.axelsson/projects/projects_2025/pacbio/paper25/output/gffcompare/merged_default_w_anno.tracking"  # File containing identifiers to filter BED rows (x=antisense
BASENAME=$(basename $BED_FILE .bed)

# Output directories and files
OUTPUT_DIR="/groups/berger/user/elin.axelsson/projects/projects_2025/pacbio/paper25/output/annos"  # Folder to store .gtf files
FILTERED_BED="$OUTPUT_DIR/${BASENAME}_antisense.bed"
ID_FILE="$OUTPUT_DIR/filtered_ids.txt"
GTF_FILE="$OUTPUT_DIR/${BASENAME}_antisense.gtf"
GFF_FILE="$OUTPUT_DIR/${BASENAME}_antisense.gff"
OUTPUT_TSV="isoform_bed.tsv"
MODIFIED_BED="$OUTPUT_DIR/${BASENAME}_modified.bed"

# Create output directories if they don't exist
mkdir -p $OUTPUT_DIR

echo "Processing tracking file: $ID_FILE..."
echo "Using $TRACKING_FILE to extract IDs with X (antisense) from column 4..."
# Extract IDs from the tracking file where column 4 is "X" (antisense)
awk -F'\t' '$4 == "x" {split($5, parts, "[:|]"); print parts[2]"\t"}' $TRACKING_FILE > $ID_FILE


echo "Filtering $BED_FILE based on identifiers in $ID_FILE..."
# Extract rows from the BED file that match identifiers in id.txt followed by tab 
# (to ensure we match the entire ID)    
grep -Ff $ID_FILE $BED_FILE > $FILTERED_BED


# Process the modified BED file to extract and split the third column
echo "Processing $FILTERED_BED to extract and split the third column..."
awk -F'\t' '{split($4, parts, ";"); print parts[1] "\t" parts[2]}' $FILTERED_BED > $OUTPUT_TSV
echo "Processing complete. Results saved to $OUTPUT_TSV."

# Modify the BED file to only contain the second part of column 4 after splitting on ";"
echo "Modifying $FILTERED_BED to only retain the second part of column 4..."
awk -F'\t' '{split($4, parts, ";"); $4 = parts[2]; print $0}' OFS='\t' $FILTERED_BED > $MODIFIED_BED
echo "Modification complete. Modified BED file saved to $MODIFIED_BED."

# Convert the modified BED file to GFF format using bed2gff
echo "Converting $MODIFIED_BED to GFF format..."
singularity run https://depot.galaxyproject.org/singularity/bed2gff:0.1.5--h9948957_1 bed2gff --bed $MODIFIED_BED --isoforms $OUTPUT_TSV --output $GFF_FILE
echo "Conversion complete. GFF file saved to $GFF_FILE."


module load 'build-env/f2022'
module load 'gffread/0.12.7-gcccore-12.2.0'

# Merge the generated GFF file with the existing GFF file
gffread $GFF_FILE -T -o $GTF_FILE

# remove the modified BED file
rm $MODIFIED_BED
rm $OUTPUT_TSV

  
