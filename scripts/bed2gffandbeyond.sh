#!/bin/bash
#SBATCH --job-name=bed_to_gtf_gff2
#SBATCH --output=bed_to_gtf_gff2_%j.out
#SBATCH --error=bed_to_gtf_gff2_%j.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G


module load 'build-env/f2022'
module load 'gffread/0.12.7-gcccore-12.2.0'

# Similar to bedtogtf.sh, but this script is used to create files for antisense and annotated sense
# It also generates the fasta file for the transcripts that is needed for salmon/kallisto.
# Input 
INPUT_BED="/groups/berger/user/elin.axelsson/projects/projects_2025/pacbio/paper25/output/annos/merged_default_antisense.bed"
EXISTING_GFF="/groups/berger/user/elin.axelsson/projects/projects_2025/pacbio/paper25/modified/Schizosaccharomyces_pombe.ASM294v2.60.pc_filtered.gff3"
GENOME_FASTA="/groups/berger/user/elin.axelsson/projects/projects_2025/pacbio/paper25/annotations/Schizosaccharomyces_pombe.ASM294v2.dna.toplevel.fa"

# and output files
MODIFIED_BED="modified_merged_default_filtered.bed"
OUTPUT_TSV="isoform_bed.tsv"
OUTPUT_GFF="merged_default_filtered.gff"
MERGED_SORTED_GFF="/groups/berger/user/elin.axelsson/projects/projects_2025/pacbio/paper25/modified/sense_antisense.gff"
MERGED_SORTED_GTF="/groups/berger/user/elin.axelsson/projects/projects_2025/pacbio/paper25/modified/sense_antisense.gtf"
TRANSCRIPTS_FASTA="/groups/berger/user/elin.axelsson/projects/projects_2025/pacbio/paper25/modified/Schizosaccharomyces_pombe.ASM294v2.60.transcripts_sas.fasta"
CLEANED_TRANSCRIPTS_FASTA="/groups/berger/user/elin.axelsson/projects/projects_2025/pacbio/paper25/modified/Schizosaccharomyces_pombe.ASM294v2.60.transcripts_sas_clean.fasta"

# Process the antisense BED file to extract and split the third column
echo "Processing $INPUT_BED to extract and split the third column..."
awk -F'\t' '{split($4, parts, ";"); print parts[1] "\t" parts[2]}' $INPUT_BED > $OUTPUT_TSV
echo "Processing complete. Results saved to $OUTPUT_TSV."

# Modify the BED file to only contain the second part of column 4 after splitting on ";"
echo "Modifying $INPUT_BED to only retain the second part of column 4..."
awk -F'\t' '{split($4, parts, ";"); $4 = parts[2]; print $0}' OFS='\t' $INPUT_BED > $MODIFIED_BED
echo "Modification complete. Modified BED file saved to $MODIFIED_BED."

# Convert the modified BED file to GFF format using bed2gff
echo "Converting $MODIFIED_BED to GFF format..."
singularity run https://depot.galaxyproject.org/singularity/bed2gff:0.1.5--h9948957_1 bed2gff --bed $MODIFIED_BED --isoforms $OUTPUT_TSV --output $OUTPUT_GFF
echo "Conversion complete. GFF file saved to $OUTPUT_GFF."

# Merge the generated GFF file with the existing GFF file

echo "Merging $OUTPUT_GFF with $EXISTING_GFF..."
cat $OUTPUT_GFF $EXISTING_GFF > merged.gff
echo "Merge complete. Temporary merged file saved as merged.gff."

# Sort the merged GFF file using AGAT
echo "Sorting the merged GFF file..."
singularity run docker://quay.io/biocontainers/agat:1.4.2--pl5321hdfd78af_0 agat_convert_sp_gxf2gxf.pl --gff merged.gff -o $MERGED_SORTED_GFF
echo "Sorting complete. Final sorted GFF file saved as $MERGED_SORTED_GFF."

gffread $MERGED_SORTED_GFF -T -o $MERGED_SORTED_GTF
#singularity run docker://quay.io/biocontainers/agat:1.4.2--pl5321hdfd78af_0 agat_convert_sp_gff2gtf.pl --gff $MERGED_SORTED_GFF -o $MERGED_SORTED_GTF

sed 's/ID=transcript://g' $OUTPUT_GFF > $OUTPUT_GFF.tmp
mv $OUTPUT_GFF.tmp $OUTPUT_GFF


# Extract transcript sequences from the genomic FASTA file
echo "Extracting transcript sequences from $GENOME_FASTA using $MERGED_SORTED_GTF..."
singularity run 'https://depot.galaxyproject.org/singularity/gffread:0.12.7--hdcf5f25_4' gffread -w $TRANSCRIPTS_FASTA -g $GENOME_FASTA $MERGED_SORTED_GTF

echo "Transcript sequences have been extracted and saved to $TRANSCRIPTS_FASTA."

# Clean the headers in the TRANSCRIPTS_FASTA file
echo "Cleaning headers in $TRANSCRIPTS_FASTA..."
awk '/^>/ {print ">" substr($1, 2)} !/^>/ {print}' $TRANSCRIPTS_FASTA > $CLEANED_TRANSCRIPTS_FASTA
echo "Headers cleaned. Cleaned FASTA file saved to $CLEANED_TRANSCRIPTS_FASTA."


# Cleanup temporary merged file
rm merged.gff
rm $OUTPUT_GFF
rm $MODIFIED_BED
rm $OUTPUT_TSV
echo "Temporary files cleaned up."

echo "All steps completed successfully."