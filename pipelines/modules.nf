#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process minimap2_transcript {
    module 'build-env/f2022:minimap2/2.24-gcccore-11.2.0:samtools/1.17-gcc-12.2.0'
    memory '32G'

    input:
    tuple val(name), path(fasta)
    path reference
    
    output:
    tuple val(name), path("${name}.sorted.bam"), path("${name}.sorted.bam.bai"), emit: bamFile
    tuple val(name), path("${name}.sorted.sam"), emit: samFile


    script:
    """
    minimap2 -ax splice:hq -uf ${reference} ${fasta} > ${name}.sam
    samtools view -S -b ${name}.sam > ${name}.bam
    samtools sort ${name}.sam -o ${name}.sorted.bam
    samtools index ${name}.sorted.bam
    samtools view ${name}.sorted.bam > ${name}.sorted.sam
    """
}
// add other settings
process tama_collapse{

    input:
    tuple val(name), path(samFile)
    path reference
    tuple val(m), val(a) 

    output:
    tuple val(name), path("*transcripts_collapsed.bed"), emit: bedFile

    script:
    """
    tama_collapse.py \
        -s ${samFile} \
        -f $reference \
        -p $name \
        -m ${m} \
        -a ${a} \
        -x no_cap 
    """
}

process tama_merge{
    publishDir "$params.out/tama", mode: 'copy'

    input:
    path(bedFile)

    output:
    path "merged_default.bed", emit: mergedBedFile

    script:
    """
    tama_merge.py \
    -f ${projectDir}/assets/col_trans.txt \
    -p merged_default \
    -d merge_dup 
    """

}
// This process filters the GFF file to retain only protein-coding genes - needed to define antisense transcripts
process modify_gff_pc {
    publishDir "$params.out/modified", mode: 'copy', pattern: "*filtered.gff*"

    input:
    path referenceGFF

    output:
    path "*pc_filtered.gff3", emit: FILTERED_GTF

    script:
    """
    BASE=\$(basename ${referenceGFF} .gff3)
    # First we need to remove all non protein coding genes from the GFF file before running gffcompare
    # This way we will find "all" antitsense transcripts (comprared to protein coding genes) in the gffcompare output
    # We will miss noverl AS to non coding genes, but that is not the focus of this analysis

    # Step 1: Extract names of ncRNA entries _ follow link!
    grep -w "ncRNA_gene" $referenceGFF | awk -F'\t' '{match(\$9, /ID=gene:([^;]+)/, arr); if (arr[1] != "") print arr[1]}' > ncRNA_names.txt
    
    # Step 2: Filter out lines containing these names in the 4th column
    grep -Fwvf ncRNA_names.txt $referenceGFF > \${BASE}_pc_filtered.gff3
    """
}


process modify_tama {
    input:
    path mergedBedFile

    output:
    path "*modified.bed", emit: modifiedBedFile

    script:
    """
    BASE=\$(basename ${mergedBedFile} .bed)
    # Modify the merged BED 
    # Change  so all G becomes sG
    sed 's/G/sG/g' ${mergedBedFile} > \${BASE}_modified.bed
    """
}


// This process runs gffcompare to compare the collapsed transcripts with the reference GTF of protein coding genes
process gffcompare {
    publishDir "$params.out/gffcompare", mode: 'copy', pattern: "*.tracking"

    input:
    path mergedBedFile
    path referenceGFF

    output:
    path "*.tracking", emit: trackingFile

    script:
    """
    # Extract the base name of the merged BED file
    BASE=\$(basename ${mergedBedFile} .bed)
    REF=\$(basename ${referenceGFF} .gff3)
    # Run gffcompare to compare the collapsed transcripts with the reference GTF
    gffcompare -R -r $referenceGFF ${mergedBedFile} -o \${BASE}_\${REF} 
    """

}

process extract_gffcompare_antisense_features {
    publishDir "$params.out/annos", mode: 'copy', pattern: "*antisense.bed"

    input:
    path trackingFile
    path mergedBedFile

    output:
    path("*antisense.bed"), emit: filteredBedFile

    script:
    """
    BASE=\$(basename ${mergedBedFile} .bed)
    # Extract IDs from the tracking file where column 4 is "X" (antisense)
    awk -F'\t' '\$4 == "x" {split(\$5, parts, "[:|]"); print parts[2]"\t"}' $trackingFile > id.txt

    # Extract rows from the BED file that match identifiers in id.txt followed by tab 
    # (to ensure we match the entire ID)    
    grep -Ff id.txt $mergedBedFile > \${BASE}_antisense.bed
    """
}


process prepare_isoform_bed {

    input:
    path modifiedBedFile

    output:
    tuple path("*.tmp.bed"), path("iso.tsv"), emit: modifiedBedFile
    script:
    """
    BASE=\$(basename ${modifiedBedFile} .bed)
    # Process the modified BED file to extract and split the third column
    awk -F'\t' '{split(\$4, parts, ";"); print parts[1] "\t" parts[2]}' ${modifiedBedFile} > iso.tsv

    # Modify the BED file to only contain the second part of column 4 after splitting on ";"
    awk -F'\t' '{split(\$4, parts, ";"); \$4 = parts[2]; print \$0}' OFS='\t' ${modifiedBedFile} > \${BASE}.tmp.bed
    """
}
 
process convert_bed_to_gff {
    publishDir "$params.out/annos", mode: 'copy'

    input:
    tuple path (modifiedBedFile), path(isoTsv)

    output:
    path "*.gff", emit: gffFile

    script:
    """
    BASE=\$(basename ${modifiedBedFile} .tmp.bed)
    # Convert the modified BED file to GFF format using bed2gff
    bed2gff --bed ${modifiedBedFile} --isoforms ${isoTsv} --output \$BASE.gff
    """
}

process gffread {
module = 'build-env/f2022:gffread/0.12.7-gcccore-12.2.0'
publishDir "$params.out/annos", mode: 'copy'

    input:
    path gffFile

    output:
    path "*.gtf", emit: gtfFile

    script:
    """
    BASE=\$(basename ${gffFile} .gff)
    # Convert the GFF file to GTF format using gffread
    gffread ${gffFile} -T -o \$BASE.gtf
    # remove transcript: and gene: from gene_id and transcript_id
    sed -i 's/transcript://g' \$BASE.gtf
    sed -i 's/gene://g' \$BASE.gtf
    """
}

process merge_gffs {
publishDir "$params.out/annos", mode: 'copy'
input:
    path outGFF
    path existingGFF

    output:
    path "sense_antisense.gff", emit: mergedGffFile

    script:
    """
    # Merge the generated GFF file with the existing GFF file and sort it using AGAT    
    cat $outGFF $existingGFF > merged.gff
    agat_convert_sp_gxf2gxf.pl --gff merged.gff -o sense_antisense.gff
    """
}

process pair_sense_as {
    module = 'build-env/f2022:r-bundle-bioconductor/3.19-foss-2023b-r-4.4.1'
    publishDir "$params.out/annos", mode: 'copy'
    input:
    path protGtf
    path asGtf

    output:
    path "sense_antisense_pairs.txt", emit: pairedOutput
    path "genes_antisense_info.txt", emit: geneInfo

    script:
    """
    Rscript ${projectDir}/bin/pair_sense_as.R ${protGtf} ${asGtf}
    """

}

process extract_fasta_seq {
    module = 'build-env/f2022:gffread/0.12.7-gcccore-12.2.0'

    input:
    path genomeFasta
    path mergedGffFile

    output:
    path "*.fasta", emit: transcriptsFasta

    script:
    """
    # Extract transcript sequences from the genomic FASTA file
    gffread -w sas.fasta -g $genomeFasta $mergedGffFile
    """
}

process kallisto_index {
   
    module = 'build-env/2020:kallisto/0.46.0-foss-2018b'

    input:
    path transcriptsFasta

    output:
    path "kallisto_index.idx", emit: kallistoIndex

    script:
    """
    # Create a kallisto index from the transcript FASTA file
    kallisto index -i kallisto_index.idx $transcriptsFasta
    """
}

process trimmomatic {
    module = 'build-env/2020:trimmomatic/0.38-java-1.8'

    input:
    tuple val(name), path(fastqFiles)

    output:
    tuple val(name), path("*_paired.fq.gz"), emit: trimmedFastq

    script:
    """
    java -jar \$EBROOTTRIMMOMATIC/trimmomatic-0.38.jar PE -threads 16\
        ${fastqFiles[0]} \
        ${fastqFiles[1]} \
        ${name}-1_paired.fq.gz \
        ${name}-1_unpaired.fq.gz \
        ${name}-2_paired.fq.gz \
        ${name}-2_unpaired.fq.gz \
        ILLUMINACLIP:${projectDir}/assets/TruSeq3-PE-2.fa:2:30:10:2:keepBothReads \
        LEADING:3 TRAILING:3 MINLEN:36 SLIDINGWINDOW:4:20
    """
      
}

process kallisto_counts {
    module = 'build-env/2020:kallisto/0.46.0-foss-2018b'
    publishDir "$params.out/kallisto", mode: 'copy'

    input:
    val(indexName) 
    path(kallistoIndex)
    tuple val(name), path(fastqFiles)

    output:
    path "${indexName}/${name}/paired", emit: kallistoOutput

    script:
    """
    mkdir -p ${indexName}/${name}/paired/
    # Run kallisto quantification on the FASTQ files using the index
    kallisto quant -i $kallistoIndex \
    -o  "${indexName}/${name}/paired/" -b 100 -t 16 --rf-stranded ${fastqFiles[0]} ${fastqFiles[1]} 
    """
}

// Processes used for supplementary analyses - not used in main workflow
process modify_gff_as {
    publishDir "$params.out/modified", mode: 'copy', pattern: "*filtered.gff*"

    input:
    path referenceGFF
    output:
    path "*as_filtered.gff3", emit: FILTERED_GTF
    script:
    """
    BASE=\$(basename ${referenceGFF} .gff3)
    # Step 1: Extract names of ncRNA entries

    grep  "antisense" $referenceGFF | awk -F'\t' '{match(\$9, /ID=gene:([^;]+)/, arr); if (arr[1] != "") print arr[1]}' > ncRNA_names.txt

    # Step 2:Keep only lines that contain the identifiers in the 4th column
    grep -Fwf ncRNA_names.txt $referenceGFF > \${BASE}_as_filtered.gff3
    """
}


process extract_gffcompare_sense_features {
    publishDir "$params.out/annos", mode: 'copy', pattern: "*sense.bed"

    input:
    path trackingFile
    path mergedBedFile

    output:
    path("*sense.bed"), emit: filteredBedFile

    script:
    """
    BASE=\$(basename ${mergedBedFile} .bed)
    # Extract IDs from the tracking file where column 4 is not "x","s","r","u" or "p"
    # use the gene name and add transcript ID

    awk -F'\t' '\$4 !~ /^[xsrup]\$/ {split(\$3,a,"|"); split(\$5,b,"|"); split(a[1],g,":"); split(b[1],q,";"); print a[1] ";" q[2]}' $trackingFile > id.txt
    #awk -F'\t' '\$4 !~ /^(x|s|r|u|p)\$/ {split(\$5, parts, "[:|]"); print parts[2]"\t"}' $trackingFile > id.txt

    awk -F'\t' 'NR==FNR {
    split(\$1,a,";");
    ids[a[2]] = \$1;    # store mapping: secondPart -> fullId
    next
    } 
    {
    split(\$4,b,";");
    if (b[2] in ids) {
        \$4 = ids[b[2]];
        OFS="\t";
        print
    }
}' id.txt $mergedBedFile > \${BASE}_sense.bed
    """
}
