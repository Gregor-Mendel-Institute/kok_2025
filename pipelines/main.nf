#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { minimap2_transcript } from './modules.nf'
include { tama_collapse } from './modules.nf'
include { tama_merge } from './modules.nf'
include { modify_gff_pc } from './modules.nf'
include { gffcompare } from './modules.nf'
include { extract_gffcompare_antisense_features } from './modules.nf'
include { prepare_isoform_bed } from './modules.nf'
include { convert_bed_to_gff } from './modules.nf'
include { gffread } from './modules.nf'
include { gffread as gffread1 } from './modules.nf'
include { merge_gffs } from './modules.nf'
include { extract_fasta_seq } from './modules.nf'
include { kallisto_index } from './modules.nf'
include { trimmomatic } from './modules.nf'
include { kallisto_counts } from './modules.nf'
include { modify_tama } from './modules.nf'
include { pair_sense_as } from './modules.nf'

workflow {
    params.fastaFiles = "/groups/berger/user/elin.axelsson/projects/projects_2025/pacbio/paper25/IsoSeq_Results/library*/*.fasta"
    params.reference = "/groups/berger/user/elin.axelsson/projects/projects_2025/pacbio/paper25/annotations/Schizosaccharomyces_pombe_all_chromosomes.fa"
    params.referenceGFF = "/groups/berger/user/elin.axelsson/projects/projects_2025/pacbio/paper25/annotations/Schizosaccharomyces_pombe.ASM294v2.60.gff3"
    params.extract = "/groups/berger/user/elin.axelsson/projects/projects_2025/pacbio/paper25/annotations/Schizosaccharomyces_pombe.ASM294v2.dna.toplevel.fa"
    params.files = '/scratch-cbe/users/elin.axelsson/stressfiles/*R{1,2}_001.fastq.gz'
    params.out = "./results"
    Channel
        .fromPath(params.fastaFiles)
        .map { file -> tuple(file.baseName, file) }
        .set { fastaFiles }
        
    Channel
        .fromPath(params.reference)
        .set { reference }

    Channel
        .fromPath(params.extract)
        .set { extract }
        // this is because of the mistake in choice of reference

    Channel
        .fromFilePairs(params.files)
        .set { files }

    // remove the _SXX from the file names
    files = files.map { k, v -> tuple(k.replaceAll(/_S\d+$/, ''), v) }
    trimmomatic(files)
    minimap2_transcript(fastaFiles, reference.first())
    modify_gff_pc(params.referenceGFF)

    // could change collapse parameters
    tama_collapse(minimap2_transcript.out.samFile, reference.first(),tuple(100,100))
    tama_merge(tama_collapse.out.bedFile.map { it[1] }.collect())
    modify_tama(tama_merge.out.mergedBedFile)
    gffcompare(tama_merge.out.mergedBedFile, modify_gff_pc.out.FILTERED_GTF)
    
    // Antisense 
    extract_gffcompare_antisense_features(gffcompare.out.trackingFile, tama_merge.out.mergedBedFile)
    prepare_isoform_bed(extract_gffcompare_antisense_features.out.filteredBedFile)
    convert_bed_to_gff(prepare_isoform_bed.out.modifiedBedFile) 
    gffread1(convert_bed_to_gff.out.gffFile)
    merge_gffs(convert_bed_to_gff.out.gffFile, modify_gff_pc.out.FILTERED_GTF)
    gffread(merge_gffs.out.mergedGffFile)
    //run pair
    pair_sense_as( modify_gff_pc.out.FILTERED_GTF, convert_bed_to_gff.out.gffFile )

    // Kallisto quantification
    extract_fasta_seq(extract, gffread.out.gtfFile )
    kallisto_index(extract_fasta_seq.out.transcriptsFasta)
    kallisto_counts("as", kallisto_index.out.kallistoIndex.first(), trimmomatic.out.trimmedFastq)
    

}
