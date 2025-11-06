#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { minimap2_transcript } from './modules.nf'
include { tama_collapse } from './modules.nf'
include { tama_merge } from './modules.nf'
include { extract_gffcompare_sense_features } from './modules.nf'
include { modify_gff_pc } from './modules.nf'
include { gffcompare } from './modules.nf'
include { prepare_isoform_bed as prepare_isoform_bed_sense } from './modules.nf'
include { convert_bed_to_gff as convert_bed_to_gff_sense } from './modules.nf'
include { gffread as gffread1 } from './modules.nf'
include { modify_tama } from './modules.nf'


workflow {
    params.fastaFiles = "/groups/berger/user/elin.axelsson/projects/projects_2025/pacbio/paper25/IsoSeq_Results/library*/*.fasta"
    params.reference = "/groups/berger/user/elin.axelsson/projects/projects_2025/pacbio/paper25/annotations/Schizosaccharomyces_pombe_all_chromosomes.fa"
    params.referenceGFF = "/groups/berger/user/elin.axelsson/projects/projects_2025/pacbio/paper25/annotations/Schizosaccharomyces_pombe.ASM294v2.60.gff3"
    params.extract = "/groups/berger/user/elin.axelsson/projects/projects_2025/pacbio/paper25/annotations/Schizosaccharomyces_pombe.ASM294v2.dna.toplevel.fa"
    params.out = "./results_review"
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

    
    minimap2_transcript(fastaFiles, reference.first())
    modify_gff_pc(params.referenceGFF)

    // could change collapse parameters
    tama_collapse(minimap2_transcript.out.samFile, reference.first(),tuple(10,10))
    tama_merge(tama_collapse.out.bedFile.map { it[1] }.collect())
    modify_tama(tama_merge.out.mergedBedFile)

// Sense specific steps
    params.extra = "/groups/berger/user/elin.axelsson/projects/projects_2025/pacbio/paper25/results/annos/sense_antisense.gtf"
    gffcompare(modify_tama.out.modifiedBedFile, params.extra )
    extract_gffcompare_sense_features(gffcompare.out.trackingFile, modify_tama.out.modifiedBedFile)
    prepare_isoform_bed_sense(extract_gffcompare_sense_features.out.filteredBedFile)
    convert_bed_to_gff_sense(prepare_isoform_bed_sense.out.modifiedBedFile) 
    gffread1(convert_bed_to_gff_sense.out.gffFile)
    
}