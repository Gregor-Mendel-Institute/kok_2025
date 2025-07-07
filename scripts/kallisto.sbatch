#!/bin/bash
#
#SBATCH --job-name=JY_mRNA_kallisto
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=0-5
#SBATCH --mem-per-cpu=64G
#SBATCH --mail-type=BEGIN,END,FAIL



#1. set directory for depositing intermediate files
work_folder=/scratch-cbe/users/elin.axelsson/nf_stress
mkdir -p ${work_folder}

#3. Sequencing file names base location
fasta=/groups/berger/user/elin.axelsson/projects/projects_2025/pacbio/paper25/modified/Schizosaccharomyces_pombe.ASM294v2.60.transcripts_sas.fasta
index=/groups/berger/user/elin.axelsson/projects/projects_2025/pacbio/paper25/modified/sas_index
bamsdir=/scratch-cbe/users/elin.axelsson/stressfiles/
#bams=( 194804 194805 194824 194825 194844 194845 194812 194832 194852 194813 194833 194853)
bams=( 194804 194824 194844 194805 194825 194845 )
#bams=( 194808 194828 194848 194809 194829 194849 194816 194836 194856 194817 194837 194857 194820 194840 194860 194821 194841 194861 )
#bams=( 194857 194841 194848 194849 )
#4. Output basename
output=${bams[$SLURM_ARRAY_TASK_ID]}
#5. Non-transient Output directory
odir=/groups/berger/user/elin.axelsson/projects/projects_2025/pacbio/paper25/output/kallisto_out/
mkdir -p $odir


ml build-env/2020
ml fastqc/0.11.8-java-1.8
ml samtools/1.9-foss-2018b
ml bedtools/2.27.1-foss-2018b
ml trimmomatic/0.38-java-1.8
ml trimmomatic/0.38-java-1.8
ml kallisto/0.46.0-foss-2018b
#ml star/2.7.1a-foss-2018b


Work_BAM1=`ls $bamsdir/${bams[$SLURM_ARRAY_TASK_ID]}_S* | egrep 'S[0-9]+_R1_.*\.fastq.gz'`
Work_BAM2=`ls $bamsdir/${bams[$SLURM_ARRAY_TASK_ID]}_S* | egrep 'S[0-9]+_R2_.*\.fastq.gz'`


set -o pipefail


#mkdir -p `echo $odir"fastQC/"$output`
fastqc -t 2 -o `echo $odir"fastQC/"$output` $Work_BAM1 $Work_BAM2

#Trim adaptors
mkdir -p `echo $work_folder"/AAAVVJTM5_0_R13470_20220422/demultiplexed/adaptor_trimmed_fqs/"`

#Note: here using the Nextera sequences. Also possible to do this with arbitrary/TruSeq adaptors.
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.38.jar PE -threads 16\
      `echo $Work_BAM1` \
      `echo $Work_BAM2` \
      `echo $work_folder"/AAAVVJTM5_0_R13470_20220422/demultiplexed/adaptor_trimmed_fqs/"$output"-1_paired.fq.gz"` \
      `echo $work_folder"/AAAVVJTM5_0_R13470_20220422/demultiplexed/adaptor_trimmed_fqs/"$output"-1_unpaired.fq.gz"` \
      `echo $work_folder"/AAAVVJTM5_0_R13470_20220422/demultiplexed/adaptor_trimmed_fqs/"$output"-2_paired.fq.gz"` \
      `echo $work_folder"/AAAVVJTM5_0_R13470_20220422/demultiplexed/adaptor_trimmed_fqs/"$output"-2_unpaired.fq.gz"` \
      ILLUMINACLIP:/groups/berger/user/elin.axelsson/projects/projects_2025/pacbio/paper25/meta/TruSeq3-PE-2.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36 \
      SLIDINGWINDOW:4:20

#Get counts per transcript from paired end reads.
mkdir -p `echo $odir"kallisto_counts/"$output"/paired"`
mkdir -p `echo $odir"kallisto_counts/"$output"/unpaired_1"`
mkdir -p `echo $odir"kallisto_counts/"$output"/unpaired_2"`

#This indexing step has already been done and does not need to be repeated. The input fasta is both the standard CDS (exons), as well as their reverse complement to try to capture some antisense transcription.
# crate index if not already done
if [ ! -f $index ]; then
    kallisto index -i $index $fasta 
else
    echo "Index already exists, skipping indexing step."
fi

kallisto quant -i $index -o `echo $odir"kallisto_counts/"$output"/paired"` -b 100 -t 16 --rf-stranded `echo $work_folder"/AAAVVJTM5_0_R13470_20220422/demultiplexed/adaptor_trimmed_fqs/"$output"-1_paired.fq.gz"` `echo $work_folder"/AAAVVJTM5_0_R13470_20220422/demultiplexed/adaptor_trimmed_fqs/"$output"-2_paired.fq.gz"`


