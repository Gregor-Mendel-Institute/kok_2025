#############################
## script to filter sense transcripts based on intron length and TSS location
## author: Elin Axelsson-Ekker
#############################

## input:
## 1, path to sense gtf from pipeline
#gtf_sense_path <- "/path_to_long_read_pipeline/results_review/annos/merged_default_modified_sense.gtf"

## 2, path to official annotation gff3
#gff_anno_gff3_path <- "/path_to_long_read_pipeline/annotations/Schizosaccharomyces_pombe.ASM294v2.60.gff3"

## 3, path to output file
#output_path <- "/path_to_long_read_pipeline/results_review/annos/merged_default_modified_sense_filtered_TSS_fixed.gtf"

library(tidyverse)
library(valr)
# more?
##############################
# read in sense gtf from pipeline
gtf_sense <- read_delim(gtf_sense_path,delim = "\t",comment = "#",col_names = FALSE)
as_was_gtf_sense <- gtf_sense # keep to get right format on output
# official annotation gff3
gff_anno <- read_delim(gff_anno_gff3_path,delim = "\t",comment = "#",col_names = FALSE)

##############################
# extract transcript name from gtf_sense
gtf_sense <- gtf_sense %>% mutate(ind=1:nrow(gtf_sense),transcript_id = 
                              str_extract(X9, "transcript_id \\\"sG[0-9]+.[0-9]+") %>% gsub("transcript_id ","",.)) %>%
  mutate(transcript_id = str_remove_all(transcript_id,"\""))

# remove some extra ;s at the end of X9
gtf_sense <- gtf_sense %>% mutate(X9 = str_remove(X9,";$"))

##############################
## Intron length filtering
##############################
# group ny transcript id and feature type (X9,X3). # there is only one transcript per X9!
gtf_sense %>% group_by(X9,X3) %>%
  summarise(lengths = sum(X5-X4)) %>% # length of per feature type (per transcript)
  filter(X3%in%c("transcript","exon")) %>%  # keep transcript and exon sum
  pivot_wider(names_from = X3, values_from = lengths) %>%
  mutate(intron_length=transcript-exon)%>% arrange(desc(intron_length)) -> intron_lengths

# plot intron lengths
intron_lengths %>% 
  ggplot(aes(x=log10(intron_length))) + geom_histogram(bins=100)

# keep when itrons lengths < 100 bp
intron_lengths_keep <- intron_lengths%>% filter(intron_length<100) 
gtf_sense_filtered <- gtf_sense %>% filter((X9 %in% intron_lengths_keep$X9))

gtf_sense_filtered <- gtf_sense_filtered %>% filter(X3=="transcript") %>%
  # gene are part of X9 between "gene_id \" and last \"
  mutate(gene = str_extract(X9, "gene_id [^;]+") %>% gsub("gene_id ","",.)) %>%
  mutate(gene = str_remove_all(gene,"\""))%>%
  mutate(transcript_id = str_extract(X9, "transcript_id \\\"sG[0-9]+.[0-9]+") %>% gsub("transcript_id ","",.)) %>%
  mutate(transcript_id = str_remove_all(transcript_id,"\"")) %>%
  dplyr::select(gene,transcript_id,X1,X4,X5,X7) %>% rename(chrom=X1,start=X4,end=X5,strand=X7)

##############################
## TSS filtering
##############################
# TSS can only be in gene body or upstreams if not overlapping (either strand) another gene - we use official annotation for this
# process official gff3 annotation
gff_anno <- gff_anno %>% filter(X3=="mRNA") %>% mutate(gene=str_extract(X9, "ID=transcript:[^;]+") %>% gsub("Name=","",.)) %>%
  dplyr::select(gene,X1,X4,X5,X7) %>% rename(chrom=X1,start=X4,end=X5,strand=X7) %>%
  mutate(gene = str_remove(gene,"ID=transcript:")) %>%
  mutate(gene = str_remove(gene,"\\.1$"))

# define TSS in our sense data
tss_sense_gtf <- gtf_sense_filtered %>% mutate(tss=ifelse(strand=="+",start,end)) %>%
  dplyr::select(chrom,tss,strand,gene,transcript_id) %>%
  mutate(start=tss-1,end=tss) %>% dplyr::select(-tss)

# find overlaps
bed_intersect(tss_sense_gtf,gff_anno) %>% filter(gene.x!=gene.y)

# Case1: If TSS of our sense overlap the same official gene as the one assigned, we will keep it
case1 = bed_intersect(tss_sense_gtf,gff_anno) %>% 
  filter(gene.x==gene.y)%>% mutate(X9 = paste0("transcript_id \"",transcript_id.x,"\"; gene_id \"",gene.x,"\"")) %>%
  dplyr::select(transcript_id.x,X9)

# remove case1 from tss_sense_gtf
gtf_TSS_left = tss_sense_gtf %>% filter(!transcript_id %in% case1$transcript_id.x) %>% filter(!str_detect(gene,"^G")) 

# Case2: Different genes in assignment of official, use official 
case2 = bed_intersect(gtf_TSS_left,gff_anno) %>% filter(strand.x==strand.y) %>%
  mutate(X9 = paste0("transcript_id \"",transcript_id.x,"\"; gene_id \"",gene.y,"\"")) %>%
  dplyr::select(transcript_id.x,X9)

# remove case2 from tss_sense_gtf
gtf_TSS_left2 = gtf_TSS_left %>% filter(!transcript_id %in% case2$transcript_id.x)

# Case3, tss in intergenic regions, those are kept for now
case3 = bed_intersect(gtf_TSS_left2,gff_anno,invert=T) %>% mutate(X9 = paste0("transcript_id \"",transcript_id,"\"; gene_id \"",gene,"\"")) %>%
  dplyr::select(transcript_id,X9) %>% dplyr::rename(transcript_id.x=transcript_id)

# remove case3 from tss_sense_gtf
# gtf_TSS_left3 = gtf_TSS_left2 %>% filter(!transcript_id %in% case3$transcript_id.x)

keep = rbind(case1,case2,case3)

## extract the cases from original gtf_sense
right_join(gtf_sense,keep,by=c("transcript_id"="transcript_id.x"),relationship = "many-to-many") %>%
  dplyr::select(-X9.x,-ind,-transcript_id) %>% rename(X9=X9.y) %>% dplyr::select(colnames(as_was_gtf_sense)) -> gtf_change

##############################
# write output
##############################
write_delim(gtf_change,output_path,
            col_names = FALSE, delim = "\t",quote = "none", escape = "none")
