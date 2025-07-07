# load libraries
library(DESeq2)
library(tidyverse)
library(tximport)


# set parameters
kallisto_path="/groups/berger/user/elin.axelsson/projects/projects_2025/pacbio/paper25/output/kallisto_out/kallisto_counts/"
min_reads=15
min_max_reads=15
padj_th=0.1
logFC = 1
logFC_high = log2(10)


data_set = data_frame(ID=c("194804","194824","194844","194805","194825","194845"),
                      genotype=factor(c("WT","WT","WT","dko","dko","dko"),levels=c("WT","dko")),
                      rep=c("R1","R2","R3","R1","R2","R3"))
rownames(data_set) <- data_set$ID

full = paste(paste(paste(kallisto_path, data_set$ID, sep = "/"),"paired/abundance.h5",sep="/"))

# read in on transcripts level first,
transcript_counts <- tximport(full, type = "kallisto", tx2gene = NULL,txOut = T)
transcript_counts <- transcript_counts$counts
#remove with less than X reads, map antisense to gene and merge.
#transcript_counts <- transcript_counts[rowSums(transcript_counts) > min_reads,]
ids <- rownames(transcript_counts) #those are the transcript with 15 or more reads

# make tx2gene table for gene level analysis
tx2gene = data.frame(tx_id = ids, gene_id = ids)
tx2gene$gene_id = str_replace(tx2gene$tx_id, "\\.[0-9]+$", "")
tx2gene$gene_id = gsub("transcript:", "", tx2gene$gene_id)

txi <- tximport(full, type = "kallisto", tx2gene = tx2gene)


########################
# DESeq genotype
########################

dds <- DESeqDataSetFromTximport(txi, data_set, ~ genotype )
dds = dds[rowMax(counts(dds))>min_max_reads, ]
dds = DESeq(dds)

res = lfcShrink(dds,coef="genotype_dko_vs_WT", type="apeglm")
res_tib= as_tibble(res) %>% 
  mutate(gene_id = rownames(res)) %>%
  mutate(isAS = str_detect(gene_id, "^G"))

# significance, AS vs Sense, up, down
res_sig = res_tib %>%
  filter(padj < padj_th,abs(log2FoldChange)>logFC) %>% 
  arrange(desc(abs(log2FoldChange)))

# create a folder if it does not exist
output_dir = "/groups/berger/user/elin.axelsson/projects/projects_2025/pacbio/paper25/output/deseq2/"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

significant_as <- res_sig %>% filter(isAS) %>% arrange(desc(abs(log2FoldChange)))
write.table(significant_as, file = "/groups/berger/user/elin.axelsson/projects/projects_2025/pacbio/paper25/output/deseq2/significant_as_genotype_mRNA.txt", sep = "\t", quote = F, row.names = F)
significant_sense<- res_sig %>% filter(!isAS) %>% arrange(desc(abs(log2FoldChange)))
write.table(significant_sense, file = "/groups/berger/user/elin.axelsson/projects/projects_2025/pacbio/paper25/output/deseq2/significant_sense_genotype_mRNA.txt", sep = "\t", quote = F, row.names = F)

# expand
significant_as_transcripts_used <- tx2gene %>% filter(gene_id %in% significant_as$gene_id) 
write.table(significant_as_transcripts_used, file = "/groups/berger/user/elin.axelsson/projects/projects_2025/pacbio/paper25/output/deseq2/significant_as_genotype_mRNA_transcripts_used.txt", sep = "\t", quote = F, row.names = F)
des = significant_as_transcripts_used 



long_read_antisense = read_delim("/groups/berger/user/elin.axelsson/projects/projects_2025/pacbio/paper25/output/annos/merged_default_antisense.gtf",col_names = F) %>%
  filter(X3=="transcript")

colnames(long_read_antisense) = c("chrom", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

protein = read_delim("/groups/berger/user/elin.axelsson/projects/projects_2025/pacbio/paper25/modified//Schizosaccharomyces_pombe.ASM294v2.60.pc_filtered.gff3", 
                     delim = "\t",col_names = F,comment = "#")
protein = protein %>% filter(X3=="mRNA")
colnames(protein) = c("chrom", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
protein = dplyr::select(protein, chrom, start, end, strand, attribute)

# which as do not start within a protein coding gene
long_read_antisense = mutate(long_read_antisense , TSS = ifelse(strand=="+", start, end))%>%
  mutate(start=TSS,end=TSS) %>% mutate(transcript_id = str_extract( attribute,"G[0-9]+\\.[0-9]+")) %>%
  dplyr::select(chrom, start, end, strand, transcript_id)

# find overlapping PC
library(valr)
inpc = valr::bed_intersect(long_read_antisense , protein) %>% filter(strand.x!=strand.y) %>% dplyr::select(transcript_id.x) %>% distinct()

setdiff(des$tx_id, long_read_antisense$transcript_id)

long_read_antisense_de = long_read_antisense%>%filter(transcript_id%in%des$tx_id) 
inpca = valr::bed_intersect(long_read_antisense_de, protein) %>% filter(strand.x!=strand.y) %>% dplyr::select(transcript_id.x) %>% distinct()

# write three numbers to summary file
summary_file = "/groups/berger/user/elin.axelsson/projects/projects_2025/pacbio/paper25/output/deseq2/summary_genotype_mRNA_t0.txt"
write(paste0("Number of AS transcripts with significant upregulation: ", nrow(long_read_antisense_de)), summary_file)
write(paste0("Number of AS transcripts with significant upregulation that overlap PC: ", nrow(inpca)), summary_file, append = T)
