# analysis of mRNA-seq wildtype mutant no stress using our sense_antisense annotation 

# load libraries
library(DESeq2)
library(tidyverse)
library(tximport)


# set parameters
kallisto_path="/groups/berger/user/elin.axelsson/projects/projects_2025/pacbio/paper25/results/kallisto/as"
path_to_sense_antisense="/groups/berger/user/elin.axelsson/projects/projects_2025/pacbio/paper25/results/annos/sense_antisense_pairs.txt"
output1 = "significant_as_mRNA_transcripts_used.txt"
output2 = "significant_as_mRNA.txt"
min_reads=1
min_max_reads=15
padj_th=0.1
logFC = 1
logFC_high = log2(10)

# read in the data

data_set = data_frame(ID=c("194804","194824","194844","194805","194825","194845"),
                  genotype=factor(c("WT","WT","WT","dko","dko","dko"),levels=c("WT","dko")),
                  rep=c("R1","R2","R3","R1","R2","R3"))
rownames(data_set) <- data_set$ID

full = paste(paste(paste(kallisto_path, data_set$ID, sep = "/"),"paired/abundance.h5",sep="/"))

# read in on transcripts level first,
transcript_counts <- tximport(full, type = "kallisto", tx2gene = NULL,txOut = T)
transcript_counts <- transcript_counts$counts
#remove with less than X reads, map antisense to gene and merge.
transcript_counts <- transcript_counts[rowSums(transcript_counts) > min_reads,]
ids <- rownames(transcript_counts) #those are the transcript with 15 or more reads

# make tx2gene table for gene level analysis
tx2gene = data.frame(tx_id = ids, gene_id = ids)
tx2gene$gene_id = str_replace(tx2gene$tx_id, "\\.[0-9]+$", "")
tx2gene$gene_id = gsub("transcript:", "", tx2gene$gene_id)

txi <- tximport(full, type = "kallisto", tx2gene = tx2gene)

sas <- sas %>% filter(antisense_transcript.anti%in%tx2gene$tx_id) %>% select(antisense_gene.anti, prot_gene.sense) %>% unique()

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

significant_as <- res_sig %>% filter(isAS) %>% arrange(desc(abs(log2FoldChange)))
# expand
significant_as_transcripts_used <- tx2gene %>% filter(gene_id %in% significant_as$gene_id) 
# write

write.table(significant_as_transcripts_used, file =output1 , sep = "\t", quote = F, row.names = F)
write.table(significant_as, file = output2, sep = "\t", quote = F, row.names = F)


