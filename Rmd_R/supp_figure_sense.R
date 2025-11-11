#############################
## script to compare correlation of read profiles between WT and hrp1hrp3 (sense vs antisense)
## also create plot S11A and S11B
## author: Elin Axelsson-Ekker
#############################
library(rtracklayer)
library(GenomicRanges)
library(valr)
library(FSA)
library(rstatix)
library(ggsignif)
library(tidyverse)


## input:
# merged bw files
#bw1 <- "/path to bigwigs/sense_mRNA/merged_bw/WT.merged.forward.bw"
#bw2 <- "/path to bigwigs/sense_mRNA/merged_bw/WT.merged.reverse.bw"
#bw3 <- "/path to bigwigs/sense_mRNA/merged_bw/hrp1hrp3.merged.forward.bw"
#bw4 <- "/path to bigwigs/sense_mRNA/merged_bw/hrp1hrp3.merged.reverse.bw"

# annotations
#anti_sense_path <- "/path to long_read_pipeline/results/annos/merged_default_antisense.gtf"
#sense_path <- "/path to long_read_pipeline/results_review/annos/merged_default_modified_sense_filtered_TSS_fixed.gtf"

# output
#output_path <- "correlation_sense_antisense_wt_hrp13.pdf"
#output_path1 <- "PRO_seq_sense.pdf"

#bw_path <- "path to bw files/merged_normed/"
#gtf_file <- sense_path

# -----------------------------------------------
# Figure S11B
# -----------------------------------------------

## function that calculates binned coverage from bigwig file
calculate_binned_coverage_bw <- function(bw_file, chromosomes =c("I","II","III"), bin_size = 10) {
  
  # Import BigWig header to get chromosome lengths
  cat("Reading BigWig header from file")
  bw_info <- seqinfo(BigWigFile(bw_file))
  chr_length <- seqlengths(bw_info)[chromosomes]
  
  
  # Create bins
  bins <- list()
  for( i in 1:length(chromosomes)){
    bin_starts <- seq(1, chr_length[chromosomes[i]], by = bin_size)
    bin_ends   <- pmin(bin_starts + bin_size - 1, chr_length[chromosomes[i]])
    bins[[i]] <- GRanges(seqnames = chromosomes[i], ranges = IRanges(bin_starts, bin_ends))
  }
  # Combine bins into a single GRanges object
  bins <- do.call(c, bins)
  
  # Import BigWig coverage for just those bins
  cat("Importing coverage from BigWig 1...\n")
  cov <- import(bw_file, which = bins, as = "RleList")
  seqlevels(bins) <- names(cov)
  
  # Summarize coverage per bin
  bin_cov <- binnedAverage(bins, cov, "mean")$mean
  
  results <- data.frame(
    chrom    = seqnames(bins),
    start     = start(bins),
    end       = end(bins),
    cov      = bin_cov
  )
  
  return(results)
}


######################################
## calculate binned coverage for each bigwig
######################################

res = list()
for ( i in c(bw1,bw2,bw3,bw4) ) {
  if ( !file.exists(i) ) {
    stop("File does not exist: ", i)
  }
  res[[i]] = calculate_binned_coverage_bw(i,bin_size=50)
}

# name apporopriately
names(res) = c("WT_forward","WT_reverse","hrp13_forward","hrp13_reverse")

######################################
# read in annotations
######################################

anti_sense <- read_delim(anti_sense_path,col_names = FALSE, 
                  delim = "\t",comment = "#")
# keep transcripts only 
anti_sense <- anti_sense %>% filter(X3=="transcript") %>% mutate(gene=str_extract(X9, "gene_id \"[^\"]+\"") %>% gsub("gene_id \"","",.)) %>% mutate(tx=str_extract(X9, "transcript_id \"[^\"]+\"") %>% gsub("transcript_id \"","",.)) %>%
  dplyr::select(gene,tx,X1,X4,X5,X7) %>% rename(chrom=X1,start=X4,end=X5,strand=X7) %>% mutate(gene = str_remove(gene,"\""),
                                                                                             tx = str_remove(tx,"\""))
sense <-  read_delim(sense_path,col_names = FALSE, 
 delim = "\t",comment = "#")

# keep transcripts only
sense <- sense %>% filter(X3=="transcript") %>% mutate(gene=str_extract(X9, "gene_id \"[^\"]+\"") %>% gsub("gene_id \"","",.)) %>% mutate(tx=str_extract(X9, "transcript_id \"[^\"]+\"") %>% gsub("transcript_id \"","",.)) %>%
  dplyr::select(gene,tx,X1,X4,X5,X7) %>% rename(chrom=X1,start=X4,end=X5,strand=X7) %>% mutate(gene = str_remove(gene,"\""),
                                                                                               tx = str_remove(tx,"\""))

######################################
# intersect coverage with annotations
######################################

overlap_function <- function(cov_data, anno_data,strand){
  # intersect using valr
  intersected <- valr::bed_intersect(cov_data,anno_data) %>%
    filter(strand.y==strand) %>%
    mutate(start = start.x, end = end.x,strand = strand.y,cov= cov.x) %>%
    dplyr::select(chrom, start, end, cov, gene=tx.y, strand, anno_start=start.y, anno_end=end.y)
  return(intersected)
}

wt_sense_forward <-overlap_function(res[["WT_forward"]],sense,"+")
wt_sense_reverse <- overlap_function(res[["WT_reverse"]],sense,"-")
wt_sense <- rbind(wt_sense_forward,wt_sense_reverse)

hrp13_sense_forward <- overlap_function(res[["hrp13_forward"]],sense,"+")
hrp13_sense_reverse <- overlap_function(res[["hrp13_reverse"]],sense,"-")
hrp13_sense <- rbind(hrp13_sense_forward,hrp13_sense_reverse)

wt_as_forward <- overlap_function(res[["WT_forward"]],anti_sense,"+")
wt_as_reverse <- overlap_function(res[["WT_reverse"]],anti_sense,"-")
wt_as <- rbind(wt_as_forward,wt_as_reverse)

hrp13_as_forward <- overlap_function(res[["hrp13_forward"]],anti_sense,"+")
hrp13_as_reverse <- overlap_function(res[["hrp13_reverse"]],anti_sense,"-")
hrp13_as <- rbind(hrp13_as_forward,hrp13_as_reverse)


######################################
# correlation
######################################
# combine mutant and WT data per annotation
as_data = rbind(wt_as %>% mutate(genotype="WT"), hrp13_as %>% mutate(genotype="hrp13"))
sense_data = rbind(wt_sense %>% mutate(genotype="WT"), hrp13_sense %>% mutate(genotype="hrp13"))


as_data %>% pivot_wider(names_from= genotype,values_from=cov) %>%
  group_by(gene) %>%
 mutate(sum_WT = sum(WT), m_WT = mean(WT), sum_hrp13 = sum(hrp13), m_hrp13 = mean(hrp13)) %>%
  mutate(prop_WT = WT/sum_WT, prop_hrp13 = hrp13/sum_hrp13) %>%
  ungroup() -> tot_as

# 
sense_data %>% distinct() %>% pivot_wider(names_from= genotype,values_from=cov) %>%
  group_by(gene) %>%
  mutate(sum_WT = sum(WT), m_WT = mean(WT), sum_hrp13 = sum(hrp13), m_hrp13 = mean(hrp13)) %>%
  mutate(prop_WT = WT/sum_WT, prop_hrp13 = hrp13/sum_hrp13) %>%
  ungroup() -> tot_sense

cor_as <- tot_as %>% 
  group_by(gene) %>%
  summarize(cor= cor(prop_WT,prop_hrp13,method="pearson")) %>% ungroup()

cor_sense <- tot_sense %>% 
  group_by(gene) %>%
  summarize(cor= cor(prop_WT,prop_hrp13,method="pearson")) %>% ungroup()

tot = rbind(cor_as %>% mutate(type="antisense"), cor_sense %>% mutate(type="sense")) 
w_results <-  wilcox_test(cor ~ type,data = tot, p.adjust.method = "holm")
stopifnot(w_results$p < 0.001) # significant difference
# change order of type 
tot$type = factor(tot$type, levels=c("sense","antisense"))

pdf(output_path,width=5,height=5)
p = tot %>%
  ggplot(aes(x=type,y=cor,fill=type)) + 
  scale_fill_manual(values=c("sense"="#1f77b4","antisense"="#d4a017")) +
  theme_classic() + theme(legend.position = "none") +
  geom_boxplot(outliers = F) +
  ylim(-1,1.2) +
  # add n per group
  geom_text(data = tot %>% group_by(type) %>% summarize(n = sum(!is.na(cor))) %>% ungroup(),
            aes(x = type, y = 1.1, label = paste0("n=", n)), vjust = -0.5) +
  ylab("Pearson correlation") + xlab("") + ggtitle("Correlation of read profiles between WT and hrp13") +
  ## Add stars for significance,  *** between all three groups =
  geom_signif(comparisons = list(c("sense", "antisense")),
              map_signif_level = TRUE,
              annotations = c("***"),vjust=0.5, y_position = c(.95,.95), tip_length = 0.03)
ggsave(output_path,width=5,height=5)
tot %>% group_by(type) %>% summarize(med = median(cor,na.rm=TRUE), n = sum(!is.na(cor)), missed = sum(is.na(cor))) -> summary_stats
wilcox_effsize(cor ~ type, data = tot)

# -----------------------------------------------
# Figures S11A
# -----------------------------------------------

library(BRGenomics)
library(tidyverse)
library(GenomicAlignments)
library(GenomicFeatures)
library(ggthemes)

## input


# Function to normalize GRanges object to RPM

RPMnorm <- function(gr) {
  gr$score <- (1e6 / sum(gr$score)) * gr$score
  return(gr)
}
totReads <- function(gr) {
  return(sum(gr$score))
}

# load merged bigWig files
PROseq_merged.lst <- list(
  "WT" = import_bigWig(file.path(bw_path, "WT_PRO_fwd.bw"), file.path(bw_path, "WT_PRO_rev.bw")),
  "hrp1" = import_bigWig(file.path(bw_path, "hrp1_PRO_fwd.bw"), file.path(bw_path, "hrp1_PRO_rev.bw")),
  "hrp3" = import_bigWig(file.path(bw_path, "hrp3_PRO_fwd.bw"), file.path(bw_path, "hrp3_PRO_rev.bw")),
  "hrp13" = import_bigWig(file.path(bw_path, "hrp13_PRO_fwd.bw"), file.path(bw_path, "hrp13_PRO_rev.bw"))
)


# Normalize merged data
PROseq_merged_normed.lst <- lapply(PROseq_merged.lst, RPMnorm)

lapply(PROseq_merged.lst, function(x) quantile(x$score))

# Load gene annotations
txdb <- makeTxDbFromGFF(gtf_file, organism = "Schizosaccharomyces pombe") 

# Get list of transcripts
genes.gr <- tidyChromosomes(transcripts(txdb), keep.X = FALSE, keep.Y = FALSE)

# Remove transcripts smaller than 600bp
genes.gr <- genes.gr[width(genes.gr) > 600]

# Define promoter regions
PR.gr <- genebodies(genes.gr, -200, 200, fix.end = "start")

# Subsample promoter regions for forward strand
PR_sub_fwd.mat <- metaSubsample(PROseq_merged_normed.lst, PR.gr, binsize = 1, first.output.xval = -200)

PR_counts = getCountsByPositions(PROseq_merged_normed.lst, PR.gr, binsize = 10)


PR_sub_fwd.mat %>%
  # refactor so that the order is WT, hrp1, hrp3, hrp13
  mutate(sample.name = factor(sample.name, levels = c("WT", "hrp1", "hrp3", "hrp13"))) %>%
  ggplot( aes(x = x, y = mean)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = sample.name), alpha = 0.5) +
  geom_line(aes(color = sample.name), size = 0.5) +
  facet_grid(sample.name ~ .) +
  scale_color_manual(values = c(WT = "darkgrey",hrp1 = "#0072B2", hrp3 = "#CC79A7",  hrp13 = "#E69F00")) +
  scale_fill_manual(values = c(WT = "darkgrey", hrp1 = "#0072B2", hrp3 = "#CC79A7", hrp13 = "#E69F00")) +
  scale_x_continuous(breaks = seq(-200, 200, by = 200), limits = c(-200, 200)) +
  scale_y_continuous(limits = c(0, 1)) + 
  theme_classic() +
  labs(x = "Distance from TSS (bp)", y = "Mean + 75% CI") +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )
ggsave(output_path1, width = 6, height = 6)
    
