library(tidyverse)

# read in data
#pairs_path = "/path_to_long_read_pipeline/results/annos/genes_antisense_info.txt"
pairs <- read_delim(pairs_path)

# plot pie chart showing proportion of genes with antisense transcripts

pairs_summary = pairs %>%
  group_by(withAnti) %>%
  summarise(n = n()) %>%
  mutate(prop = n / sum(n),
         label = ifelse(withAnti == TRUE, "With antisense", "Without antisense"))
S3A = ggplot(pairs_summary, aes(x = "", y = prop, fill = label)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  theme_void() +
  geom_text(aes(label = scales::percent(prop, accuracy = 0.1)), 
            position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = c("With antisense" = "salmon", "Without antisense" = "grey")) +
  labs(title = "Proportion of genes with antisense transcripts") +

  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))


S4E =  ggplot(pairs,aes(x=(length), fill=withAnti)) + 
  geom_density(position="identity", alpha=0.5) + 
  theme_minimal() + 
  labs(title="Gene length distribution", x="Gene length", y="Density") + 
  scale_fill_manual(values=c("grey", "salmon"), labels=c("Without antisense", "With antisense"), name="Gene type") +
  xlim(0, 10000) 
  
  
  
# save plots
ggsave("supp_figure_S3A_gene_antisense_piechart.pdf", S3A, width = 5, height = 5)
ggsave("supp_figure_S4E_gene_length_distribution.pdf", S4E, width = 7, height = 5)