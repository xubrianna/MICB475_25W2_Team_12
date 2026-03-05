library(phyloseq)
library(tidyverse)

#load and preprocess phyloseq object
phyloseq <- readRDS("data/phyloseq_filtered.rds")

#remove zero-abundance taxa
phyloseq <- prune_taxa(taxa_sums(phyloseq) > 0, phyloseq)

#filter low-prevalence taxa (present in at least 10% of samples)
prev_threshold <- 0.1 * nsamples(phyloseq)
phyloseq_filtered <- prune_taxa(
  apply(otu_table(phyloseq), 1, function(x) sum(x > 0)) > prev_threshold,
  phyloseq
)

#transform counts to relative abundance
phyloseq_RA <- transform_sample_counts(phyloseq_filtered, function(x) x / sum(x))

#agglomerate at Genus level
phyloseq_genus <- tax_glom(phyloseq_RA, taxrank = "Genus")

#keep top 15 genera, group the rest as "Other"
top_genera <- names(sort(taxa_sums(phyloseq_genus), decreasing = TRUE))[1:15]
phyloseq_top <- prune_taxa(top_genera, phyloseq_genus)

#melt to long format for aggregation 
phy_melt <- psmelt(phyloseq_top)

#aggregate by Host_disease and collection_method
phy_agg <- phy_melt %>%
  group_by(Host_disease, collection_method, Genus) %>%
  summarise(RelAbundance = mean(Abundance), .groups = "drop")

#plot taxa barplot
gg_agg <- ggplot(phy_agg, aes(x = Host_disease, y = RelAbundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~collection_method, scales = "free_x") +
  theme_bw() +
  scale_fill_brewer(palette = "Set3") +
  labs(
    title = "Mean Taxonomic Composition by Disease and Body Site",
    x = "Disease Status",
    y = "Relative Abundance"
  )

#view plot
gg_agg

#save plot
ggsave(
  filename = "results/aim2/02_taxa_barplot.png",
  plot = gg_agg,
  width = 12,
  height = 8,
  dpi = 300
)

