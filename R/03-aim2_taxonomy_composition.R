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

#agglomerate at Phyla level
phyloseq_phylum <- tax_glom(phyloseq_RA, taxrank = "Phylum")

#melt to long format
phy_melt <- psmelt(phyloseq_genus)

#clean genus names
phy_melt$Genus <- gsub("^g_+", "", phy_melt$Genus)

#identify top 15 genera
top_genera <- phy_melt %>%
  group_by(Genus) %>%
  summarise(total = sum(Abundance)) %>%
  arrange(desc(total)) %>%
  slice(1:15) %>%
  pull(Genus)

#group remaining genera as "Other"
phy_melt$Genus <- ifelse(phy_melt$Genus %in% top_genera, phy_melt$Genus, "Other")

#aggregate by Host_disease and collection_method
phy_agg <- phy_melt %>%
  group_by(Host_disease, collection_method, Genus) %>%
  summarise(RelAbundance = mean(Abundance), .groups = "drop") %>%
  group_by(Host_disease, collection_method) %>%
  mutate(RelAbundance = RelAbundance / sum(RelAbundance))

#plot taxa barplot
gg_agg <- ggplot(phy_agg, aes(x = Host_disease, y = RelAbundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~collection_method, scales = "free_x") +
  theme_bw() +
  scale_y_continuous(
    limits = c(0, 1.05),
    breaks = seq(0, 1, by = 0.25),
    expand = c(0,0)
  ) +
  labs(
    title = "Mean Taxonomic Composition by Disease and Body Site",
    x = "Disease Status",
    y = "Relative Abundance"
  )

#view plot
gg_agg

#save plot
ggsave(
  filename = "results/aim2/03-tax_composition/03_taxa_barplot.png",
  plot = gg_agg,
  width = 12,
  height = 8,
  dpi = 300
)

