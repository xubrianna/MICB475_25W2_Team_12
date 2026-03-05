library(phyloseq)
library(tidyverse)

#load phyloseq object
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

#create taxa bar plot by disease and body site
gg_taxa <- plot_bar(phyloseq_top, fill = "Genus") +
  facet_grid(collection_method ~ Host_disease, scales = "free_x") +
  theme_bw() +
  labs(
    title = "Taxonomic Composition Across Disease Status and Body Site",
    x = "Sample",
    y = "Relative Abundance"
  )

#view plot
gg_taxa

#save plot
ggsave(
  filename = "results/aim1/02_taxa_barplot.png",
  plot = gg_taxa,
  width = 12,
  height = 8,
  dpi = 300
)
