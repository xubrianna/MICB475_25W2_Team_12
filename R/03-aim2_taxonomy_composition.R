library(phyloseq)
library(tidyverse)

#load data
phyloseq <- readRDS("data/phyloseq_filtered.rds")

#remove taxa with zero abundance
phyloseq <- prune_taxa(taxa_sums(phyloseq) > 0, phyloseq)

#filter low-prevalence taxa (present in <10% of samples)
prev_threshold <- 0.1 * nsamples(phyloseq)

phyloseq_filtered <- prune_taxa(
  apply(otu_table(phyloseq), 1, function(x) sum(x > 0)) > prev_threshold,
  phyloseq
)

#transform to relative abundance
phyloseq_RA <- transform_sample_counts(phyloseq_filtered, function(x) x / sum(x))

#aggregate at genus level
phyloseq_genus <- tax_glom(phyloseq_RA, taxrank = "Genus")

#plot taxonomic composition
taxa_barplot <- plot_bar(phyloseq_genus, fill = "Genus") +
  facet_grid(collection_method ~ Host_disease, scales = "free_x") +
  theme_bw() +
  labs(
    title = "Taxonomic Composition Across Disease Status and Body Site",
    x = "Sample",
    y = "Relative Abundance"
  )

#print
taxa_barplot

#save
ggsave(
  "results/aim2/03_taxa_barplot.png",
  plot = taxa_barplot,
  width = 10,
  height = 5,
  units = "in",
  dpi = 300
)





