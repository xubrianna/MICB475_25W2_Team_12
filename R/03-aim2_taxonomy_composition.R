library(phyloseq)
library(tidyverse)

#load and preprocess phyloseq object
phyloseq <- readRDS("data/phyloseq_filtered.rds")

#remove zero-abundance taxa
phyloseq <- prune_taxa(taxa_sums(phyloseq) > 0, phyloseq)

#agglomerate at Genus level
phyloseq_genus <- tax_glom(phyloseq, taxrank = "Genus")

#filter low-prevalence taxa at genus level (present in at least 10% of samples)
#set prevalence threshold (as whole number)
prev_threshold <- ceiling(0.1 * nsamples(phyloseq_genus))
phyloseq_filtered <- prune_taxa(
  apply(otu_table(phyloseq_genus), 1, function(x) sum(x > 0)) >= prev_threshold,
  phyloseq_genus
)

#transform counts to relative abundance
phyloseq_genus <- transform_sample_counts(phyloseq_filtered, function(x) x / sum(x))

#keep top 15 genera across ALL SITES/CONDITIONS, group the rest as "Other"
top_genera <- names(sort(taxa_sums(phyloseq_RA), decreasing = TRUE))[1:15]
phyloseq_top <- prune_taxa(top_genera, phyloseq_genus)
                  
#melt to long format for aggregation 
phy_melt <- psmelt(phyloseq_top)
phy_melt$Genus <- gsub("^g_+", "", phy_melt$Genus)

#aggregate by Host_disease and collection_method
phy_agg <- phy_melt %>%
  group_by(Host_disease, collection_method, Genus) %>%
  summarise(RelAbundance = mean(Abundance), .groups = "drop")

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

ggsave(
  filename = "results/aim2/03_taxa_barplot.png",
  plot = gg_agg,
  width = 12,
  height = 8,
  dpi = 300
)



# label "Other" genera
other_genera <- setdiff(taxa_names(phyloseq_genus), top_genera)
tax_tab <- tax_table(phyloseq_genus)
tax_tab[other_genera, "Genus"] <- "Other"

# update phyloseq's taxa table
tax_table(phyloseq_genus) <- tax_tab       


rect <- subset_samples(phyloseq_genus, env_medium == 'rectal')
rect <- prune_taxa(taxa_sums(rect) > 0, rect)
rect <- prune_samples(sample_sums(rect) > 0, rect)

r <- plot_bar(rect, fill = "Genus") +
  facet_grid( . ~ Host_disease, scales = "free_x") +
  theme_classic() +
  theme(
    plot.margin = margin(15, 15, 15, 15),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 13),
    axis.title.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    plot.title = element_text(size = 20),
    strip.text = element_text(size = 14)) +
  labs(
    title = "Taxonomic Composition Across Disease Status: Rectal",
    x = "Sample",
    y = "Relative Abundance"
  )

ggsave(
  filename = "results/aim2/01_taxa_barplot_rectal.png",
  plot = r,
  width = 15,
  height = 10,
  dpi = 300
)

vag <- subset_samples(phyloseq_genus, env_medium == 'vaginal')
vag <- prune_taxa(taxa_sums(vag) > 0, vag)
vag <- prune_samples(sample_sums(vag) > 0, vag)

v <- plot_bar(vag, fill = "Genus") +
  facet_grid( . ~ Host_disease, scales = "free_x") +
  theme_classic() +
  theme(
    plot.margin = margin(15, 15, 15, 15),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 13),
    axis.title.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    plot.title = element_text(size = 20),
    strip.text = element_text(size = 14)) +
  labs(
    title = "Taxonomic Composition Across Disease Status: Vaginal",
    x = "Sample",
    y = "Relative Abundance"
  )

ggsave(
  filename = "results/aim2/01_taxa_barplot_vaginal.png",
  plot = v,
  width = 15,
  height = 10,
  dpi = 300
)



