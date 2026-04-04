library(phyloseq)
library(tidyverse)

# load phyloseq object
phyloseq <- readRDS("data/phyloseq_filtered.rds")

# remove zero-abundance taxa
phyloseq <- prune_taxa(taxa_sums(phyloseq) > 0, phyloseq)

# agglomerate at Phylum level
phyloseq <- tax_glom(phyloseq, taxrank = "Phylum")

# filter low-prevalence taxa (present in at least 10% of samples)
prev_threshold <- ceiling(0.1 * nsamples(phyloseq))

phyloseq_filtered <- prune_taxa(
  apply(otu_table(phyloseq), 1, function(x) sum(x > 0)) >= prev_threshold,
  phyloseq
)

# transform counts to relative abundance
phyloseq_RA <- transform_sample_counts(phyloseq_filtered, function(x) x / sum(x))



# melt to long format
phy_melt <- psmelt(phyloseq_RA)

# clean phylum names
phy_melt$Phylum <- gsub("^p_+", "", phy_melt$Phylum)

top_phyla <- phy_melt %>%
  group_by(Phylum) %>%
  summarise(total = sum(Abundance)) %>%
  arrange(desc(total)) %>%
  pull(Phylum)

# group remaining phyla as "Other"
## phy_melt$Phylum <- ifelse(phy_melt$Phylum %in% top_phyla, phy_melt$Phylum, "Other")

# aggregate by Host_disease and env_medium
phy_agg <- phy_melt %>%
  group_by(Host_disease, env_medium, Phylum) %>%
  summarise(RelAbundance = mean(Abundance), .groups = "drop") %>%
  group_by(Host_disease, env_medium) %>%
  mutate(RelAbundance = RelAbundance / sum(RelAbundance))

# plot taxa barplot
gg_agg <- ggplot(phy_agg, aes(x = Host_disease, y = RelAbundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~env_medium, scales = "free_x") +
  theme_classic() +
  scale_y_continuous(
    limits = c(0, 1.05),
    breaks = seq(0, 1, by = 0.25),
    expand = c(0, 0)
  ) +
  labs(
    x = "Disease Status",
    y = "Relative Abundance"
  ) +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 18),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 18),
    strip.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    plot.margin = margin(15, 15, 15, 15)
  )

# view plot
gg_agg

# save plot
ggsave(
  filename = "results/aim2/03-tax_composition/03_taxa_phylum_barplot.png",
  plot = gg_agg,
  width = 12,
  height = 8,
  dpi = 300
)


####### Genus
# load phyloseq object
phyloseq <- readRDS("data/phyloseq_filtered.rds")

# remove zero-abundance taxa
phyloseq <- prune_taxa(taxa_sums(phyloseq) > 0, phyloseq)

# agglomerate at Genus level
phyloseq <- tax_glom(phyloseq, taxrank = "Genus")

# filter low-prevalence taxa (present in at least 10% of samples)
prev_threshold <- ceiling(0.1 * nsamples(phyloseq))

phyloseq_filtered <- prune_taxa(
  apply(otu_table(phyloseq), 1, function(x) sum(x > 0)) >= prev_threshold,
  phyloseq
)

# transform counts to relative abundance
phyloseq_RA <- transform_sample_counts(phyloseq_filtered, function(x) x / sum(x))


# melt to long format
phy_melt <- psmelt(phyloseq_RA)

# clean phylum names
phy_melt$Genus <- gsub("^g_+", "", phy_melt$Genus)

top_genera <- phy_melt %>%
  group_by(Genus) %>%
  summarise(total = sum(Abundance)) %>%
  arrange(desc(total)) %>%
  slice_head(n = 10) %>%
  pull(Genus)


phy_melt$Genus <- ifelse(phy_melt$Genus %in% top_genera, phy_melt$Genus, "Other")

# aggregate by Host_disease and env_medium
phy_agg <- phy_melt %>%
  group_by(Host_disease, env_medium, Genus) %>%
  summarise(RelAbundance = mean(Abundance), .groups = "drop") %>%
  group_by(Host_disease, env_medium) %>%
  mutate(RelAbundance = RelAbundance / sum(RelAbundance))

# plot taxa barplot
gg_agg <- ggplot(phy_agg, aes(x = Host_disease, y = RelAbundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~env_medium, scales = "free_x") +
  theme_classic() +
  scale_y_continuous(
    limits = c(0, 1.05),
    breaks = seq(0, 1, by = 0.25),
    expand = c(0, 0)
  ) +
  labs(
    x = "Disease Status",
    y = "Relative Abundance"
  ) +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 18),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 18),
    strip.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    plot.margin = margin(15, 15, 15, 15)
  )

# view plot
gg_agg



# save plot
ggsave(
  filename = "results/aim2/03-tax_composition/03_taxa_genera_barplot.png",
  plot = gg_agg,
  width = 12,
  height = 8,
  dpi = 300
)

