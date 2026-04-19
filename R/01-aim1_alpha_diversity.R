##### load packages ####
library(phyloseq)
library(vegan)
library(ggplot2)
library(ape)
library(tidyverse)
library(FSA)
library(picante)
library(ggpubr)
library(patchwork)
library(cowplot)
library(rstatix)
set.seed(2026)

##### read data ####
phylo <- readRDS('data/phyloseq_filtered.rds')

otu <- read_delim("data/feature-table.txt", delim = "\t", skip = 1)
otu_mat <- as.matrix(otu[, -1])
rownames(otu_mat) <- otu$`#OTU ID`

sample_ids <- colnames(otu_mat)
colors <- rainbow(length(sample_ids))

##### rarefaction curve ####
png("results/aim1/00-rarefaction_curve.png")
rarecurve(t(otu_mat), step = 20, col = colors, label = FALSE, cex = 0.6)
abline(v = 43000, col = "red", lty = 2, lwd = 2)
dev.off()

##### rarefy samples ####
ps_rarefied <- rarefy_even_depth(
  phylo, sample.size = 43000, rngseed = 2026,
  replace = FALSE, verbose = TRUE
)
saveRDS(ps_rarefied, "data/phyloseq_rarefied.rds")

##### diversity metrics ####
otu_mat <- as(otu_table(ps_rarefied), "matrix")
if(taxa_are_rows(ps_rarefied)) otu_mat <- t(otu_mat)

alpha_faith <- pd(otu_mat, phy_tree(ps_rarefied), include.root = TRUE)
alpha_faith$sample <- rownames(alpha_faith)

meta_df <- as(sample_data(ps_rarefied), "data.frame") %>%
  rownames_to_column("sample")
alpha_faith <- alpha_faith %>%
  left_join(meta_df, by = "sample")

alpha_faith$env_medium <- gsub('rectal', 'Rectal',alpha_faith$env_medium)
alpha_faith$env_medium <- gsub('vaginal', 'Vaginal',alpha_faith$env_medium)
alpha_faith$env_medium <- factor(alpha_faith$env_medium, levels = c("Rectal","Vaginal"))



##### Faith's PD boxplot ####
alpha_faith$Host_disease = factor(alpha_faith$Host_disease, levels= c('Control', 'CPP Only', 'CPP Endo'))
group_cols <- c(
  "Control" = "#c7e9b4",   # green
  "CPP Only" = "#41b6c4",       # turquoise
  "CPP Endo" = "#225ea8"   # dark blue
)


                                  
p <- ggplot(alpha_faith, aes(x = Host_disease, y = PD, fill = Host_disease)) +
  geom_boxplot(width = 0.6, linewidth = 0.8,outlier.shape = NA) +
  geom_jitter(width = 0.15,
              size = 2.5,
              color = 'black') +
  facet_wrap(~env_medium) +
  theme_classic() +
  ylab("Faith's Phylogenetic Diversity") +
  xlab("Host Disease") +
  scale_fill_manual(values = group_cols) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 14),
    axis.title.x = element_text(size = 18, margin = margin(t = 20)),
    axis.title.y = element_text(size = 18, margin = margin(r = 20)),
    strip.text = element_text(size = 18),
    strip.background = element_rect(fill = "grey80", color = "black"),
    plot.margin = margin(15, 15, 15, 15),
    panel.border = element_rect(color = "black", fill = NA),
  ) +
  scale_y_continuous(limits = c(0, 25), breaks = seq(0, 25, by = 2)) 
  

print(p)
ggsave("results/aim1/01_faith_PD_boxplot.png", plot = p, width = 10, height = 5, units = "in", dpi = 300)

##### Faith's PD violin plot ####


##### Kruskal-Wallis per site (assigned to objects) ####
kruskal_vaginal <- kruskal.test(PD ~ Host_disease, data = subset(alpha_faith, env_medium=="Vaginal"))
kruskal_rectal  <- kruskal.test(PD ~ Host_disease, data = subset(alpha_faith, env_medium=="Rectal"))


dunn_results <- alpha_faith %>%
  group_by(env_medium) %>%
  dunn_test(PD ~ Host_disease, p.adjust.method = "BH") %>%
  mutate(y.position = c(22, 24, 26, 21, 23, 25))  # 3 comparisons per site


dunn_results %>% select(env_medium, group1, group2, y.position)

p2 <- ggplot(alpha_faith, aes(x = Host_disease, y = PD, fill = Host_disease)) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  # geom_jitter(width = 0.07,
  #             size = 2.5,
  #             color = 'black') +
  geom_boxplot(width = 0.07)+
  facet_wrap(~env_medium) +
  theme_classic() +
  ylab("Faith's Phylogenetic Diversity") +
  scale_fill_manual(values = group_cols) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 18, margin = margin(r = 20)),
    strip.text = element_text(size = 18),
    plot.margin = margin(15, 15, 15, 15),
    strip.background = element_rect( color = "black"),
    panel.border = element_rect(color = "black", fill = NA),
  ) +
  scale_y_continuous(limits = c(0, 28), breaks = seq(0, 28, by = 2)) +
  stat_pvalue_manual(
    dunn_results,
    label = "p.adj.signif", 
    tip.length = 0.01,
    facet.by = "env_medium"   # add this
  )

print(p2)
ggsave("results/aim1/01_faith_PD_violin.png", plot = p2, width = 10, height = 5, units = "in", dpi = 300)

##### Clean Dunn's test with BH FDR correction ####
sites <- c("Vaginal", "Rectal")
dunn_results <- list()

for(site in sites) {
  df <- subset(alpha_faith, env_medium == site)
  dunn_res <- dunnTest(PD ~ Host_disease, data = df, method = "none")
  
  pairwise <- dunn_res$res %>%
    select(Comparison, P.unadj) %>%
    mutate(
      P.adj = p.adjust(P.unadj, method = "BH"),
      Significant = P.adj < 0.05
    )
  
  dunn_results[[site]] <- pairwise  # store in list
  cat("\nDunn's test for", site, "samples (BH-adjusted):\n")
  print(pairwise)
}

kruskal_vaginal   # Kruskal-Wallis test for vaginal samples
kruskal_rectal    # Kruskal-Wallis test for rectal samples
dunn_results[["vaginal"]]  # Dunn's test table for vaginal samples
dunn_results[["rectal"]]   # Dunn's test table for rectal samples

##### NEW: FACET BY DISEASE, COMPARE RECTAL vs VAGINAL #####

# Stats per disease
site_stats <- alpha_faith %>%
  group_by(Host_disease) %>%
  summarise(
    p = wilcox.test(PD ~ env_medium)$p.value
  ) %>%
  mutate(
    label = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE ~ "ns"
    ),
    y_pos = 25
  )

# New plot
p3 <- ggplot(alpha_faith, aes(x = env_medium, y = PD, fill = env_medium)) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  geom_boxplot(width = 0.07)+
  # geom_point(position = position_jitter(width = 0.15),
  #            size = 2.5, alpha = 0.9, color = "black") +
  facet_wrap(~Host_disease) +
  theme_classic() +
  ylab("Faith's Phylogenetic Diversity") +
  scale_fill_manual(values = c("Rectal" = "#f1a340", "Vaginal" = "#998ec3")) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 18),
    strip.text = element_text(size = 18),
    axis.title.x = element_blank(),
    strip.background = element_rect( color = "black"),
    panel.border = element_rect(color = "black", fill = NA),
  ) +
  scale_y_continuous(limits = c(0, 27), breaks = seq(0, 25, by = 2)) +
  stat_compare_means(
    comparisons = list(c("Rectal", "Vaginal")),
    method = "wilcox.test",
    label = "p.signif",
    label.y = 25
  )

print(p3)
ggsave("results/aim1/02_faith_PD_by_site_within_disease.png",
       plot = p3, width = 10, height = 5, units = "in", dpi = 300)




alpha_plots <- plot_grid(p2, p3, ncol=1, labels = c('A', 'B'), label_size=20)
alpha_plots





####################################################################################
### BETA DIVERSITY #################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################

#####load packages####
library(phyloseq)
library(vegan)
library(ggplot2)
library(ape)
library(tidyverse)
library(FSA)
library(picante)
set.seed(2026)




ps_rarefied <- readRDS('data/phyloseq_rarefied.rds')
#beta diversity
unifrac_dm <- phyloseq::distance(ps_rarefied, method = "unifrac")  # unweighted by default

#PCoA ordination
pcoa_uni <- ordinate(ps_rarefied,
                     method = "PCoA",
                     distance = "unifrac")

#prepare metadata
meta_df <- as(sample_data(ps_rarefied), "data.frame")

sample_data(ps_rarefied)$Host_disease <- factor(meta_df$Host_disease,
                                                levels = c("Control","CPP Only","CPP Endo"))
meta_df$env_medium <- gsub("rectal", "Rectal", meta_df$env_medium)
meta_df$env_medium <- gsub("vaginal", "Vaginal", meta_df$env_medium)
sample_data(ps_rarefied)$env_medium <- factor(meta_df$env_medium,
                                              levels = c("Rectal","Vaginal"))  




group_cols <- c(
  "Control" = "#c7e9b4",   # green
  "CPP Only" = "#41b6c4",       # turquoise
  "CPP Endo" = "#225ea8"   # dark blue
)

####first attempt####
gg_uni <- plot_ordination(ps_rarefied, pcoa_uni, type="samples",
                          color="Host_disease") +
  geom_point(size=3.5, alpha=0.9) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 18),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 18),
    strip.background = element_rect(fill = "grey80", color = "black"),
    plot.margin = margin(15, 15, 15, 15),
    panel.border = element_rect(color = "black", fill = NA),
  ) +
  facet_wrap(~ env_medium) +
  scale_color_manual(values = group_cols) +
  labs(color="Host Disease", shape="Site")

gg_uni

ggsave("results/aim1/02-pcoa.png", 
       plot =gg_uni, 
       width = 10, height = 5, units = "in", dpi = 300)

# PERMANOVA overall
permanova_all <- adonis2(unifrac_dm ~ env_medium + Host_disease,
                         data = meta_df,
                         permutations = 999,
                         by = "margin")   # tests each term controlling for the other
permanova_all

# PERMANOVA stratified by site

meta_feces <- subset(meta_df, env_medium == "Rectal")
meta_vagina <- subset(meta_df, env_medium == "Vaginal")

dm_feces  <- as.dist(as.matrix(unifrac_dm)[rownames(meta_feces), rownames(meta_feces)])
dm_vagina <- as.dist(as.matrix(unifrac_dm)[rownames(meta_vagina), rownames(meta_vagina)])

set.seed(2026)
permanova_feces <- adonis2(dm_feces ~ Host_disease,
                           data = meta_feces,
                           permutations = 999)
permanova_feces

set.seed(2026)
permanova_vagina <- adonis2(dm_vagina ~ Host_disease,
                            data = meta_vagina,
                            permutations = 999)
permanova_vagina

####revised####
# Extract % variation for axes
eig_vals <- pcoa_uni$values$Eigenvalues
var_explained <- round(eig_vals / sum(eig_vals) * 100, 1)

# PERMANOVA p-value (overall)
pval <- permanova_all$`Pr(>F)`[2]  # Host_disease term

# Create new PCoA plot object
gg_uni_ellipse <- (
  plot_ordination(ps_rarefied, pcoa_uni, type = "samples", color = "Host_disease")
  + geom_point(aes(shape = env_medium), size = 3.5, alpha = 0.9)
  + stat_ellipse(aes(group = interaction(Host_disease, env_medium), color = Host_disease),
                 level = 0.95, linewidth = 1)
  + theme_classic()
  + theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 18),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.margin = margin(15, 15, 15, 15)
  )
  + scale_color_manual(values = group_cols)
  + scale_shape_manual(values = c(16, 17)) 
  + labs(
    color = "Host Disease",
    shape = "Body Site",
    x = sprintf("PCoA Axis 1 (%.1f%% of community dissimilarity)", var_explained[1]),
    y = sprintf("PCoA Axis 2 (%.1f%% of community dissimilarity)", var_explained[2])
  )
  + scale_x_continuous(
    expand = expansion(mult = 0.25),
    breaks = pretty(pcoa_uni$vectors[,1], n = 5)
  )
  + scale_y_continuous(
    limits = c(floor(min(pcoa_uni$vectors[,2])*15)/15 - 0.5,
               ceiling(max(pcoa_uni$vectors[,2])*15)/15 + 0.5), 
    breaks = seq(floor(min(pcoa_uni$vectors[,2])*15)/15 - 0.5,
                 ceiling(max(pcoa_uni$vectors[,2])*15)/15 + 0.5,
                 by = 0.25)  # start seq at padded bottom
  )
)

# Display the new plot
gg_uni_ellipse

# Save new plot
ggsave("results/aim1/02-pcoa_ellipse.png",
       plot = gg_uni_ellipse,
       width = 10, height = 7, units = "in", dpi = 300)




####################################################################################
### TAXONOMY COMP #################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################

library(RColorBrewer)
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

# capitalize env_medium labels
phy_melt$env_medium <- gsub("rectal", "Rectal", phy_melt$env_medium)
phy_melt$env_medium <- gsub("vaginal", "Vaginal", phy_melt$env_medium)

# clean phylum names
phy_melt$Phylum <- gsub("^p_+", "", phy_melt$Phylum)

top_phyla <- phy_melt %>%
  group_by(Phylum) %>%
  summarise(total = sum(Abundance)) %>%
  arrange(desc(total)) %>%
  slice_head(n = 10) %>%
  pull(Phylum)

# group remaining phyla as "Other"
phy_melt$Phylum <- ifelse(phy_melt$Phylum %in% top_phyla, phy_melt$Phylum, "Other")

# aggregate by Host_disease and env_medium
phy_agg <- phy_melt %>%
  group_by(Host_disease, env_medium, Phylum) %>%
  summarise(RelAbundance = mean(Abundance), .groups = "drop") %>%
  group_by(Host_disease, env_medium) %>%
  mutate(RelAbundance = RelAbundance / sum(RelAbundance))


phy_agg$Host_disease <- factor(phy_agg$Host_disease, levels = c('Control', 'CPP Only', 'CPP Endo'))



phy_agg$Phylum <- factor(phy_agg$Phylum, levels = c("Actinomycetota",
                                                    "Bacillota",               
                                                    "Bacteroidota",
                                                    "Campylobacterota",       
                                                    "Fusobacteriota",
                                                    "Pseudomonadota",
                                                    "Synergistota",
                                                    "Thermodesulfobacteriota",
                                                    "Verrucomicrobiota",      
                                                    "Other"))

# plot taxa barplot
gg_agg_phy <- ggplot(phy_agg, aes(x = Host_disease, y = RelAbundance, fill = Phylum)) +
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
    strip.background = element_rect( color = "black"),
    panel.border = element_rect(color = "black", fill = NA),
    axis.title.x = element_text(size = 18, margin = margin(t = 20)),
    plot.margin = margin(15, 15, 15, 20)
  ) +
  scale_fill_brewer('Blues', direction = -1)

# view plot
print(gg_agg_phy)

# save plot
ggsave(
  filename = "results/aim2/03-tax_composition/03_taxa_phylum_barplot.png",
  plot = gg_agg_phy,
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

# capitalize env_medium labels
phy_melt$env_medium <- gsub("rectal", "Rectal", phy_melt$env_medium)
phy_melt$env_medium <- gsub("vaginal", "Vaginal", phy_melt$env_medium)

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


phy_agg$Genus <- factor(phy_agg$Genus, levels = c("Bacteroides",
                                                  "Dialister",
                                                  "Escherichia-Shigella",
                                                  "Faecalibacterium",    
                                                  "Fannyhessea",
                                                  "Gardnerella",
                                                  "Hoylesella",
                                                  "Lactobacillus",
                                                  "Peptoniphilus",
                                                  "Prevotella",  
                                                  "Other"))


phy_agg$Host_disease <- factor(phy_agg$Host_disease, levels = c('Control', 'CPP Only', 'CPP Endo'))

blues_11 <- colorRampPalette(rev(brewer.pal(9, "Blues")))(11)


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
    axis.title.x = element_blank(),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 18),
    # strip.background = element_rect( color = "black"),
    # panel.border = element_rect(color = "black", fill = NA),
    plot.margin = margin(15, 15, 15, 20)
  ) + 
  scale_fill_manual(values = blues_11)

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


beta_plots <- plot_grid(gg_uni_ellipse, gg_agg, ncol=1, labels = c('C', 'D'), label_size=20)
beta_plots

final <- plot_grid(alpha_plots, beta_plots, ncol=2, label_size = 20)
final

# save plot
ggsave(
  filename = "results/aim1/03-final.png",
  plot = final,
  width = 20,
  height = 15,
  dpi = 300
)



