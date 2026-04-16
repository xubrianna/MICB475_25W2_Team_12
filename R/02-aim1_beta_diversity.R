
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
unifrac_dm <- distance(ps_rarefied, method = "unifrac")  # unweighted by default

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

