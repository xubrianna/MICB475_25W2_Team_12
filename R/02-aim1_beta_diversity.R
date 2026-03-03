
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
                                                levels = c("Control","CPP","CPP Endo"))
sample_data(ps_rarefied)$env_medium <- factor(meta_df$env_medium,
                                              levels = c("rectal","vaginal"))  




group_cols <- c(
  "Control"  = "steelblue1",
  "CPP"      = "mediumpurple1",
  "CPP Endo" = "lightpink"
)

gg_uni <- plot_ordination(ps_rarefied, pcoa_uni, type="samples",
                          color="Host_disease") +
  geom_point(size=3.5, alpha=0.9) +
  theme_bw() +
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

meta_feces <- subset(meta_df, env_medium == "rectal")
meta_vagina <- subset(meta_df, env_medium == "vaginal")

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

