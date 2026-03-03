
#####load packages####
library(phyloseq)
library(vegan)
library(ggplot2)
library(ape)
library(tidyverse)
library(FSA)
library(picante)
set.seed(2026)

phylo <- readRDS('data/phyloseq_filtered.rds')

otu <- read_delim("data/feature-table.txt", delim = "\t", skip = 1)
otu_mat <- as.matrix(otu[, -1])
rownames(otu_mat) <- otu$`#OTU ID`

sample_ids <- colnames(otu_mat)
colors <- rainbow(length(sample_ids))

png("results/aim1/00-rarefaction_curve.png")

rarecurve(
  t(otu_mat),
  step = 20,
  col = colors,
  label = FALSE,                  
  cex = 0.6
)

abline(v = 43000, col = "red", lty = 2, lwd = 2)

dev.off()

#rarefy samples
ps_rarefied <- rarefy_even_depth(phylo, sample.size = 43000, rngseed = 2026,
                                 replace = FALSE, verbose = TRUE)

saveRDS(ps_rarefied, "data/phyloseq_rarefied.rds")

# `set.seed(2026)` was used to initialize repeatable random subsampling.
# Please record this for your records so others can reproduce.
# Try `set.seed(2026); .Random.seed` for the full vector
# ...
# 36 samples removedbecause they contained fewer reads than `sample.size`.
# Up to first five removed samples are: 
#   
#   SRR27830968SRR27830980SRR27830990SRR27831003SRR27831008
# ...
# 119OTUs were removed because they are no longer 
# present in any sample after random subsampling
# 
# ...

####diversity metrics####
#extract OTU matrix
otu_mat <- as(otu_table(ps_rarefied), "matrix")

##samples as rows
if(taxa_are_rows(ps_rarefied)) {
  otu_mat <- t(otu_mat)
}

#calculate Faith's PD
alpha_faith <- pd(otu_mat, phy_tree(ps_rarefied), include.root = TRUE)

#add sample names
alpha_faith$sample <- rownames(alpha_faith)

#merge with metadata (use your existing META)
meta_df <- as(sample_data(ps_rarefied), "data.frame") %>%
  rownames_to_column("sample")

alpha_faith <- alpha_faith %>%
  left_join(meta_df, by = "sample")

#check which samples remain
table(alpha_faith$env_medium, useNA="ifany")

# 6. Factor env_medium for consistent facet order

alpha_faith$env_medium <- factor(alpha_faith$env_medium, levels = c("rectal","vaginal"))

#plot Faith's PD by Host Disease, faceted by env_medium
p <- ggplot(alpha_faith, aes(x=Host_disease, y=PD, fill=Host_disease)) +
  geom_boxplot(lwd=0.8) +
  geom_jitter(width=0.2, size=3) +
  facet_wrap(~env_medium) +
  theme_bw() +
  ylab("Faith's Phylogenetic Diversity") +
  xlab("Host Disease") +
  scale_fill_manual(values=c(
    "Control"="steelblue1",
    "CPP"="mediumpurple1",
    "CPP Endo"="lightpink"
  ))

print(p)
ggsave("results/aim1/01_faith_PD_boxplot.png", 
       plot =gg_uni, 
       width = 10, height = 5, units = "in", dpi = 300)

#Kruskal-Wallis per env_medium
kruskal.test(PD ~ Host_disease, data = subset(alpha_faith, env_medium=="vaginal"))
kruskal.test(PD ~ Host_disease, data = subset(alpha_faith, env_medium=="rectal"))


#Dunn's test pairwise comparisons 
dunnTest(PD ~ Host_disease, data = subset(alpha_faith, env_medium=="vaginal"), method="bh") 
dunnTest(PD ~ Host_disease, data = subset(alpha_faith, env_medium=="rectal"), method="bh")

#Faith's PD as violin
p2 <- ggplot(alpha_faith, aes(x = Host_disease, y = PD, fill = Host_disease)) +
  geom_violin(trim = FALSE, alpha = 0.8) +                 
  geom_jitter(width = 0.15, size = 3, alpha = 0.9, color = "black") + 
  facet_wrap(~env_medium) +
  theme_bw() +
  ylab("Faith's Phylogenetic Diversity") +
  xlab("Host Disease") +
  scale_fill_manual(values = group_cols) +                 
  theme(
    legend.position = "none",                              
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  ) +
  scale_y_continuous(
    limits = {
      y_range <- max(alpha_faith$PD) - min(alpha_faith$PD)
      y_mid <- mean(alpha_faith$PD)
      c(y_mid - 0.5 * y_range, y_mid + 0.5 * y_range)
    },
    breaks = seq(
      floor(min(alpha_faith$PD)), 
      ceiling(max(alpha_faith$PD)), 
      by = 1  
    )
  )

# Display plot
print(p2)

# Save plot
ggsave("results/aim1/01_faith_PD_violin.png",
       plot = p2,
       width = 10, height = 5, units = "in", dpi = 300)

# [TODO]
# Multiple testing correction: Benjamini–Hochberg false discovery rate (FDR) 
# 
# Significance cutoff: 
#   FDR-adjusted p < 0.05 

