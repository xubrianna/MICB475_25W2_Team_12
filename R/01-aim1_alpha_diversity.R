##### load packages ####
library(phyloseq)
library(vegan)
library(ggplot2)
library(ape)
library(tidyverse)
library(FSA)
library(picante)
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

alpha_faith$env_medium <- factor(alpha_faith$env_medium, levels = c("rectal","vaginal"))

group_cols <- c(
  "Control" = "steelblue1",
  "CPP" = "mediumpurple1",
  "CPP Endo" = "lightpink"
)

##### Faith's PD boxplot ####
p <- ggplot(alpha_faith, aes(x = Host_disease, y = PD, fill = Host_disease)) +
  geom_boxplot(lwd = 0.8) +
  geom_jitter(width = 0.2, size = 3, alpha = 0.9, color = "black") +
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
    limits = c(0, 25),
    breaks = seq(0, 25, by = 2)
  )

print(p)
ggsave("results/aim1/01_faith_PD_boxplot.png", plot = p, width = 10, height = 5, units = "in", dpi = 300)

##### Faith's PD violin plot ####
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
  scale_y_continuous(limits = c(0, 25), breaks = seq(0, 25, by = 2))

print(p2)
ggsave("results/aim1/01_faith_PD_violin.png", plot = p2, width = 10, height = 5, units = "in", dpi = 300)

##### Kruskal-Wallis per site (assigned to objects) ####
kruskal_vaginal <- kruskal.test(PD ~ Host_disease, data = subset(alpha_faith, env_medium=="vaginal"))
kruskal_rectal  <- kruskal.test(PD ~ Host_disease, data = subset(alpha_faith, env_medium=="rectal"))

##### Clean Dunn's test with BH FDR correction ####
sites <- c("vaginal", "rectal")
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

