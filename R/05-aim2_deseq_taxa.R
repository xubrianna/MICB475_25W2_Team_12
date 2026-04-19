library(DESeq2)
library(phyloseq)
library(EnhancedVolcano)
library(ggrepel)
library(patchwork)
library(dplyr)
library(tibble)

phylo <- readRDS('data/phyloseq_filtered.rds')

# subset into body site 
phylo_rectal <- subset_samples(phylo, env_medium == "rectal")
phylo_vaginal <- subset_samples(phylo, env_medium == "vaginal")

# datasets contain zeroes
phylo_rectal <- transform_sample_counts(phylo_rectal, function(x) x+1)
phylo_vaginal <- transform_sample_counts(phylo_vaginal, function(x) x+1)

# transform into DeSeq object
deseq_rectal <- phyloseq_to_deseq2(phylo_rectal, ~`Host_disease`)
deseq_vaginal <- phyloseq_to_deseq2(phylo_vaginal, ~`Host_disease`)

# run DESeq model
DESEQ_rect <- DESeq(deseq_rectal)
DESEQ_vag <- DESeq(deseq_vaginal)

### CASE-CONTROL COMPARISIONS: sanity check 
## CPP Only vs Control
res_cpp_control_rect <- results(DESEQ_rect, tidy=TRUE, 
                    contrast = c("Host_disease","CPP Only","Control"))

res_cpp_control_vag <- results(DESEQ_vag, tidy=TRUE, 
                   contrast = c("Host_disease","CPP Only","Control"))

df_res_cpp_control_rect <- as.data.frame(res_cpp_control_rect)
df_res_cpp_control_vag <- as.data.frame(res_cpp_control_vag)

taxdf_rectal <- as.data.frame(tax_table(phylo_rectal))

taxdf_vaginal <- as.data.frame(tax_table(phylo_vaginal))

df_cpp_control_rectal_taxa <- cbind(df_res_cpp_control_rect, taxdf_rectal)
df_cpp_control_vaginal_taxa <- cbind(df_res_cpp_control_vag, taxdf_vaginal)

## write_tsv(df_cpp_control_rectal_taxa, "results/aim2/05-deseq2/05-cpp_control_rect_results.tsv")
## write_tsv(df_cpp_control_vaginal_taxa, "results/aim2/05-deseq2/05-cpp_control_vag_results.tsv")


ggplot(res_cpp_control_rect) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj))) +
  theme_classic()
  

ggplot(res_cpp_control_vag) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj))) +
  theme_classic()

cpp_control_rect <- res_cpp_control_rect %>%
  mutate(significant = ifelse(padj<0.05 & abs(log2FoldChange)>2, "Significant", "Non-significant")) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant)) +
  scale_color_manual(values = c("Significant" = "red", "Non-significant" = "grey70")) +
  theme_classic() +
  xlim(-7.5, 7.5) + ylim(0, 20)+
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 18),
    plot.margin = margin(15, 15, 15, 15),
    legend.position = "none",
    strip.text = element_text(size = 18),
    panel.border = element_rect(color = "black", fill = NA),
  ) +
  geom_vline(xintercept = 2, linetype = "dashed")+
  geom_vline(xintercept = -2, linetype = "dashed") +
  geom_abline(slope = 0, intercept = -log10(0.05),linetype = "dashed") +
  labs(
    x = "Log2 Fold Change",
    y = "-Log10 (Adjusted P-value)",
  ) 

cpp_control_vag <- res_cpp_control_vag %>%
  mutate(significant = ifelse(padj<0.05 & abs(log2FoldChange)>2, "Significant", "Non-significant")) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant)) +
  scale_color_manual(values = c("Significant" = "red", "Non-significant" = "grey70")) +
  theme_classic() +
  xlim(-7.5, 7.5) + ylim(0, 20)+
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 18),
    plot.margin = margin(15, 15, 15, 15),
    legend.position = "none",
    strip.text = element_text(size = 18),
    panel.border = element_rect(color = "black", fill = NA),
    ) +
  geom_vline(xintercept = 2, linetype = "dashed")+
  geom_vline(xintercept = -2, linetype = "dashed") +
  geom_abline(slope = 0, intercept = -log10(0.05),linetype = "dashed") +
  labs(
    x = "Log2 Fold Change",
    y = "-Log10 (Adjusted P-value)",
  ) 


cpp_control_vag <- cpp_control_vag + ggtitle("Vaginal") + theme(plot.title = element_text(hjust = 0.5,  size= 18))
cpp_control_rect <- cpp_control_rect + ggtitle("Rectal") + theme(plot.title = element_text(hjust = 0.5, size= 18))


bar_panel <- cpp_control_rect + cpp_control_vag
bar_panel

ggsave(plot= bar_panel,"results/aim2/05-deseq2/05-cpp_control_contrast.png",
       width = 10, height = 7, units = "in", dpi = 300)

## CPP Endo vs Control
res_cpp_endo_control_rect <- results(DESEQ_rect, tidy=TRUE, 
                                contrast = c("Host_disease","CPP Endo","Control"))
res_cpp_endo_control_vag <- results(DESEQ_vag, tidy=TRUE, 
                               contrast = c("Host_disease","CPP Endo","Control"))


df_res_cpp_endo_control_rect <- as.data.frame(res_cpp_endo_control_rect)
df_res_cpp_endo_control_vag <- as.data.frame(res_cpp_endo_control_vag)

df_cpp_endo_control_rectal_taxa <- cbind(df_res_cpp_endo_control_rect, taxdf_rectal)
df_cpp_endo_control_vaginal_taxa <- cbind(df_res_cpp_endo_control_vag, taxdf_vaginal)


# write_tsv(df_cpp_endo_control_rectal_taxa, "results/aim2/05-deseq2/05-cpp_endo_control_rect_results.tsv")
# write_tsv(df_cpp_endo_control_vaginal_taxa, "results/aim2/05-deseq2/05-cpp_endo_control_vag_results.tsv")


ggplot(res_cpp_endo_control_rect) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj)))+
  theme_classic()

ggplot(res_cpp_endo_control_vag) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj))) +
  theme_classic()


cpp_endo_control_rect <- res_cpp_endo_control_rect %>%
  mutate(significant = ifelse(padj<0.05 & abs(log2FoldChange)>2, "Significant", "Non-significant")) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant)) +
  scale_color_manual(values = c("Significant" = "red", "Non-significant" = "grey70")) +
  theme_classic() +
  xlim(-7.5, 7.5) + ylim(0, 20)+
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 18),
    plot.margin = margin(15, 15, 15, 15),
    legend.position = "none",
    strip.text = element_text(size = 18),
    panel.border = element_rect(color = "black", fill = NA),
  ) +
  geom_vline(xintercept = 2, linetype = "dashed")+
  geom_vline(xintercept = -2, linetype = "dashed") +
  geom_abline(slope = 0, intercept = -log10(0.05),linetype = "dashed")+
  labs(
    x = "Log2 Fold Change",
    y = "-Log10 (Adjusted P-value)",
  ) 


cpp_endo_control_vag <- res_cpp_endo_control_vag %>%
  mutate(significant = ifelse(padj<0.05 & abs(log2FoldChange)>2, "Significant", "Non-significant")) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant)) +
  scale_color_manual(values = c("Significant" = "red", "Non-significant" = "grey70")) +
  theme_classic() +
  xlim(-7.5, 7.5) + ylim(0, 20)+
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 18),
    plot.margin = margin(15, 15, 15, 15),
    legend.position = "none",
    strip.text = element_text(size = 18),
    panel.border = element_rect(color = "black", fill = NA),
  ) +
  geom_vline(xintercept = 2, linetype = "dashed")+
  geom_vline(xintercept = -2, linetype = "dashed") +
  geom_abline(slope = 0, intercept = -log10(0.05),linetype = "dashed") +
  labs(
    x = "Log2 Fold Change",
    y = "-Log10 (Adjusted P-value)",
  ) 


cpp_endo_control_vag <- cpp_endo_control_vag + ggtitle("Vaginal") + theme(plot.title = element_text(hjust = 0.5, size= 18))
cpp_endo_control_rect <- cpp_endo_control_rect + ggtitle("Rectal") + theme(plot.title = element_text(hjust = 0.5, size= 18))

bar_panel_endo_control <- cpp_endo_control_rect + cpp_endo_control_vag
bar_panel_endo_control

ggsave(plot= bar_panel_endo_control,"results/aim2/05-deseq2/05-cpp_endo_control_contrast.png",
       width = 10, height = 7, units = "in", dpi = 300)


### CPP - CPP Endo COMPARISONS: 
res_rect <- results(DESEQ_rect, tidy=TRUE, 
                    contrast = c("Host_disease","CPP Endo","CPP Only"))
res_vag <- results(DESEQ_vag, tidy=TRUE, 
                   contrast = c("Host_disease","CPP Endo","CPP Only"))

df_res_rect <- as.data.frame(res_rect)
df_res_vag <- as.data.frame(res_vag)

df_res_rect_taxa <- cbind(df_res_rect, taxdf_rectal)
df_res_vag_taxa <- cbind(df_res_vag, taxdf_vaginal)


# write_tsv(df_res_rect_taxa, "results/aim2/05-deseq2/05-cpp_cpp_endo_rect_results.tsv")
# write_tsv(df_res_vag_taxa, "results/aim2/05-deseq2/05-cpp_cpp_endo_vag_results.tsv")


ggplot(res_rect) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj))) +
  theme_classic()


ggplot(res_vag) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj))) +
  theme_classic()


# Highlight unique significant genera in CPP-endo/CPP Only (Rectal) 
# Get significant genera for each contrast
sig_genera_cpp_control_rect <- df_cpp_control_rectal_taxa %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 2) %>%
  pull(Genus) %>% unique()
sig_genera_cpp_endo_control_rect <- df_cpp_endo_control_rectal_taxa %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 2) %>%
  pull(Genus) %>% unique()
sig_genera_rect <- df_res_rect_taxa %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 2) %>%
  pull(Genus) %>% unique()

# Add highlight column
df_res_rect_taxa <- df_res_rect_taxa %>%
  mutate(
    significant = ifelse(padj < 0.05 & abs(log2FoldChange) > 2, "Significant", "Non-significant"),
    unique_genus = ifelse(
      significant == "Significant" &
      !(Genus %in% c(sig_genera_cpp_control_rect, sig_genera_cpp_endo_control_rect)),
      "Unique", "NotUnique"
    )
  )

cpp_endo_cpp_rect <- ggplot(df_res_rect_taxa) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj),
                 col=case_when(
                   unique_genus == "Unique" ~ "Unique",
                   significant == "Significant" ~ "Significant",
                   TRUE ~ "Non-significant"
                 ))) +
  xlim(-7.5, 7.5) + ylim(0, 20)+
  geom_text_repel(
    data = df_res_rect_taxa %>% 
      filter(unique_genus == "Unique") %>% 
      mutate(Genus = gsub("g__", "", Genus)),
    aes(
      x = log2FoldChange,
      y = -log10(padj),
      label = Genus
    ),
    colour = "blue",
    size = 3,
    max.overlaps = Inf
  ) +
  scale_color_manual(values = c("Unique" = "blue", "Significant" = "red", "Non-significant" = "grey70")) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 18),
    plot.margin = margin(15, 15, 15, 15),
    legend.position = "none",
    strip.text = element_text(size = 18),
    panel.border = element_rect(color = "black", fill = NA),
  ) +
  geom_vline(xintercept = 2, linetype = "dashed")+
  geom_vline(xintercept = -2, linetype = "dashed") +
  geom_abline(slope = 0, intercept = -log10(0.05),linetype = "dashed") +
  labs(
    x = "Log2 Fold Change",
    y = "-Log10 (Adjusted P-value)",
  ) 


# Highlight unique significant genera in CPP-endo/CPP Only (Vaginal) 
sig_genera_cpp_control_vag <- df_cpp_control_vaginal_taxa %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 2) %>%
  pull(Genus) %>% unique()

sig_genera_cpp_endo_control_vag <- df_cpp_endo_control_vaginal_taxa %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 2) %>%
  pull(Genus) %>% unique()
sig_genera_vag <- df_res_vag_taxa %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 2) %>%
  pull(Genus) %>% unique()

df_res_vag_taxa <- df_res_vag_taxa %>%
  mutate(
    significant = ifelse(padj < 0.05 & abs(log2FoldChange) > 2, "Significant", "Non-significant"),
    unique_genus = ifelse(
      significant == "Significant" &
      !(Genus %in% c(sig_genera_cpp_control_vag, sig_genera_cpp_endo_control_vag)),
      "Unique", "NotUnique"
    )
  )

cpp_endo_cpp_vag <- ggplot(df_res_vag_taxa) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj),
                 col=case_when(
                   unique_genus == "Unique" ~ "Unique",
                   significant == "Significant" ~ "Significant",
                   TRUE ~ "Non-significant"
                 ))) +
  geom_text_repel(
    data = df_res_vag_taxa %>% 
      filter(unique_genus == "Unique") %>% 
      mutate(Genus = gsub("g__", "", Genus)),
    aes(
      x = log2FoldChange,
      y = -log10(padj),
      label = Genus
    ),
    colour = "blue",
    size = 3,
    max.overlaps = Inf
  ) +
  xlim(-7.5, 7.5) + ylim(0, 20)+
  scale_color_manual(values = c("Unique" = "blue", "Significant" = "red", "Non-significant" = "grey70")) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 18),
    plot.margin = margin(15, 15, 15, 15),
    legend.position = "none",
    strip.text = element_text(size = 18),
    panel.border = element_rect(color = "black", fill = NA),
  ) +
  geom_vline(xintercept = 2, linetype = "dashed")+
  geom_vline(xintercept = -2, linetype = "dashed") +
  geom_abline(slope = 0, intercept = -log10(0.05),linetype = "dashed") +
  labs(
    x = "Log2 Fold Change",
    y = "-Log10 (Adjusted P-value)",
  ) 

cpp_endo_cpp_vag <- cpp_endo_cpp_vag + ggtitle("Vaginal") + theme(plot.title = element_text(hjust = 0.5, size= 18))
cpp_endo_cpp_rect <- cpp_endo_cpp_rect + ggtitle("Rectal") + theme(plot.title = element_text(hjust = 0.5,  size= 18))


bar_panel_endo_cpp <- cpp_endo_cpp_rect  + cpp_endo_cpp_vag
bar_panel_endo_cpp

ggsave(plot= bar_panel_endo_cpp,"results/aim2/05-deseq2/05-cpp_cpp_endo_contrast.png",
       width = 10, height = 7, units = "in", dpi = 300)


## panel for final plot
p <- plot_grid(bar_panel, bar_panel_endo_control, bar_panel_endo_cpp, ncol = 1, labels = 'AUTO', label_size = 20)

ggsave(
  "results/aim2/05-deseq2/05-final_figure.png",
  plot = p,
  width = 13, height = 13, units = "in", dpi = 300)



# To get table of results
sigASVs_rect <- res_rect %>% 
  filter(padj<0.05 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)

sigASVs_vag <- res_vag %>% 
  filter(padj<0.05 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)


# Get only asv names
sigASVs_vec_rect <- sigASVs_rect %>%
  pull(ASV)

sigASVs_vec_vag <- sigASVs_vag %>%
  pull(ASV)


write_tsv(sigASVs_rect, "results/aim2/05-deseq2/05-sig_rect_results.tsv")
write_tsv(sigASVs_vag, "results/aim2/05-deseq2/05-sig_vag_results.tsv")


### PRUNE PHYLOSEQ FILE
# RECTAL
DESeq_pruned_rect <- prune_taxa(sigASVs_vec_rect, phylo_rectal)

sigASVs_rect <- tax_table(DESeq_pruned_rect) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_rect) %>%
  arrange(log2FoldChange) %>%
  mutate(Phylum = gsub("p__", "", Phylum)) %>%
  mutate(Phylum = make.unique(Phylum)) %>%
  mutate(Phylum = factor(Phylum, levels=unique(Phylum)))


# Compute shared fill range for unified legend
shared_lfc_range <- range(c(sigASVs_rect$log2FoldChange, sigASVs_vag$log2FoldChange), na.rm = TRUE)

bar_rect <- ggplot(sigASVs_rect) +
  geom_bar(aes(x=log2FoldChange, y=Phylum, fill=log2FoldChange), stat="identity")+
  geom_errorbar(aes(y=Phylum, xmin=log2FoldChange-lfcSE, xmax=log2FoldChange+lfcSE)) +
  scale_fill_gradient(low = "#1B98E026", high = "#225ea8", limits = shared_lfc_range) +
  theme_classic() +
  theme(axis.text.y = element_text(hjust=1, vjust=0.5),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 16),
        axis.title.x = element_text(size = 18, margin = margin(t = 20)),
        axis.title.y = element_text(size = 18, margin = margin(r = 20)),
        legend.text = element_text(size = 14),
        strip.text = element_text(size = 18),
        strip.background = element_rect(fill = "grey80", color = "black"),
        panel.border = element_rect(color = "black", fill = NA),
        plot.margin = margin(15, 15, 15, 15),
        plot.title = element_text(size = 18, hjust = 0.5),
        ) +
  labs(
    x = "Log2 Fold Change",
    title = "Rectal",
    fill= "Log2 Fold Change"
  
  ) 
 
bar_rect

ggsave(plot = bar_rect, "results/aim2/05-deseq2/05-sigASVs_rect.png",
       width = 10, height = 20, units = "in", dpi = 300)


# VAGINAL
DESeq_pruned_vag <- prune_taxa(sigASVs_vec_vag, phylo_vaginal)

sigASVs_vag <- tax_table(DESeq_pruned_vag) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_vag) %>%
  arrange(log2FoldChange) %>%
  mutate(Phylum = gsub("p__", "", Phylum)) %>%
  mutate(Phylum = make.unique(Phylum)) %>%
  mutate(Phylum = factor(Phylum, levels=unique(Phylum)))


bar_vag <- ggplot(sigASVs_vag) +
  geom_bar(aes(x=log2FoldChange, y=Phylum, fill=log2FoldChange), stat="identity")+
  geom_errorbar(aes(y=Phylum, xmin=log2FoldChange-lfcSE, xmax=log2FoldChange+lfcSE)) +
  scale_fill_gradient(low = "#1B98E026", high = "#225ea8", limits = shared_lfc_range) +
  theme_classic()+
  theme(axis.text.y = element_text(hjust=1, vjust=0.5),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 16),
        axis.title.x = element_text(size = 18, margin = margin(t = 20)),
        axis.title.y = element_text(size = 18, margin = margin(r = 20)),
        legend.text = element_text(size = 14),
        strip.text = element_text(size = 18),
        strip.background = element_rect(fill = "grey80", color = "black"),
        panel.border = element_rect(color = "black", fill = NA),
        plot.margin = margin(15, 15, 15, 15),
        plot.title = element_text(size = 18, hjust = 0.5)) +
  labs(
    x = "Log2 Fold Change",
    title = "Vaginal",
    fill= "Log2 Fold Change"
  ) 
bar_vag

ggsave(plot = bar_vag, "results/aim2/05-deseq2/05-sigASVs_vag.png",
       width = 10, height = 20, units = "in", dpi = 300)

# Combined bar plot panel with shared legend
bar_combined <- bar_rect + bar_vag + plot_layout(guides = "collect") &
  theme(legend.position = "right")

ggsave(plot = bar_combined, "results/aim2/05-deseq2/05-sigASVs_combined.png",
       width = 20, height = 20, units = "in", dpi = 300)



