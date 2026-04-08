library(tidyverse)
library(phyloseq)
library(ggpicrust2)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(patchwork)
set.seed(2026)


### ADJUST THRESHOLDS
logFC_threshold = 4
p_threshold = 0.01

ps <- readRDS("data/phyloseq_filtered.rds")

meta <- sample_data(ps) %>%
  data.frame() %>%
  rownames_to_column("sample_name")

ko <- read.delim("data/picrust_out/KO_metagenome_out/pred_metagenome_unstrat.tsv", row.names = 1)

rownames(ko) <- gsub("ko:", "", rownames(ko))

rectal_meta <- meta %>%
  dplyr::filter(collection_method == "rectal_swab")
vaginal_meta <- meta %>%
  dplyr::filter(collection_method == "vaginal_swab")
rectal_ko <- ko %>%
  dplyr::select(all_of(rectal_meta$sample_name))
vaginal_ko <- ko %>%
  dplyr::select(all_of(vaginal_meta$sample_name))



rectal_meta$Host_disease  <- factor(rectal_meta$Host_disease)
vaginal_meta$Host_disease <- factor(vaginal_meta$Host_disease)

rectal_meta$Host_disease <- relevel(factor(rectal_meta$Host_disease), ref = "Control")
vaginal_meta$Host_disease <- relevel(factor(vaginal_meta$Host_disease), ref = "Control")

rectal_mat  <- as.matrix(rectal_ko)
storage.mode(rectal_mat) <- "numeric"

vaginal_mat <- as.matrix(vaginal_ko)
storage.mode(vaginal_mat) <- "numeric"

#rectal 

dds_rectal <- DESeqDataSetFromMatrix(
  countData = round(rectal_mat),
  colData   = rectal_meta,
  design    = ~ Host_disease
)

dds_rectal <- DESeq(dds_rectal)
res_rectal <- results(dds_rectal)

res_rectal_df <- as.data.frame(res_rectal)
res_rectal_df$KO <- rownames(res_rectal_df)
res_rectal_df <- res_rectal_df[order(res_rectal_df$padj), ]

sig_rectal <- subset(res_rectal_df, !is.na(padj) & padj < p_threshold)
nrow(sig_rectal)
head(sig_rectal, 10)

#vaginal

dds_vaginal <- DESeqDataSetFromMatrix(
  countData = round(vaginal_mat),
  colData   = vaginal_meta,
  design    = ~ Host_disease
)

dds_vaginal <- DESeq(dds_vaginal)
res_vaginal <- results(dds_vaginal)

res_vaginal_df <- as.data.frame(res_vaginal)
res_vaginal_df$KO <- rownames(res_vaginal_df)
res_vaginal_df <- res_vaginal_df[order(res_vaginal_df$padj), ]


sig_vaginal <- subset(res_vaginal_df, !is.na(padj) & padj < p_threshold)
nrow(sig_vaginal)
head(sig_vaginal, 10)


###################################PCA plots################################################
group_cols <- c(
  "Control" = "#c7e9b4",
  "CPP" = "#41b6c4",
  "CPP Endo" = "#225ea8"
)

rectal_mat_filtered <- rectal_mat[
  apply(rectal_mat, 1, var) != 0,
]

rectal_pca <- prcomp(t(log1p(rectal_mat_filtered)), scale. = TRUE)

rectal_pca_df <- data.frame(rectal_pca$x)
rectal_pca_df$Host_disease <- rectal_meta$Host_disease

rectal_var_explained <- summary(rectal_pca)$importance[2,]

rectal_KO_PCA <- ggplot(
  rectal_pca_df,
  aes(x = PC1, y = PC2, color = Host_disease, fill = Host_disease)
) +
  stat_ellipse(
    geom = "polygon",
    alpha = 0.15,
    color = NA,
    level = 0.95
  ) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95, linewidth = 1) +
  scale_color_manual(values = group_cols) +
  scale_fill_manual(values = group_cols) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 18),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 18),
    strip.background = element_rect(fill = "grey80", color = "black"),
    panel.border = element_rect(color = "black", fill = NA),
    plot.margin = margin(15, 15, 15, 15)
  ) +
  labs(
    x = paste0("PC1 (", round(rectal_var_explained[1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(rectal_var_explained[2] * 100, 1), "%)"),
    color = "Host disease",
    fill = "Host disease"
  )

vaginal_mat_filtered <- vaginal_mat[
  apply(vaginal_mat, 1, var) != 0,
]

vaginal_pca <- prcomp(t(log1p(vaginal_mat_filtered)), scale. = TRUE)

vaginal_pca_df <- data.frame(vaginal_pca$x)
vaginal_pca_df$Host_disease <- vaginal_meta$Host_disease

vaginal_var_explained <- summary(vaginal_pca)$importance[2,]

vaginal_KO_PCA <- ggplot(
  vaginal_pca_df,
  aes(x = PC1, y = PC2, color = Host_disease, fill = Host_disease)
) +
  stat_ellipse(
    geom = "polygon",
    alpha = 0.15,
    color = NA,
    level = 0.95
  ) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95, linewidth = 1) +
  scale_color_manual(values = group_cols) +
  scale_fill_manual(values = group_cols) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 18),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 18),
    strip.background = element_rect(fill = "grey80", color = "black"),
    panel.border = element_rect(color = "black", fill = NA),
    plot.margin = margin(15, 15, 15, 15)
  ) +
  labs(
    x = paste0("PC1 (", round(vaginal_var_explained[1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(vaginal_var_explained[2] * 100, 1), "%)"),
    color = "Host disease",
    fill = "Host disease"
  )
rectal_KO_PCA
vaginal_KO_PCA

combined_KO_PCA <- rectal_KO_PCA + vaginal_KO_PCA
combined_KO_PCA

ggsave(
  "results/aim3/KO_rectal_vaginal_PCA_combined.png",
  plot = combined_KO_PCA,
  width = 16, height = 7, units = "in", dpi = 300
)

ggsave(
  "results/aim3/KO_rectal_vaginal_PCA_combined.png",
  plot = combined_KO_PCA,
  width = 16, height = 7, units = "in", dpi = 300
)

ggsave("results/aim3/KO_rectal_pcoa_ellipse.png",
       plot = rectal_KO_PCA,
       width = 10, height = 7, units = "in", dpi = 300)
ggsave("results/aim3/KO_vaginal_pcoa_ellipse.png",
       plot = vaginal_KO_PCA,
       width = 10, height = 7, units = "in", dpi = 300)


######################################volcano plot#######

#rectal
res_rectal_CPP_vs_Control <- results(
  dds_rectal,
  contrast = c("Host_disease", "CPP", "Control")
)

res_rectal_Endo_vs_Control <- results(
  dds_rectal,
  contrast = c("Host_disease", "CPP Endo", "Control")
)

res_rectal_Endo_vs_CPP <- results(
  dds_rectal,
  contrast = c("Host_disease", "CPP Endo", "CPP")
)

rectal_CPP_vs_Control_df <- as.data.frame(res_rectal_CPP_vs_Control)
rectal_CPP_vs_Control_df$KO <- rownames(rectal_CPP_vs_Control_df)

rectal_Endo_vs_Control_df <- as.data.frame(res_rectal_Endo_vs_Control)
rectal_Endo_vs_Control_df$KO <- rownames(rectal_Endo_vs_Control_df)

rectal_Endo_vs_CPP_df <- as.data.frame(res_rectal_Endo_vs_CPP)
rectal_Endo_vs_CPP_df$KO <- rownames(rectal_Endo_vs_CPP_df)

rectal_CPP_vs_Control_df$significance <- "Not Significant"
rectal_CPP_vs_Control_df$significance[
  !is.na(rectal_CPP_vs_Control_df$padj) &
    rectal_CPP_vs_Control_df$padj < p_threshold &
    abs(rectal_CPP_vs_Control_df$log2FoldChange) > logFC_threshold
] <- "Significant"

rectal_Endo_vs_Control_df$significance <- "Not Significant"
rectal_Endo_vs_Control_df$significance[
  !is.na(rectal_Endo_vs_Control_df$padj) &
    rectal_Endo_vs_Control_df$padj < p_threshold &
    abs(rectal_Endo_vs_Control_df$log2FoldChange) > logFC_threshold
] <- "Significant"

rectal_Endo_vs_CPP_df$significance <- "Not Significant"
rectal_Endo_vs_CPP_df$significance[
  !is.na(rectal_Endo_vs_CPP_df$padj) &
    rectal_Endo_vs_CPP_df$padj < p_threshold &
    abs(rectal_Endo_vs_CPP_df$log2FoldChange) > logFC_threshold
] <- "Significant"

# write df results as tsv
write_tsv(rectal_Endo_vs_CPP_df |> filter(significance == 'Significant'),
          'results/aim3/06-cpp_cpp_endo_significant_KOs_rectal.tsv')
write_tsv(rectal_Endo_vs_Control_df |> filter(significance == 'Significant'),
          'results/aim3/06-cpp_endo_control_significant_KOs_rectal.tsv')
write_tsv(rectal_CPP_vs_Control_df |> filter(significance == 'Significant'),
          'results/aim3/06-cpp_control_significant_KOs_rectal.tsv')



top_labels <- rectal_CPP_vs_Control_df %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::arrange(padj) %>%
  head(10)

volcano_rectal_CPP_vs_Control <- ggplot(
  rectal_CPP_vs_Control_df,
  aes(x = log2FoldChange, y = -log10(padj), color = significance)
) +
  geom_point(alpha = 0.7) +
  geom_text_repel(
    data = top_labels,
    aes(label = KO),
    size = 3
  ) +
  geom_hline(yintercept = -log10(p_threshold), linetype = "dashed") +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed") +
  scale_color_manual(values = c("grey70", "red")) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 18),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 18),
    strip.background = element_rect(fill = "grey80", color = "black"),
    panel.border = element_rect(color = "black", fill = NA),
    plot.margin = margin(15, 15, 15, 15),
    legend.position = "none"
  ) +
  labs(
    x = "Log2 Fold Change",
    y = "-Log10 (Adjusted P-value)"
  )

volcano_rectal_CPP_vs_Control

top_labels_Endo_vs_Control <- rectal_Endo_vs_Control_df %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::arrange(padj) %>%
  head(10)

volcano_rectal_Endo_vs_Control <- ggplot(
  rectal_Endo_vs_Control_df,
  aes(x = log2FoldChange, y = -log10(padj), color = significance)
) +
  geom_point(alpha = 0.7) +
  geom_text_repel(
    data = top_labels_Endo_vs_Control,
    aes(label = KO),
    size = 3
  ) +
  geom_hline(yintercept = -log10(p_threshold), linetype = "dashed") +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed") +
  scale_color_manual(values = c("grey70", "red")) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 18),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 18),
    strip.background = element_rect(fill = "grey80", color = "black"),
    panel.border = element_rect(color = "black", fill = NA),
    plot.margin = margin(15, 15, 15, 15),
    legend.position = "none"
  ) +
  labs(
    x = "Log2 Fold Change",
    y = "-Log10 (Adjusted P-value)"
  )

volcano_rectal_Endo_vs_Control

rectal_sig_CPP_vs_Control <- rectal_CPP_vs_Control_df %>%
  dplyr::filter(significance == "Significant") %>%
  dplyr::pull(KO)

rectal_sig_Endo_vs_Control <- rectal_Endo_vs_Control_df %>%
  dplyr::filter(significance == "Significant") %>%
  dplyr::pull(KO)

rectal_sig_Endo_vs_CPP <- rectal_Endo_vs_CPP_df %>%
  dplyr::filter(significance == "Significant") %>%
  dplyr::pull(KO)

rectal_unique_Endo_vs_CPP <- setdiff(
  rectal_sig_Endo_vs_CPP,
  union(rectal_sig_CPP_vs_Control, rectal_sig_Endo_vs_Control)
)

rectal_Endo_vs_CPP_df$plot_group <- "Not Significant"

rectal_Endo_vs_CPP_df$plot_group[
  rectal_Endo_vs_CPP_df$significance == "Significant"
] <- "Significant (Overlapping)"

rectal_Endo_vs_CPP_df$plot_group[
  rectal_Endo_vs_CPP_df$KO %in% rectal_unique_Endo_vs_CPP
] <- "Unique to Endo vs CPP"


# label top 10 UNIQUE ones in dark blue
top_unique_labels_Endo_vs_CPP <- rectal_Endo_vs_CPP_df %>%
  dplyr::filter(KO %in% rectal_unique_Endo_vs_CPP, !is.na(padj)) %>%
  dplyr::arrange(padj) %>%
  head(10)

# label top 10 overlapping significant ones in red
top_overlap_labels_Endo_vs_CPP <- rectal_Endo_vs_CPP_df %>%
  dplyr::filter(
    significance == "Significant",
    !(KO %in% rectal_unique_Endo_vs_CPP),
    !is.na(padj)
  ) %>%
  dplyr::arrange(padj) %>%
  head(10)

volcano_rectal_Endo_vs_CPP <- ggplot(
  rectal_Endo_vs_CPP_df,
  aes(x = log2FoldChange, y = -log10(padj), color = plot_group)
) +
  geom_point(alpha = 0.7) +
  geom_text_repel(
    data = top_unique_labels_Endo_vs_CPP,
    aes(label = KO),
    color = "blue",
    size = 3
  ) +
  geom_hline(yintercept = -log10(p_threshold), linetype = "dashed") +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed") +
  scale_color_manual(values = c(
    "Not Significant" = "grey70",
    "Significant (Overlapping)" = "red",
    "Unique to Endo vs CPP" = "blue"
  )) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 18),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 18),
    strip.background = element_rect(fill = "grey80", color = "black"),
    panel.border = element_rect(color = "black", fill = NA),
    plot.margin = margin(15, 15, 15, 15),
    legend.position = "none"
  ) +
  labs(
    x = "Log2 Fold Change",
    y = "-Log10 (Adjusted P-value)",
    color = NULL
  )

volcano_rectal_Endo_vs_CPP

#vaginal

res_vaginal_CPP_vs_Control <- results(
  dds_vaginal,
  contrast = c("Host_disease", "CPP", "Control")
)

res_vaginal_Endo_vs_Control <- results(
  dds_vaginal,
  contrast = c("Host_disease", "CPP Endo", "Control")
)

res_vaginal_Endo_vs_CPP <- results(
  dds_vaginal,
  contrast = c("Host_disease", "CPP Endo", "CPP")
)

vaginal_CPP_vs_Control_df <- as.data.frame(res_vaginal_CPP_vs_Control)
vaginal_CPP_vs_Control_df$KO <- rownames(vaginal_CPP_vs_Control_df)

vaginal_Endo_vs_Control_df <- as.data.frame(res_vaginal_Endo_vs_Control)
vaginal_Endo_vs_Control_df$KO <- rownames(vaginal_Endo_vs_Control_df)

vaginal_Endo_vs_CPP_df <- as.data.frame(res_vaginal_Endo_vs_CPP)
vaginal_Endo_vs_CPP_df$KO <- rownames(vaginal_Endo_vs_CPP_df)

res_vaginal_CPP_vs_Control <- results(
  dds_vaginal,
  contrast = c("Host_disease", "CPP", "Control")
)

res_vaginal_Endo_vs_Control <- results(
  dds_vaginal,
  contrast = c("Host_disease", "CPP Endo", "Control")
)

res_vaginal_Endo_vs_CPP <- results(
  dds_vaginal,
  contrast = c("Host_disease", "CPP Endo", "CPP")
)

vaginal_CPP_vs_Control_df <- as.data.frame(res_vaginal_CPP_vs_Control)
vaginal_CPP_vs_Control_df$KO <- rownames(vaginal_CPP_vs_Control_df)

vaginal_Endo_vs_Control_df <- as.data.frame(res_vaginal_Endo_vs_Control)
vaginal_Endo_vs_Control_df$KO <- rownames(vaginal_Endo_vs_Control_df)

vaginal_Endo_vs_CPP_df <- as.data.frame(res_vaginal_Endo_vs_CPP)
vaginal_Endo_vs_CPP_df$KO <- rownames(vaginal_Endo_vs_CPP_df)

vaginal_CPP_vs_Control_df$significance <- "Not Significant"
vaginal_CPP_vs_Control_df$significance[
  !is.na(vaginal_CPP_vs_Control_df$padj) &
    vaginal_CPP_vs_Control_df$padj < p_threshold &
    abs(vaginal_CPP_vs_Control_df$log2FoldChange) > logFC_threshold
] <- "Significant"

vaginal_Endo_vs_Control_df$significance <- "Not Significant"
vaginal_Endo_vs_Control_df$significance[
  !is.na(vaginal_Endo_vs_Control_df$padj) &
    vaginal_Endo_vs_Control_df$padj < p_threshold &
    abs(vaginal_Endo_vs_Control_df$log2FoldChange) > logFC_threshold
] <- "Significant"

vaginal_Endo_vs_CPP_df$significance <- "Not Significant"
vaginal_Endo_vs_CPP_df$significance[
  !is.na(vaginal_Endo_vs_CPP_df$padj) &
    vaginal_Endo_vs_CPP_df$padj < p_threshold &
    abs(vaginal_Endo_vs_CPP_df$log2FoldChange) > logFC_threshold
] <- "Significant"

# write df results as tsv
write_tsv(vaginal_Endo_vs_CPP_df |> filter(significance == 'Significant'),
          'results/aim3/06-cpp_cpp_endo_significant_KOs_vaginal.tsv')

write_tsv(vaginal_Endo_vs_Control_df |> filter(significance == 'Significant'),
          'results/aim3/06-cpp_endo_control_significant_KOs_vaginal.tsv')

write_tsv(vaginal_CPP_vs_Control_df |> filter(significance == 'Significant'),
          'results/aim3/06-cpp_control_significant_KOs_vaginal.tsv')


top_labels_vaginal_CPP_vs_Control <- vaginal_CPP_vs_Control_df %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::arrange(padj) %>%
  head(10)

volcano_vaginal_CPP_vs_Control <- ggplot(
  vaginal_CPP_vs_Control_df,
  aes(x = log2FoldChange, y = -log10(padj), color = significance)
) +
  geom_point(alpha = 0.7) +
  geom_text_repel(
    data = top_labels_vaginal_CPP_vs_Control,
    aes(label = KO),
    size = 3
  ) +
  geom_hline(yintercept = -log10(p_threshold), linetype = "dashed") +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed") +
  scale_color_manual(values = c("grey70", "red")) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 18),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 18),
    strip.background = element_rect(fill = "grey80", color = "black"),
    panel.border = element_rect(color = "black", fill = NA),
    plot.margin = margin(15, 15, 15, 15),
    legend.position = "none"
  ) +
  labs(
    x = "Log2 Fold Change",
    y = "-Log10 (Adjusted P-value)"
  )

top_labels_vaginal_Endo_vs_Control <- vaginal_Endo_vs_Control_df %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::arrange(padj) %>%
  head(10)

volcano_vaginal_Endo_vs_Control <- ggplot(
  vaginal_Endo_vs_Control_df,
  aes(x = log2FoldChange, y = -log10(padj), color = significance)
) +
  geom_point(alpha = 0.7) +
  geom_text_repel(
    data = top_labels_vaginal_Endo_vs_Control,
    aes(label = KO),
    size = 3
  ) +
  geom_hline(yintercept = -log10(p_threshold), linetype = "dashed") +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed") +
  scale_color_manual(values = c("grey70", "red")) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 18),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 18),
    strip.background = element_rect(fill = "grey80", color = "black"),
    panel.border = element_rect(color = "black", fill = NA),
    plot.margin = margin(15, 15, 15, 15),
    legend.position = "none"
  ) +
  labs(
    x = "Log2 Fold Change",
    y = "-Log10 (Adjusted P-value)"
  )

vaginal_sig_CPP_vs_Control <- vaginal_CPP_vs_Control_df %>%
  dplyr::filter(significance == "Significant") %>%
  dplyr::pull(KO)

vaginal_sig_Endo_vs_Control <- vaginal_Endo_vs_Control_df %>%
  dplyr::filter(significance == "Significant") %>%
  dplyr::pull(KO)

vaginal_sig_Endo_vs_CPP <- vaginal_Endo_vs_CPP_df %>%
  dplyr::filter(significance == "Significant") %>%
  dplyr::pull(KO)

vaginal_unique_Endo_vs_CPP <- setdiff(
  vaginal_sig_Endo_vs_CPP,
  union(vaginal_sig_CPP_vs_Control, vaginal_sig_Endo_vs_Control)
)

vaginal_Endo_vs_CPP_df$plot_group <- "Not Significant"

vaginal_Endo_vs_CPP_df$plot_group[
  vaginal_Endo_vs_CPP_df$significance == "Significant"
] <- "Significant (Overlapping)"

vaginal_Endo_vs_CPP_df$plot_group[
  vaginal_Endo_vs_CPP_df$KO %in% vaginal_unique_Endo_vs_CPP
] <- "Unique to Endo vs CPP"


# label all UNIQUE ones in dark blue
top_unique_labels_vaginal_Endo_vs_CPP <- vaginal_Endo_vs_CPP_df %>%
  dplyr::filter(
    KO %in% vaginal_unique_Endo_vs_CPP,
    !is.na(padj)
  )

# label top 10 overlapping significant ones in red
top_overlap_labels_vaginal_Endo_vs_CPP <- vaginal_Endo_vs_CPP_df %>%
  dplyr::filter(
    significance == "Significant",
    !(KO %in% vaginal_unique_Endo_vs_CPP),
    !is.na(padj)
  ) %>%
  dplyr::arrange(padj) %>%
  head(10)

volcano_vaginal_Endo_vs_CPP <- ggplot(
  vaginal_Endo_vs_CPP_df,
  aes(x = log2FoldChange, y = -log10(padj), color = plot_group)
) +
  geom_point(alpha = 0.7) +
  geom_text_repel(
    data = top_unique_labels_vaginal_Endo_vs_CPP,
    aes(label = KO),
    color = "blue",
    size = 3
  ) +
  geom_hline(yintercept = -log10(p_threshold), linetype = "dashed") +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed") +
  scale_color_manual(values = c(
    "Not Significant" = "grey70",
    "Significant (Overlapping)" = "red",
    "Unique to Endo vs CPP" = "blue"
  )) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 18),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 18),
    strip.background = element_rect(fill = "grey80", color = "black"),
    panel.border = element_rect(color = "black", fill = NA),
    plot.margin = margin(15, 15, 15, 15),
    legend.position = "none"
  ) +
  labs(
    x = "Log2 Fold Change",
    y = "-Log10 (Adjusted P-value)",
    color = NULL
  )

volcano_vaginal_Endo_vs_CPP
# Combine rectal + vaginal panels per contrast
panel_CPP_vs_Control <- volcano_rectal_CPP_vs_Control + volcano_vaginal_CPP_vs_Control 

panel_Endo_vs_Control <- volcano_rectal_Endo_vs_Control + volcano_vaginal_Endo_vs_Control 
#  plot_annotation(title = "KO Differential Abundance (CPP Endo vs Control)") &
#  theme(plot.title = element_text(size = 25))

panel_Endo_vs_CPP <- volcano_rectal_Endo_vs_CPP + volcano_vaginal_Endo_vs_CPP 

panel_CPP_vs_Control
panel_Endo_vs_Control
panel_Endo_vs_CPP

ggsave(
  "results/aim3/CPP_vs_Control_volcano.png",
  plot = panel_CPP_vs_Control,
  width = 10, height = 7, units = "in", dpi = 300
)

ggsave(
  "results/aim3/Endo_vs_Control_volcano.png",
  plot = panel_Endo_vs_Control,
  width = 10, height = 7, units = "in", dpi = 300
)

ggsave(
  "results/aim3/Endo_vs_CPP_volcano.png",
  plot = panel_Endo_vs_CPP,
  width = 10, height = 7, units = "in", dpi = 300
)

##########Bar###############
df <- data.frame(
  pathway = c(
    "Metabolic pathways",
    "ABC transporters",
    "Lipopolysaccharide biosynthesis",
    "Biosynthesis of unsaturated fatty acids",
    "Cationic antimicrobial peptide (CAMP) resistance",
    "Bacterial secretion system",
    "Biosynthesis of cofactors",
    "Thiamine metabolism",
    "Protein export"
  ),
  count = c(2, 2, 1, 1, 1, 1, 1, 1, 1)
)

df <- df %>% arrange(count)

bar <- ggplot(df, aes(x = count, y = reorder(pathway, count), fill = count)) +
  geom_col() +
  geom_text(aes(label = count), hjust = -0.15, size = 4) +
  scale_fill_gradient(low = "#1B98E026", high="#225ea8") +
  labs(
    x = "Number of associated enriched KOs",
    y = NULL,
    fill = "Count"
  ) +
  theme_classic(base_size = 12) +
  xlim(0, max(df$count) + 0.5)

ggsave(
  "results/aim3/bar.png",
  plot = bar,
  width = 10, height = 5, units = "in", dpi = 300
)
bar