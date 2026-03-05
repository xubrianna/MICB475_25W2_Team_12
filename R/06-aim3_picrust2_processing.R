install.packages("devtools")
devtools::install_github("cafferychen777/ggpicrust2")
install.packages("MicrobiomeStat")
install.packages("GGally")
install.packages("phyloseq")

library(tidyverse)
library(phyloseq)
library(ggpicrust2)
library(DESeq2)

ps <- readRDS("data/phyloseq_filtered.rds")

meta <- sample_data(ps) %>%
  data.frame() %>%
  rownames_to_column("sample_name")
ko <- read.delim("data/pred_metagenome_unstrat.tsv", row.names = 1)
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

sig_rectal <- subset(res_rectal_df, !is.na(padj) & padj < 0.05)
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


sig_vaginal <- subset(res_vaginal_df, !is.na(padj) & padj < 0.05)
nrow(sig_vaginal)
head(sig_vaginal, 10)