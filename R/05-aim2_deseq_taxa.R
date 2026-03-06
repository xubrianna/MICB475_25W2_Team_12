#load packages
library(phyloseq)
library(DESeq2)
library(dplyr)
library(ggplot2)

#load and preprocess
phyloseq <- readRDS("data/phyloseq_filtered.rds")
phyloseq <- prune_taxa(taxa_sums(phyloseq) > 0, phyloseq)
prev_threshold <- 0.1 * nsamples(phyloseq)
phyloseq_filtered <- prune_taxa(apply(otu_table(phyloseq), 1, function(x) sum(x > 0)) > prev_threshold, phyloseq)
phyloseq_genus_counts <- tax_glom(phyloseq_filtered, taxrank = "Genus")
phyloseq_genus_counts <- prune_taxa(taxa_sums(phyloseq_genus_counts) > 5, phyloseq_genus_counts)

#DESeq2 wrapper for a comparison
run_deseq_pair <- function(phylo_site, group1, group2, site_name){
  phylo_site <- prune_taxa(taxa_sums(phylo_site) > 0, phylo_site)
  sample_data(phylo_site)$Host_disease <- make.names(sample_data(phylo_site)$Host_disease)
  dds <- phyloseq_to_deseq2(phylo_site, ~ Host_disease)
  dds <- estimateSizeFactors(dds, type="poscounts")
  dds <- DESeq(dds)
  
  res <- results(dds, contrast = c("Host_disease", make.names(group1), make.names(group2)))
  res_df <- as.data.frame(res)
  res_df$Genus <- gsub("^g_+", "", tax_table(phylo_site)[rownames(res_df), "Genus"])
  res_df <- res_df[!is.na(res_df$Genus), ]
  res_df <- res_df[order(res_df$padj), ]
  return(res_df)
}

#volcano plot in ALDEx2 style
plot_volcano_fixed_shape <- function(res_df, site_name, lfc_thresh = 1, padj_thresh = 0.05,
                                     x_lim = NULL, y_lim = NULL) {
  
  res_df <- res_df %>% filter(!is.na(padj)) %>%
    mutate(Significant = ifelse(padj <= padj_thresh & abs(log2FoldChange) >= lfc_thresh,
                                "Significant","Not_Significant"))
  res_df$Significant <- factor(res_df$Significant, levels = c("Not_Significant","Significant"))
  res_df$negLog10P <- -log10(res_df$padj)
  if(is.null(x_lim)) x_lim <- max(abs(res_df$log2FoldChange)) * 1.05
  if(is.null(y_lim)) y_lim <- max(res_df$negLog10P) * 1.05
  
  ggplot(res_df, aes(x = log2FoldChange, y = negLog10P, color = Significant)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), color = "red", linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_thresh), color = "red", linetype = "dashed") +
    scale_color_manual(values = c("Significant"="red","Not_Significant"="gray")) +
    theme_bw() +
    labs(title = paste0("Differentially Abundant Taxa (", site_name, ")"),
         x = "Log2 Fold Change",
         y = "-log10 Adjusted p-value") +
    theme(plot.title = element_text(face="bold", size=14),
          legend.position = "none") +
    scale_x_continuous(limits = c(-x_lim, x_lim)) +
    scale_y_continuous(limits = c(0, y_lim))
}

#rectal comparisons
phy_rectal <- subset_samples(phyloseq_genus_counts, collection_method=="rectal_swab")
res_cpp_vs_endo <- run_deseq_pair(phy_rectal, "CPP Endo", "CPP", "Rectal_CPP_vs_CPPEndo")
res_cpp_vs_ctrl <- run_deseq_pair(phy_rectal, "CPP", "Control", "Rectal_CPP_vs_Control")
res_endo_vs_ctrl <- run_deseq_pair(phy_rectal, "CPP Endo", "Control", "Rectal_CPPEndo_vs_Control")

all_res <- bind_rows(res_cpp_vs_endo, res_cpp_vs_ctrl, res_endo_vs_ctrl)
global_x_lim <- max(abs(all_res$log2FoldChange)) * 1.05
global_y_lim <- max(-log10(all_res$padj), na.rm = TRUE) * 1.05

volcano_cpp_vs_endo <- plot_volcano_fixed_shape(res_cpp_vs_endo, "Rectal CPP Endo vs CPP",
                                                x_lim = global_x_lim, y_lim = global_y_lim)
volcano_cpp_vs_ctrl <- plot_volcano_fixed_shape(res_cpp_vs_ctrl, "Rectal CPP vs Control",
                                                x_lim = global_x_lim, y_lim = global_y_lim)
volcano_endo_vs_ctrl <- plot_volcano_fixed_shape(res_endo_vs_ctrl, "Rectal CPP Endo vs Control",
                                                 x_lim = global_x_lim, y_lim = global_y_lim)

#print
volcano_cpp_vs_endo
volcano_cpp_vs_ctrl
volcano_endo_vs_ctrl

ggsave("results/aim2/volcano_rectal_cpp_vs_endo.png", volcano_cpp_vs_endo, width=8, height=6)
ggsave("results/aim2/volcano_rectal_cpp_vs_ctrl.png", volcano_cpp_vs_ctrl, width=8, height=6)
ggsave("results/aim2/volcano_rectal_endo_vs_ctrl.png", volcano_endo_vs_ctrl, width=8, height=6)

# Vaginal comparisons using fixed-shape volcano plots
phy_vaginal <- subset_samples(phyloseq_genus_counts, collection_method=="vaginal_swab")
res_cpp_vs_endo_vag <- run_deseq_pair(phy_vaginal, "CPP Endo", "CPP", "Vaginal_CPP_vs_CPPEndo")
res_cpp_vs_ctrl_vag <- run_deseq_pair(phy_vaginal, "CPP", "Control", "Vaginal_CPP_vs_Control")
res_endo_vs_ctrl_vag <- run_deseq_pair(phy_vaginal, "CPP Endo", "Control", "Vaginal_CPPEndo_vs_Control")

# Compute global limits across vaginal comparisons
all_res_vag <- bind_rows(res_cpp_vs_endo_vag, res_cpp_vs_ctrl_vag, res_endo_vs_ctrl_vag)
global_x_lim_vag <- max(abs(all_res_vag$log2FoldChange)) * 1.05
global_y_lim_vag <- max(-log10(all_res_vag$padj), na.rm = TRUE) * 1.05

#generate fixed-shape volcano plots
volcano_cpp_vs_endo_vag <- plot_volcano_fixed_shape(res_cpp_vs_endo_vag, "Vaginal CPP Endo vs CPP",
                                                    x_lim = global_x_lim_vag, y_lim = global_y_lim_vag)
volcano_cpp_vs_ctrl_vag <- plot_volcano_fixed_shape(res_cpp_vs_ctrl_vag, "Vaginal CPP vs Control",
                                                    x_lim = global_x_lim_vag, y_lim = global_y_lim_vag)
volcano_endo_vs_ctrl_vag <- plot_volcano_fixed_shape(res_endo_vs_ctrl_vag, "Vaginal CPP Endo vs Control",
                                                     x_lim = global_x_lim_vag, y_lim = global_y_lim_vag)
#print plots
volcano_cpp_vs_endo_vag
volcano_cpp_vs_ctrl_vag
volcano_endo_vs_ctrl_vag

#save plots
ggsave("results/aim2/volcano_vaginal_cpp_vs_endo.png", volcano_cpp_vs_endo_vag, width=8, height=6)
ggsave("results/aim2/volcano_vaginal_cpp_vs_ctrl.png", volcano_cpp_vs_ctrl_vag, width=8, height=6)
ggsave("results/aim2/volcano_vaginal_endo_vs_ctrl.png", volcano_endo_vs_ctrl_vag, width=8, height=6)

#thresholds
padj_thresh <- 0.05
lfc_thresh <- 1

# Rectal CPP Endo vs CPP
sig_rectal <- res_cpp_vs_endo %>%
  filter(!is.na(padj)) %>%
  filter(padj <= padj_thresh & abs(log2FoldChange) >= lfc_thresh) %>%
  arrange(padj)
print("Rectal CPP Endo vs CPP:")
print(sig_rectal$Genus)

# Rectal CPP vs Control
sig_cpp_vs_ctrl <- res_cpp_vs_ctrl %>%
  filter(!is.na(padj)) %>%
  filter(padj <= padj_thresh & abs(log2FoldChange) >= lfc_thresh) %>%
  arrange(padj)
print("Rectal CPP vs Control:")
print(sig_cpp_vs_ctrl$Genus)

# Rectal CPP Endo vs Control
sig_endo_vs_ctrl <- res_endo_vs_ctrl %>%
  filter(!is.na(padj)) %>%
  filter(padj <= padj_thresh & abs(log2FoldChange) >= lfc_thresh) %>%
  arrange(padj)
print("Rectal CPP Endo vs Control:")
print(sig_endo_vs_ctrl$Genus)

# Vaginal CPP Endo vs CPP
sig_vag_cppendo_vs_cpp <- res_cpp_vs_endo_vag %>%
  filter(!is.na(padj)) %>%
  filter(padj <= padj_thresh & abs(log2FoldChange) >= lfc_thresh) %>%
  arrange(padj)
print("Vaginal CPP Endo vs CPP:")
print(sig_vag_cppendo_vs_cpp$Genus)

# Vaginal CPP vs Control
sig_vag_cpp_vs_ctrl <- res_cpp_vs_ctrl_vag %>%
  filter(!is.na(padj)) %>%
  filter(padj <= padj_thresh & abs(log2FoldChange) >= lfc_thresh) %>%
  arrange(padj)
print("Vaginal CPP vs Control:")
print(sig_vag_cpp_vs_ctrl$Genus)

# Vaginal CPP Endo vs Control
sig_vag_endo_vs_ctrl <- res_endo_vs_ctrl_vag %>%
  filter(!is.na(padj)) %>%
  filter(padj <= padj_thresh & abs(log2FoldChange) >= lfc_thresh) %>%
  arrange(padj)
print("Vaginal CPP Endo vs Control:")
print(sig_vag_endo_vs_ctrl$Genus)

