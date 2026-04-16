library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)
library(cowplot)
library(ggpubr)

####LOAD AND MODIFY PHYLOSEQ####
#load filtered phyloseq without rarification (expects uneven sequencing depths)
phyloseq <- readRDS('data/phyloseq_filtered.rds')

# Transform OTU table into relative abundance from absolute counts
CPP_RA <- transform_sample_counts(phyloseq, fun = function(x) x / sum(x))

#### STRATIFY BY HOST DISEASE ####
CPP_only <- subset_samples(CPP_RA, `Host_disease` == "CPP Only")
CPP_endo <- subset_samples(CPP_RA, `Host_disease` == "CPP Endo")
CPP_healthy <- subset_samples(CPP_RA, `Host_disease` == "Control")

#### RECTAL CORE MICROBIOME ####
# Filter by rectal samples
CPP_only_rectal <- subset_samples(CPP_only, `collection_method` == "rectal_swab")
CPP_endo_rectal <- subset_samples(CPP_endo, `collection_method` == "rectal_swab")
CPP_healthy_rectal <- subset_samples(CPP_healthy, `collection_method` == "rectal_swab")

# Detection 0.001 (filters out rare taxa) , prevalence 0.5
CPP_only_ASVs_rectal_0.001 <- core_members(CPP_only_rectal, detection = 0.001, prevalence = 0.5)
CPP_endo_ASVs_rectal_0.001 <- core_members(CPP_endo_rectal, detection = 0.001, prevalence = 0.5)
CPP_healthy_ASVs_rectal_0.001 <- core_members(CPP_healthy_rectal, detection = 0.001, prevalence = 0.5)

CPP_list_rectal_0.001 <- list(
  "Control" = CPP_healthy_ASVs_rectal_0.001,
  "CPP Only" = CPP_only_ASVs_rectal_0.001,
  "CPP Endo" = CPP_endo_ASVs_rectal_0.001
)

venn_all_diseases_rectal_0.001 <- ggVennDiagram(
  x = CPP_list_rectal_0.001,
  category.names = c("Control", "CPP Only", "CPP Endo"),
  label_alpha = 0) + 
  scale_fill_gradient(low = "#1B98E026", high = "#225ea8") +
  theme(legend.position = 'none') 

venn_all_diseases_rectal_0.001 <- ggVennDiagram(
  x = CPP_list_rectal_0.001,
  category.names = c("Control", "CPP Only", "CPP Endo"),
  label_alpha = 0
) +
  theme(legend.position = "none")
  # scale_fill_manual(c("Control" = "#c7e9b4", "CPP Only" = "#41b6c4", "CPP Endo" = "#225ea8"))
venn_all_diseases_rectal_0.001

ggsave("results/aim2/04-core_microbiome/venn_all_diseases_rectal_0.001.png",
       venn_all_diseases_rectal_0.001, width = 7, height = 7)

#get the taxonomic information for CPP_endo disease group and find the ASVs that are unique to the core of CPP-endo 
tax_mat_rectal <- tax_table(CPP_endo_rectal)

core_taxonomy_rectal <- as.data.frame(tax_mat_rectal[CPP_endo_ASVs_rectal_0.001, ])
print(head(core_taxonomy_rectal))

unique_to_CPP_endo_rectal<- setdiff(CPP_endo_ASVs_rectal_0.001, union(CPP_only_ASVs_rectal_0.001, CPP_healthy_ASVs_rectal_0.001))
print(unique_to_CPP_endo_rectal)

#### VAGINAL COREMICROBIOME #### 

#stratifying by the vaginal sample instead now  
CPP_only_vaginal <- subset_samples(CPP_only,`collection_method` == "vaginal_swab")
CPP_endo_vaginal <- subset_samples(CPP_endo,`collection_method` == "vaginal_swab")
CPP_healthy_vaginal <-subset_samples(CPP_healthy, `collection_method` == "vaginal_swab")

CPP_only_ASVs_vaginal_0.001 <- core_members(CPP_only_vaginal, detection=0.001, prevalence = 0.5)
CPP_endo_ASVs_vaginal_0.001 <- core_members(CPP_endo_vaginal, detection=0.001, prevalence = 0.5)
CPP_healthy_ASVs_vaginal_0.001 <- core_members(CPP_healthy_vaginal, detection=0.001, prevalence = 0.5)

CPP_list_vaginal_0.001 <- list(Control = CPP_healthy_ASVs_vaginal_0.001, CPP_Only = CPP_only_ASVs_vaginal_0.001, CPP_Endo = CPP_endo_ASVs_vaginal_0.001)

venn_all_diseases_vaginal_0.001 <- ggVennDiagram(x = CPP_list_vaginal_0.001, category.names = c("Control", "CPP Only", "CPP Endo"), label_alpha = 0)+ scale_fill_gradient(low = "#1B98E026", high="#225ea8") + 
  theme(legend.position = 'none')
venn_all_diseases_vaginal_0.001

ggsave("results/aim2/04-core_microbiome/venn_all_diseases_vaginal_0.001.png", venn_all_diseases_vaginal_0.001, width = 7, height = 7)


#### MAKE PLOTS WITH SAME FILL SCALE

# build processed venn objects
rectal_pd  <- process_data(Venn(CPP_list_rectal_0.001))
vaginal_pd <- process_data(Venn(CPP_list_vaginal_0.001))

# get shared max count across both
shared_max <- max(
  rectal_pd$regionData$count,
  vaginal_pd$regionData$count
)

venn_all_diseases_rectal_0.001 <- ggVennDiagram(
  x = CPP_list_rectal_0.001,
  category.names = c("Control", "CPP Only", "CPP Endo"),
  label_alpha = 0
) +
  scale_fill_gradient(
    low = "#1B98E026",
    high = "#225ea8",
    limits = c(0, shared_max)
  ) +
  theme(legend.position = "none",
        plot.margin = margin(t = 0, 
                             r = 15,
                             b = 15,  
                             l = 15,  
                             unit = "pt"
        ))

venn_all_diseases_vaginal_0.001 <- ggVennDiagram(
  x = CPP_list_vaginal_0.001,
  category.names = c("Control", "CPP Only", "CPP Endo"),
  label_alpha = 0
) +
  scale_fill_gradient(
    low = "#1B98E026",
    high = "#225ea8",
    limits = c(0, shared_max),
    name = 'Count'
  ) +
  theme(legend.position = 'none',
        plot.margin = margin(t = 0, 
                             r = 15,
                             b = 15,  
                             l = 15,  
                             unit = "pt"))
  
#prune 
prune_taxa(CPP_only_ASVs_rectal_0.001,CPP_RA) %>%
  tax_table()

prune_taxa(CPP_endo_ASVs_rectal_0.001,CPP_RA) %>%
  tax_table()

prune_taxa(CPP_healthy_ASVs_rectal_0.001,CPP_RA) %>%
  tax_table()

#### MEAN RELATIVE ABUNDANCE PLOTS FOR CORE TAXA ####

# Helper: extract mean RA data from pruned phyloseq
get_mean_ra_df <- function(asv_list, ps) {
  ps_pruned <- prune_taxa(asv_list, ps)
  df <- psmelt(ps_pruned)
  df$Genus <- gsub("^g__", "", df$Genus)
  df %>%
    group_by(Host_disease, Genus) %>%
    summarise(mean_Abundance = mean(Abundance), .groups = "drop")
}

# Helper: extract RA data from pruned phyloseq
get_ra_df <- function(asv_list, ps) {
  ps_pruned <- prune_taxa(asv_list, ps)
  df <- psmelt(ps_pruned)
  df$Genus <- gsub("^g__", "", df$Genus)
  df %>%
    group_by(Host_disease, Genus) %>%
    select(Host_disease,Genus, Sample, Abundance)
    
}

# Get mean RA data for all 6 plots
df_rectal_healthy  <- get_mean_ra_df(CPP_healthy_ASVs_rectal_0.001, CPP_RA)
df_rectal_cpp_only <- get_mean_ra_df(CPP_only_ASVs_rectal_0.001, CPP_RA)
df_rectal_cpp_endo <- get_mean_ra_df(CPP_endo_ASVs_rectal_0.001, CPP_RA)


# Get RA data for all 6 plots
df_rectal_healthy_sample  <- get_ra_df(CPP_healthy_ASVs_rectal_0.001, CPP_RA) # %>% filter(Host_disease == 'Control')
df_rectal_cpp_only_sample  <- get_ra_df(CPP_only_ASVs_rectal_0.001, CPP_RA)  # %>% filter(Host_disease == 'CPP Only')
df_rectal_cpp_endo_sample  <- get_ra_df(CPP_endo_ASVs_rectal_0.001, CPP_RA) # %>% filter(Host_disease == 'CPP Endo')


#prune vaginal
prune_taxa(CPP_only_ASVs_vaginal_0.001,CPP_RA) %>%
  tax_table()

prune_taxa(CPP_endo_ASVs_vaginal_0.001,CPP_RA) %>%
  tax_table()

prune_taxa(CPP_healthy_ASVs_vaginal_0.001,CPP_RA) %>%
  tax_table()

df_vaginal_healthy  <- get_mean_ra_df(CPP_healthy_ASVs_vaginal_0.001, CPP_RA)
df_vaginal_cpp_only <- get_mean_ra_df(CPP_only_ASVs_vaginal_0.001, CPP_RA)
df_vaginal_cpp_endo <- get_mean_ra_df(CPP_endo_ASVs_vaginal_0.001, CPP_RA)

df_vaginal_healthy_sample  <- get_ra_df(CPP_healthy_ASVs_vaginal_0.001, CPP_RA) # %>% filter(Host_disease == 'Control')
df_vaginal_cpp_only_sample  <- get_ra_df(CPP_only_ASVs_vaginal_0.001, CPP_RA) # %>% filter(Host_disease == 'CPP Only')
df_vaginal_cpp_endo_sample  <- get_ra_df(CPP_endo_ASVs_vaginal_0.001, CPP_RA)#  %>% filter(Host_disease == 'CPP Endo')

# Consistent genus color palette across all plots
all_genera <- unique(c(
  df_rectal_healthy$Genus, df_rectal_cpp_only$Genus, df_rectal_cpp_endo$Genus,
  df_vaginal_healthy$Genus, df_vaginal_cpp_only$Genus, df_vaginal_cpp_endo$Genus
))

genus_colors <- setNames(
  colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(length(all_genera)),
  sort(all_genera)
)

# Helper: make mean RA stacked bar plot
make_core_plot <- function(df, title = "") {
  df$Host_disease <- factor(df$Host_disease, levels = c("Control", "CPP Only", "CPP Endo"))
  df$Genus <- factor(df$Genus, levels = sort(all_genera))
  ggplot(df, aes(x = Host_disease, y = mean_Abundance, fill = Genus)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = genus_colors) +
    theme_classic() +
    labs(x = "Host Disease", y = "Mean Relative Abundance", title = title) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.margin = margin(15,15,15,15,))
}

make_core_plot_by_sample <- function(df, title = "") {
  df$Host_disease <- factor(df$Host_disease, levels = c("Control", "CPP Only", "CPP Endo"))
  df$Genus <- factor(df$Genus, levels = sort(all_genera))
  ggplot(df, aes(x = Sample, y = Abundance, fill = Genus)) +
    facet_wrap(~Host_disease, scales = "free_x" )+
    geom_bar(stat = "identity") +
    scale_fill_manual(values = genus_colors) +
    theme_classic() +
    labs(x = "Sample", y = "Relative Abundance", title = title) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          plot.margin = margin(15,15,15,15,))
}

# Create all 6 plots
prune_taxa_rectal_CPP_healthy_0.001 <- make_core_plot(df_rectal_healthy)
prune_taxa_rectal_CPP_only_0.001    <- make_core_plot(df_rectal_cpp_only)
prune_taxa_rectal_CPP_endo_0.001    <- make_core_plot(df_rectal_cpp_endo)
prune_taxa_vaginal_CPP_healthy_0.001 <- make_core_plot(df_vaginal_healthy)
prune_taxa_vaginal_CPP_only_0.001    <- make_core_plot(df_vaginal_cpp_only)
prune_taxa_vaginal_CPP_endo_0.001    <- make_core_plot(df_vaginal_cpp_endo)

# Create all 6 plots
prune_taxa_rectal_CPP_healthy_0.001_sample <- make_core_plot_by_sample(df_rectal_healthy_sample)
prune_taxa_rectal_CPP_only_0.001_sample    <- make_core_plot_by_sample(df_rectal_cpp_only_sample)
prune_taxa_rectal_CPP_endo_0.001_sample    <- make_core_plot_by_sample(df_rectal_cpp_endo_sample)
prune_taxa_vaginal_CPP_healthy_0.001_sample <- make_core_plot_by_sample(df_vaginal_healthy_sample)
prune_taxa_vaginal_CPP_only_0.001_sample    <- make_core_plot_by_sample(df_vaginal_cpp_only_sample)
prune_taxa_vaginal_CPP_endo_0.001_sample    <- make_core_plot_by_sample(df_vaginal_cpp_endo_sample)

ggsave("results/aim2/04-core_microbiome/mean_prune_taxa_rectal_CPP_only_0.001.png", prune_taxa_rectal_CPP_only_0.001, width = 7, height = 7)
ggsave("results/aim2/04-core_microbiome/mean_prune_taxa_rectal_CPP_endo_0.001.png", prune_taxa_rectal_CPP_endo_0.001, width = 7, height = 7)
ggsave("results/aim2/04-core_microbiome/mean_prune_taxa_rectal_CPP_healthy_0.001.png", prune_taxa_rectal_CPP_healthy_0.001, width = 7, height = 7)
ggsave("results/aim2/04-core_microbiome/mean_prune_taxa_vaginal_CPP_only_0.001.png", prune_taxa_vaginal_CPP_only_0.001, width = 7, height = 7)
ggsave("results/aim2/04-core_microbiome/mean_prune_taxa_vaginal_CPP_endo_0.001.png", prune_taxa_vaginal_CPP_endo_0.001, width = 7, height = 7)
ggsave("results/aim2/04-core_microbiome/mean_prune_taxa_vaginal_CPP_healthy_0.001.png", prune_taxa_vaginal_CPP_healthy_0.001, width = 7, height = 7)

#### COMBINED PRUNE TAXA PANEL ####
prune_taxa_panel_sample <- plot_grid(
  prune_taxa_rectal_CPP_healthy_0.001_sample, prune_taxa_vaginal_CPP_healthy_0.001_sample,
  prune_taxa_rectal_CPP_only_0.001_sample, prune_taxa_vaginal_CPP_only_0.001_sample,
  prune_taxa_rectal_CPP_endo_0.001_sample, prune_taxa_vaginal_CPP_endo_0.001_sample,
  ncol = 2, nrow = 3,
  labels = c("A", "B", "C", "D", "E", "F"),
  label_size = 16,
  rel_widths = c(1.25, 1)
)
prune_taxa_panel_sample
ggsave("results/aim2/04-core_microbiome/prune_taxa_panel_sample.png",
       prune_taxa_panel_sample, width = 14, height = 15, dpi = 300)

prune_taxa_panel <- plot_grid(
  prune_taxa_rectal_CPP_healthy_0.001, prune_taxa_vaginal_CPP_healthy_0.001,
  prune_taxa_rectal_CPP_only_0.001, prune_taxa_vaginal_CPP_only_0.001,
  prune_taxa_rectal_CPP_endo_0.001, prune_taxa_vaginal_CPP_endo_0.001,
  ncol = 2, nrow = 3,
  labels = c("A", "B", "C", "D", "E", "F"),
  label_size = 16,
  rel_widths = c(1.25, 1)
)
prune_taxa_panel
ggsave("results/aim2/04-core_microbiome/prune_taxa_panel.png",
       prune_taxa_panel, width = 14, height = 15, dpi = 300)

#### CAMPYLOBACTER RELATIVE ABUNDANCE ####
# Subset to rectal samples and transform to relative abundance
phylo_rectal <- subset_samples(phyloseq, collection_method == "rectal_swab")
phylo_rectal_RA <- transform_sample_counts(phylo_rectal, function(x) x / sum(x))

# Identify all ASVs classified as Campylobacter at genus level
tax <- as.data.frame(tax_table(phylo_rectal_RA))
campy_asvs <- rownames(tax[grepl("Campylobacter", tax$Genus, ignore.case = TRUE), ])

# Sum relative abundance of all Campylobacter ASVs per sample
otu_mat <- as.data.frame(as.matrix(otu_table(phylo_rectal_RA)))
campy_abund <- colSums(otu_mat[campy_asvs, , drop = FALSE])

# Build data frame with disease group info
campy_df <- data.frame(
  Sample = names(campy_abund),
  Campy_RA = campy_abund
) %>%
  left_join(
    as.data.frame(sample_data(phylo_rectal_RA)) %>%
      rownames_to_column("Sample"),
    by = "Sample"
  ) %>%
  filter(!is.na(Host_disease))

## add significances
comparisons <- list(
  c("Control", "CPP Only"),
  c("Control", "CPP Endo"),
  c("CPP Only", "CPP Endo")
)

# Mean relative abundance per disease group
campy_summary <- campy_df %>%
  group_by(Host_disease) %>%
  summarise(Mean_RA = mean(Campy_RA), .groups = "drop")
print(campy_summary)

campy_df$Host_disease <- factor(campy_df$Host_disease, levels = c("Control", "CPP Only", "CPP Endo"))
# Plot
plot_campylobacter <- ggplot(campy_df, aes(x = Host_disease, y = Campy_RA, fill = Host_disease)) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2.5, color = "black") +
  scale_fill_manual(values = c("Control" = "#c7e9b4", "CPP Only" = "#41b6c4", "CPP Endo" = "#225ea8")) +
  theme_classic() +
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    label = "p.signif",
    hide.ns = FALSE
  ) +
  ylab("Campylobacter Relative Abundance") +
  xlab("") +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 18),
    plot.margin = margin(15, 15, 15, 15),
    axis.title.x = element_text(size = 18, margin = margin(t = 20)),
    axis.title.y = element_text(size = 18, margin = margin(r = 20)),
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", colour = NA),
    legend.background = element_rect(fill = "transparent"),
    legend.box.background = element_rect(fill = "transparent")
  )

ggsave("results/aim2/04-core_microbiome/plot_campylobacter.png",
       plot_campylobacter, width = 7, height = 7)



###################################
##### panel for final plot



venn_all_diseases_vaginal_0.001 <- venn_all_diseases_vaginal_0.001 + ggtitle("Vaginal") + theme(plot.title = element_text(hjust = 0.5, size= 18))
venn_all_diseases_rectal_0.001 <- venn_all_diseases_rectal_0.001 + ggtitle("Rectal") + theme(plot.title = element_text(hjust = 0.5,  size= 18))


venn_legend <- get_legend(venn_all_diseases_vaginal_0.001 + theme(legend.position = 'right'))


# Put vaginal plot + legend together, controlling widths
venn_vaginal_with_legend <- plot_grid(
  venn_vaginal, venn_legend,
  ncol = 2,
  rel_widths = c(1, 0.25)
)

p <- plot_grid(
  venn_all_diseases_rectal_0.001,
  venn_all_diseases_vaginal_0.001,
  venn_legend,
  labels = c("", "", ""),
  ncol = 3,
  label_size = 20,
  rel_widths = c(1, 1, 0.25)
)

p <- plot_grid(p, plot_campylobacter, labels = c("A", "B"),
               ncol = 1, label_size = 20)

ggsave(
  "results/aim2/04-core_microbiome/04-final_figure.png",
  plot = p,
  width = 12, height = 13, units = "in", dpi = 300)


### March 28 trying to find what is unqiue to healthy in rectal coremicrobiome
#get the taxonomic information for CPP_endo disease group and find the ASVs that are unique to the core of CPP-endo 
tax_mat_rectal_2 <- tax_table(CPP_endo_rectal)

core_taxonomy_rectal_2 <- as.data.frame(tax_mat_rectal_2[CPP_healthy_ASVs_rectal_0.001, ])
print(head(core_taxonomy_rectal_2))
      
unique_to_CPP_endo_rectal_2<- setdiff(CPP_healthy_ASVs_rectal_0.001, union(CPP_only_ASVs_rectal_0.001, CPP_endo_ASVs_rectal_0.001))
print(unique_to_CPP_endo_rectal_2)


##April 3 - trying to find the relative abundance of ASVs in all of the samples that are campylobacter 
#try filtering the full taxatable from the phyloseq object to be campylobacter genus only just to see how many there are
filtered_phyloseq_campylobacter <- subset_taxa(phyloseq, Genus == "g__Campylobacter")
ntaxa(filtered_phyloseq_campylobacter)


#try pruning 
ps_A <- prune_taxa(taxa_sums(ps_A) > 0, ps_A)

#try filtering the phyloseq to only have campylobacter that are in that genus = just for cpp_only rectal 
CPP_only_rectal_campylobacter_genus <- subset_taxa(CPP_only_rectal_campylobacter, Genus == "g__Campylobacter")
#count the number of total ASVs in the group 
ntaxa(CPP_only_rectal_campylobacter)
#count the number of ASvs that are campylobacter
ntaxa(CPP_only_rectal_campylobacter_genus)
#divide 
(12/2330)*100



#try pruning for cpp-only when its campylobacter filtered 
CPP_only_rectal_campylobacter_genus <- subset_taxa(CPP_only_rectal_campylobacter, Genus == "g__Campylobacter")
CPP_only_rectal_campylobacter_genus <- prune_taxa(taxa_sums(CPP_only_rectal_campylobacter_genus) > 0, CPP_only_rectal_campylobacter_genus)
ntaxa(CPP_only_rectal_campylobacter_genus)

#try pruning for cpp-endo when its campylobacter filtered 
CPP_endo_rectal_campylobacter_genus <- subset_taxa(CPP_endo_rectal_campylobacter, Genus == "g__Campylobacter")
CPP_endo_rectal_campylobacter_genus <- prune_taxa(taxa_sums(CPP_endo_rectal_campylobacter_genus) > 0, CPP_endo_rectal_campylobacter_genus)
ntaxa(CPP_endo_rectal_campylobacter_genus)

#try pruning for cpp-healthy when its campylobacter filtered 
CPP_healthy_rectal_campylobacter_genus <- subset_taxa(CPP_healthy_rectal_campylobacter, Genus == "g__Campylobacter")
CPP_healthy_rectal_campylobacter_genus <- prune_taxa(taxa_sums(CPP_healthy_rectal_campylobacter_genus) > 0, CPP_healthy_rectal_campylobacter_genus)
ntaxa(CPP_healthy_rectal_campylobacter_genus)

#try pruning entire taxa table for cpp_only
CPP_only_rectal_campylobacter_prune <- prune_taxa(taxa_sums(CPP_only_rectal_campylobacter) > 0, CPP_only_rectal_campylobacter)
ntaxa(CPP_only_rectal_campylobacter_prune)

#try pruning entire taxa table for cpp_endo 
CPP_endo_rectal_campylobacter_prune <- prune_taxa(taxa_sums(CPP_endo_rectal_campylobacter) > 0, CPP_endo_rectal_campylobacter)
ntaxa(CPP_endo_rectal_campylobacter_prune)

#try pruning entire taxa table for cpp_healthy
CPP_healthy_rectal_campylobacter_prune <- prune_taxa(taxa_sums(CPP_healthy_rectal_campylobacter) > 0, CPP_healthy_rectal_campylobacter)
ntaxa(CPP_healthy_rectal_campylobacter_prune)


#try just getting abundance 
taxasum <- taxa_sums(CPP_only_rectal_campylobacter)["g__Campylobacter"]


#trying again using another approach 

relative_abundance <- transform_sample_counts(phyloseq, function(x) x / sum(x))
ps_genus <- subset_taxa(relative_abundance, Genus == "g__Campylobacter")
ps_genus <- prune_taxa(taxa_sums(ps_genus) > 0, ps_genus)
df <- psmelt(ps_genus)

ggplot(df, aes(x = Host_disease, y = Abundance)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  theme_minimal()


#try again with only ASv
relative_abundance <- transform_sample_counts(phyloseq, function(x) x / sum(x))
ps_asv <- prune_taxa(taxa_names(relative_abundance) == "8365c56e7d1f57079c8821b426206d94", relative_abundance)
ps_genus <- prune_taxa(taxa_sums(ps_asv) > 0, ps_asv)
df <- psmelt(ps_genus)

ggplot(df, aes(x = Host_disease, y = Abundance)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  stat_compare_means(comparisons = comparisons,
                     method = "wilcox.test")
  theme_minimal()

#do stats
kruskal.test(Abundance ~ Host_disease, data = df)
comparisons <- pairwise.wilcox.test(df$Abundance, df$Host_disease, p.adjust.method = "BH")
