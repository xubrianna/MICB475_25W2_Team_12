library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)

####LOAD AND MODIFY PHYLOSEQ####
#load filtered phyloseq without rarification (expects uneven sequencing depths)
phyloseq <- readRDS('data/phyloseq_filtered.rds')

# Transform OTU table into relative abundance from absolute counts
CPP_RA <- transform_sample_counts(phyloseq, fun = function(x) x / sum(x))

#### STRATIFY BY HOST DISEASE ####
CPP_only <- subset_samples(CPP_RA, `Host_disease` == "CPP")
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
  Healthy = CPP_healthy_ASVs_rectal_0.001,
  "CPP Only" = CPP_only_ASVs_rectal_0.001,
  "CPP Endo" = CPP_endo_ASVs_rectal_0.001
)

venn_all_diseases_rectal_0.001 <- ggVennDiagram(
  x = CPP_list_rectal_0.001,
  category.names = c("Healthy", "CPP Only", "CPP Endo")
) + scale_fill_gradient(low = "#1B98E026", high = "#225ea8")
venn_all_diseases_rectal_0.001

ggsave("results/aim2/04-core_microbiome/venn_all_diseases_rectal_0.001.png",
       venn_all_diseases_rectal_0.001, width = 7, height = 7)

#get the taxonomic information for CPP_endo disease group and find the ASVs that are unique to the core of CPP-endo 
tax_mat_rectal <- tax_table(CPP_endo_rectal)

core_taxonomy_rectal <- as.data.frame(tax_mat_rectal[CPP_endo_ASVs_rectal_0.001, ])
print(head(core_taxonomy_rectal)

unique_to_CPP_endo_rectal<- setdiff(CPP_endo_ASVs_rectal_0.001, union(CPP_only_ASVs_rectal_0.001, CPP_healthy_ASVs_rectal_0.001))
print(unique_to_CPP_endo_rectal)

#prune 
prune_taxa(CPP_only_ASVs_rectal_0.001,CPP_RA) %>%
  tax_table()

prune_taxa(CPP_endo_ASVs_rectal_0.001,CPP_RA) %>%
  tax_table()

prune_taxa(CPP_healthy_ASVs_rectal_0.001,CPP_RA) %>%
  tax_table()

#plot at the genus level - gives the samples within each group that contain the ASVs only present within CPP only and abundance 
prune_taxa_rectal_CPP_only_0.001 <- prune_taxa(CPP_only_ASVs_rectal_0.001,CPP_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`Host_disease`, scales ="free")

ggsave("results/aim2/04-core_microbiome/prune_taxa_rectal_CPP_only_0.001.png", prune_taxa_rectal_CPP_only_0.001, width = 7, height = 7)

#plot at the genus level - gives the samples within each group that contain ASVs only present within CPP Endo and abundance 
prune_taxa_rectal_CPP_endo_0.001 <- prune_taxa(CPP_endo_ASVs_rectal_0.001,CPP_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`Host_disease`, scales ="free")

ggsave("results/aim2/04-core_microbiome/prune_taxa_rectal_CPP_endo_0.001.png", prune_taxa_rectal_CPP_endo_0.001, width = 7, height = 7)

##plot at the genus level - gives the samples within each group that contain ASVs only present within CPP healthy and abundance 
prune_taxa_rectal_CPP_healthy_0.001 <- prune_taxa(CPP_healthy_ASVs_rectal_0.001,CPP_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`Host_disease`, scales ="free")

ggsave("results/aim2/04-core_microbiome/prune_taxa_rectal_CPP_healthy_0.001.png", prune_taxa_rectal_CPP_healthy_0.001, width = 7, height = 7)

#try detection of 0 (absence/presence) with the 0.5 prevalence  
CPP_only_ASVs_rectal_0 <- core_members(CPP_only_rectal, detection=0, prevalence = 0.5)
CPP_endo_ASVs_rectal_0 <- core_members(CPP_endo_rectal, detection=0, prevalence = 0.5)
CPP_healthy_ASVs_rectal_0 <- core_members(CPP_healthy_rectal, detection=0, prevalence = 0.5)

CPP_list_rectal_0 <- list(Healthy = CPP_healthy_ASVs_rectal_0, CPP_Only = CPP_only_ASVs_rectal_0, CPP_Endo = CPP_endo_ASVs_rectal_0)

venn_all_diseases_rectal_0<- ggVennDiagram(x = CPP_list_rectal_0)+ scale_fill_gradient(low ="#c7e9b4", high = "#225ea8")
venn_all_diseases_rectal_0

ggsave("results/aim2/04-core_microbiome/venn_all_diseases_rectal_0.png", venn_all_diseases_rectal_0, width = 7, height = 7)

#try detection of 0.01 with 0.5 prevalence to keep only abundant taxa 
CPP_only_ASVs_rectal_0.01 <- core_members(CPP_only_rectal, detection=0.01, prevalence = 0.5)
CPP_endo_ASVs_rectal_0.01 <- core_members(CPP_endo_rectal, detection=0.01, prevalence = 0.5)
CPP_healthy_ASVs_rectal_0.01 <- core_members(CPP_healthy_rectal, detection=0.01, prevalence = 0.5)

CPP_list_rectal_0.01 <- list(Healthy = CPP_healthy_ASVs_rectal_0.01, CPP_Only = CPP_only_ASVs_rectal_0.01, CPP_Endo = CPP_endo_ASVs_rectal_0.01)

venn_all_diseases_rectal_0.01<- ggVennDiagram(x = CPP_list_rectal_0.01)+ scale_fill_gradient(low ="#c7e9b4", high = "#225ea8")
venn_all_diseases_rectal_0.01

ggsave("results/aim2/04-core_microbiome/venn_all_diseases_rectal_0.01.png", venn_all_diseases_rectal_0.01, width = 7, height = 7)

#try doing an abundance of 0.3 with 0.001 prevalence threshold 
CPP_only_ASVs_rectal_0.3 <- core_members(CPP_only_rectal, detection=0.001, prevalence = 0.3)
CPP_endo_ASVs_rectal_0.3 <- core_members(CPP_endo_rectal, detection=0.001, prevalence = 0.3)
CPP_healthy_ASVs_rectal_0.3 <- core_members(CPP_healthy_rectal, detection=0.001, prevalence = 0.3)

CPP_list_rectal_0.3 <- list(Healthy = CPP_healthy_ASVs_rectal_0.3, CPP_Only = CPP_only_ASVs_rectal_0.3, CPP_Endo = CPP_endo_ASVs_rectal_0.3)

venn_all_diseases_rectal_0.3 <- ggVennDiagram(x = CPP_list_rectal_0.3)+ scale_fill_gradient(low ="#c7e9b4", high = "#225ea8")
venn_all_diseases_rectal_0.3

ggsave("results/aim2/04-core_microbiome/venn_all_diseases_rectal_0.3.png", venn_all_diseases_rectal_0.3, width = 7, height = 7)

#### VAGINAL COREMICROBIOME #### 

#stratifying by the vaginal sample instead now  
CPP_only_vaginal <- subset_samples(CPP_only,`collection_method` == "vaginal_swab")
CPP_endo_vaginal <- subset_samples(CPP_endo,`collection_method` == "vaginal_swab")
CPP_healthy_vaginal <-subset_samples(CPP_healthy, `collection_method` == "vaginal_swab")

CPP_only_ASVs_vaginal_0.001 <- core_members(CPP_only_vaginal, detection=0.001, prevalence = 0.5)
CPP_endo_ASVs_vaginal_0.001 <- core_members(CPP_endo_vaginal, detection=0.001, prevalence = 0.5)
CPP_healthy_ASVs_vaginal_0.001 <- core_members(CPP_healthy_vaginal, detection=0.001, prevalence = 0.5)

CPP_list_vaginal_0.001 <- list(Healthy = CPP_healthy_ASVs_vaginal_0.001, CPP_Only = CPP_only_ASVs_vaginal_0.001, CPP_Endo = CPP_endo_ASVs_vaginal_0.001)

venn_all_diseases_vaginal_0.001 <- ggVennDiagram(x = CPP_list_vaginal_0.001, category.names = c("Healthy", "CPP Only", "CPP Endo"))+ scale_fill_gradient(low = "#1B98E026", high="#225ea8")
venn_all_diseases_vaginal_0.001

ggsave("results/aim2/04-core_microbiome/venn_all_diseases_vaginal_0.001.png", venn_all_diseases_vaginal_0.001, width = 7, height = 7)
 
#prune 
prune_taxa(CPP_only_ASVs_vaginal_0.001,CPP_RA) %>%
  tax_table()

prune_taxa(CPP_endo_ASVs_vaginal_0.001,CPP_RA) %>%
  tax_table()

prune_taxa(CPP_healthy_ASVs_vaginal_0.001,CPP_RA) %>%
  tax_table()

#plot at the genus level - gives the samples within each group that contain the ASVs only present within CPP only and abundance 
prune_taxa_vaginal_CPP_only_0.001 <- prune_taxa(CPP_only_ASVs_vaginal_0.001,CPP_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`Host_disease`, scales ="free")

ggsave("results/aim2/04-core_microbiome/prune_taxa_vaginal_CPP_only_0.001.png", prune_taxa_vaginal_CPP_only_0.001, width = 7, height = 7)

#plot at the genus level - gives the samples within each group that contain ASVs only present within CPP Endo and abundance 
prune_taxa_vaginal_CPP_endo_0.001 <- prune_taxa(CPP_endo_ASVs_vaginal_0.001,CPP_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`Host_disease`, scales ="free")

ggsave("results/aim2/04-core_microbiome/prune_taxa_vaginal_CPP_endo_0.001.png", prune_taxa_vaginal_CPP_endo_0.001, width = 7, height = 7)

##plot at the genus level - gives the samples within each group that contain ASVs only present within CPP healthy and abundance 
prune_taxa_vaginal_CPP_healthy_0.001 <- prune_taxa(CPP_healthy_ASVs_vaginal_0.001,CPP_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`Host_disease`, scales ="free")

ggsave("results/aim2/04-core_microbiome/prune_taxa_vaginal_CPP_healthy_0.001.png", prune_taxa_vaginal_CPP_healthy_0.001, width = 7, height = 7)

#try detection of 0 with 0.5 prevalence
CPP_only_ASVs_vaginal_0 <- core_members(CPP_only_vaginal, detection=0, prevalence = 0.5)
CPP_endo_ASVs_vaginal_0 <- core_members(CPP_endo_vaginal, detection=0, prevalence = 0.5)
CPP_healthy_ASVs_vaginal_0 <- core_members(CPP_healthy_vaginal, detection=0, prevalence = 0.5)

CPP_list_vaginal_0 <- list(Healthy = CPP_healthy_ASVs_vaginal_0, CPP_Only = CPP_only_ASVs_vaginal_0, CPP_Endo = CPP_endo_ASVs_vaginal_0)

venn_all_diseases_vaginal_0 <- ggVennDiagram(x = CPP_list_vaginal_0)+ scale_fill_gradient(low ="#c7e9b4", high = "#225ea8")
venn_all_diseases_vaginal_0

ggsave("results/aim2/04-core_microbiome/venn_all_diseases_vaginal_0.png", venn_all_diseases_vaginal_0, width = 7, height = 7)

#get the taxonomic information for CPP_endo disease group and find the ASVs that are unique to the core of CPP-endo 
tax_mat_vaginal <- tax_table(CPP_endo_vaginal)

core_taxonomy_vaginal <- as.data.frame(tax_mat[CPP_endo_ASVs_vaginal_0 , ])
print(head(core_taxonomy))

unique_to_CPP_endo_vaginal <- setdiff(CPP_endo_ASVs_vaginal_0, union(CPP_only_ASVs_vaginal_0, CPP_healthy_ASVs_vaginal_0))
print(unique_to_CPP_endo_vaginal)

#try detection of 0.01 with 0.5 prevalence
CPP_only_ASVs_vaginal_0.01 <- core_members(CPP_only_vaginal, detection=0.01, prevalence = 0.5)
CPP_endo_ASVs_vaginal_0.01 <- core_members(CPP_endo_vaginal, detection=0.01, prevalence = 0.5)
CPP_healthy_ASVs_vaginal_0.01 <- core_members(CPP_healthy_vaginal, detection=0.01, prevalence = 0.5)

CPP_list_vaginal_0.01 <- list(Healthy = CPP_healthy_ASVs_vaginal_0.01, CPP_Only = CPP_only_ASVs_vaginal_0.01, CPP_Endo = CPP_endo_ASVs_vaginal_0.01)

venn_all_diseases_vaginal_0.01 <- ggVennDiagram(x = CPP_list_vaginal_0.01)+ scale_fill_gradient(low ="#c7e9b4", high = "#225ea8")
venn_all_diseases_vaginal_0.01

ggsave("results/aim2/04-core_microbiome/venn_all_diseases_vaginal_0.01.png", venn_all_diseases_vaginal_0.01, width = 7, height = 7)

#try doing an abundnace of 0.3 with detection of 0.001
CPP_only_ASVs_vaginal_0.3 <- core_members(CPP_only_vaginal, detection=0.001, prevalence = 0.3)
CPP_endo_ASVs_vaginal_0.3 <- core_members(CPP_endo_vaginal, detection=0.001, prevalence = 0.3)
CPP_healthy_ASVs_vaginal_0.3 <- core_members(CPP_healthy_vaginal, detection=0.001, prevalence = 0.3)

print(CPP_only_ASVs_vaginal_0.3)
print (CPP_endo_ASVs_vaginal_0.3)
print (CPP_healthy_ASVs_vaginal_0.3)

CPP_list_vaginal_0.3 <- list(Healthy = CPP_healthy_ASVs_vaginal_0.3, CPP_Only = CPP_only_ASVs_vaginal_0.3, CPP_Endo = CPP_endo_ASVs_vaginal_0.3)

venn_all_diseases_vaginal_0.3 <- ggVennDiagram(x = CPP_list_vaginal_0.3)+ scale_fill_gradient(low ="#c7e9b4", high = "#225ea8")
venn_all_diseases_vaginal_0.3

ggsave("results/aim2/04-core_microbiome/venn_all_diseases_vaginal_0.3.png", venn_all_diseases_vaginal_0.3, width = 7, height = 7)




#testing a prevalence of 0.2 and a detection of 0.001

CPP_only_ASVs_vaginal_prevalence_0.2 <- core_members(CPP_only_vaginal, detection=0.001, prevalence = 0.2)
CPP_endo_ASVs_vaginal_prevalence_0.2 <- core_members(CPP_endo_vaginal, detection=0.001, prevalence = 0.2)
CPP_healthy_ASVs_vaginal_prevalence_0.2 <- core_members(CPP_healthy_vaginal, detection=0.001, prevalence = 0.2)

CPP_list_vaginal_prevalence_0.2 <- list(Healthy = CPP_healthy_ASVs_vaginal_prevalence_0.2, CPP_Only = CPP_only_ASVs_vaginal_prevalence_0.2, CPP_Endo = CPP_endo_ASVs_vaginal_prevalence_0.2)

venn_all_diseases_vaginal_prevalence_0.2 <- ggVennDiagram(x = CPP_list_vaginal_prevalence_0.2)+ scale_fill_gradient(low ="#c7e9b4", high = "#225ea8")
venn_all_diseases_vaginal_prevalence_0.2
ggsave("results/aim2/04-core_microbiome/venn_all_diseases_vaginal_prevalence_0.2.png", venn_all_diseases_vaginal_prevalence_0.2, width = 7, height = 7)



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
  )

# Mean relative abundance per disease group
campy_summary <- campy_df %>%
  group_by(Host_disease) %>%
  summarise(Mean_RA = mean(Campy_RA), .groups = "drop")
print(campy_summary)

# Plot
plot_campylobacter <- ggplot(campy_df, aes(x = Host_disease, y = Campy_RA, fill = Host_disease)) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2.5, color = "black") +
  scale_fill_manual(values = c("Control" = "#c7e9b4", "CPP" = "#41b6c4", "CPP Endo" = "#225ea8")) +
  theme_classic() +
  xlab("Host Disease") +
  ylab("Campylobacter Relative Abundance") +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 18),
    plot.margin = margin(15, 15, 15, 15)
  )

ggsave("results/aim2/04-core_microbiome/plot_campylobacter.png",
       plot_campylobacter, width = 7, height = 7)

### March 28 trying to find what is unqiue to healthy in rectal coremicrobiome
#get the taxonomic information for CPP_endo disease group and find the ASVs that are unique to the core of CPP-endo 
tax_mat_rectal_2 <- tax_table(CPP_endo_rectal)

core_taxonomy_rectal_2 <- as.data.frame(tax_mat_rectal_2[CPP_healthy_ASVs_rectal_0.001, ])
print(head(core_taxonomy_rectal_2)
      
unique_to_CPP_endo_rectal_2<- setdiff(CPP_healthy_ASVs_rectal_0.001, union(CPP_only_ASVs_rectal_0.001, CPP_endo_ASVs_rectal_0.001))
print(unique_to_CPP_endo_rectal_2)
