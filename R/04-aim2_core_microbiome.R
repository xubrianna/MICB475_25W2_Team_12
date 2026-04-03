
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)

install.packages("ggpubr")   # only once
library(ggpubr)

####LOAD AND MODIFY PHYLOSEQ####
#load filtered phyloseq without rarification (expects uneven sequencing depths)
phyloseq <- readRDS('data/phyloseq_filtered.rds')

#transform OTU table into relative abundance from absolute counts
CPP_RA <- transform_sample_counts(phyloseq, fun=function(x) x/sum(x))

#view 
otu_table(phyloseq)
sample_data(phyloseq)

#### STRATIFY BY HOST DISEASE #### 

#try stratifying by the rectal sample after stratify by disease 
CPP_only <- subset_samples(CPP_RA, `Host_disease`=="CPP")
CPP_endo <- subset_samples(CPP_RA, `Host_disease`=="CPP Endo")
CPP_healthy <-subset_samples(CPP_RA, `Host_disease`=="Control")

#### RECTAL COREMICROBIOME #### 
#filter by rectal samples 
CPP_only_rectal <- subset_samples(CPP_only,`collection_method` == "rectal_swab")
CPP_endo_rectal <- subset_samples(CPP_endo,`collection_method` == "rectal_swab")
CPP_healthy_rectal <-subset_samples(CPP_healthy, `collection_method` == "rectal_swab")

#assign ASVs using 0.001 detection (filters out rare taxa) and 0.5 prevalence
CPP_only_ASVs_rectal_0.001 <- core_members(CPP_only_rectal, detection=0.001, prevalence = 0.5)
CPP_endo_ASVs_rectal_0.001 <- core_members(CPP_endo_rectal, detection=0.001, prevalence = 0.5)
CPP_healthy_ASVs_rectal_0.001 <- core_members(CPP_healthy_rectal, detection=0.001, prevalence = 0.5)

CPP_list_rectal_0.001 <- list(Healthy = CPP_healthy_ASVs_rectal_0.001, CPP_Only = CPP_only_ASVs_rectal_0.001, CPP_Endo = CPP_endo_ASVs_rectal_0.001)

venn_all_diseases_rectal_0.001 <- ggVennDiagram(x = CPP_list_rectal_0.001, category.names = c("Healthy", "CPP Only", "CPP Endo"))+ scale_fill_gradient(low = "#1B98E026", high="#225ea8") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(plot.title = element_text(size=20, face="bold"))
venn_all_diseases_rectal_0.001

ggsave("results/aim2/04-core_microbiome/venn_all_diseases_rectal_0.001.png", venn_all_diseases_rectal_0.001, width = 7, height = 7)

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

venn_all_diseases_rectal_0<- ggVennDiagram(x = CPP_list_rectal_0)+ scale_fill_gradient(low ="#c7e9b4", high = "#225ea8") +
  labs(title = "Rectal Samples") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(plot.title = element_text(size=20, face="bold"))
venn_all_diseases_rectal_0

ggsave("results/aim2/04-core_microbiome/venn_all_diseases_rectal_0.png", venn_all_diseases_rectal_0, width = 7, height = 7)

#try detection of 0.01 with 0.5 prevalence to keep only abundant taxa 
CPP_only_ASVs_rectal_0.01 <- core_members(CPP_only_rectal, detection=0.01, prevalence = 0.5)
CPP_endo_ASVs_rectal_0.01 <- core_members(CPP_endo_rectal, detection=0.01, prevalence = 0.5)
CPP_healthy_ASVs_rectal_0.01 <- core_members(CPP_healthy_rectal, detection=0.01, prevalence = 0.5)

CPP_list_rectal_0.01 <- list(Healthy = CPP_healthy_ASVs_rectal_0.01, CPP_Only = CPP_only_ASVs_rectal_0.01, CPP_Endo = CPP_endo_ASVs_rectal_0.01)

venn_all_diseases_rectal_0.01<- ggVennDiagram(x = CPP_list_rectal_0.01)+ scale_fill_gradient(low ="#c7e9b4", high = "#225ea8") +
  labs(title = "Rectal Samples") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(plot.title = element_text(size=20, face="bold"))
venn_all_diseases_rectal_0.01

ggsave("results/aim2/04-core_microbiome/venn_all_diseases_rectal_0.01.png", venn_all_diseases_rectal_0.01, width = 7, height = 7)

#try doing an abundance of 0.3 with 0.001 prevalence threshold 
CPP_only_ASVs_rectal_0.3 <- core_members(CPP_only_rectal, detection=0.001, prevalence = 0.3)
CPP_endo_ASVs_rectal_0.3 <- core_members(CPP_endo_rectal, detection=0.001, prevalence = 0.3)
CPP_healthy_ASVs_rectal_0.3 <- core_members(CPP_healthy_rectal, detection=0.001, prevalence = 0.3)

CPP_list_rectal_0.3 <- list(Healthy = CPP_healthy_ASVs_rectal_0.3, CPP_Only = CPP_only_ASVs_rectal_0.3, CPP_Endo = CPP_endo_ASVs_rectal_0.3)

venn_all_diseases_rectal_0.3 <- ggVennDiagram(x = CPP_list_rectal_0.3)+ scale_fill_gradient(low ="#c7e9b4", high = "#225ea8")+
  labs(title = "Rectal Samples") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(plot.title = element_text(size=20, face="bold"))
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

venn_all_diseases_vaginal_0.001 <- ggVennDiagram(x = CPP_list_vaginal_0.001, category.names = c("Healthy", "CPP Only", "CPP Endo"))+ scale_fill_gradient(low = "#1B98E026", high="#225ea8")+
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(plot.title = element_text(size=20, face="bold"))
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

venn_all_diseases_vaginal_0 <- ggVennDiagram(x = CPP_list_vaginal_0)+ scale_fill_gradient(low ="#c7e9b4", high = "#225ea8")+
  labs(title = "Vaginal Samples") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(plot.title = element_text(size=20, face="bold"))
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

venn_all_diseases_vaginal_0.01 <- ggVennDiagram(x = CPP_list_vaginal_0.01)+ scale_fill_gradient(low ="#c7e9b4", high = "#225ea8")+
  labs(title = "Vaginal Samples") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(plot.title = element_text(size=20, face="bold"))
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

venn_all_diseases_vaginal_0.3 <- ggVennDiagram(x = CPP_list_vaginal_0.3)+ scale_fill_gradient(low ="#c7e9b4", high = "#225ea8")+
  labs(title = "Vaginal Samples") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(plot.title = element_text(size=20, face="bold"))
venn_all_diseases_vaginal_0.3

ggsave("results/aim2/04-core_microbiome/venn_all_diseases_vaginal_0.3.png", venn_all_diseases_vaginal_0.3, width = 7, height = 7)




#testing a prevalence of 0.2 and a detection of 0.001

CPP_only_ASVs_vaginal_prevalence_0.2 <- core_members(CPP_only_vaginal, detection=0.001, prevalence = 0.2)
CPP_endo_ASVs_vaginal_prevalence_0.2 <- core_members(CPP_endo_vaginal, detection=0.001, prevalence = 0.2)
CPP_healthy_ASVs_vaginal_prevalence_0.2 <- core_members(CPP_healthy_vaginal, detection=0.001, prevalence = 0.2)

CPP_list_vaginal_prevalence_0.2 <- list(Healthy = CPP_healthy_ASVs_vaginal_prevalence_0.2, CPP_Only = CPP_only_ASVs_vaginal_prevalence_0.2, CPP_Endo = CPP_endo_ASVs_vaginal_prevalence_0.2)

venn_all_diseases_vaginal_prevalence_0.2 <- ggVennDiagram(x = CPP_list_vaginal_prevalence_0.2)+ scale_fill_gradient(low ="#c7e9b4", high = "#225ea8")+
  labs(title = "Vaginal Samples") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(plot.title = element_text(size=20, face="bold"))
venn_all_diseases_vaginal_prevalence_0.2
ggsave("results/aim2/04-core_microbiome/venn_all_diseases_vaginal_prevalence_0.2.png", venn_all_diseases_vaginal_prevalence_0.2, width = 7, height = 7)



#### March 24 finding percentage of samples with campylobacter total #### 
otu_table(phyloseq)
ASV_table <- as.data.frame(as.matrix(otu_table(phyloseq)))
ASV_table_flipped <- as.data.frame(t(ASV_table))
campylobacter <- "8365c56e7d1f57079c8821b426206d94"
samples_with_campylobacter <- sum(ASV_table_flipped[, campylobacter] > 0)
total_samples <- nrow(ASV_table_flipped)
percentage <- (samples_with_campylobacter / total_samples) * 100
print(percentage)

#now trying to stratify using subsetted phyloseqs 

CPP_only_campylobacter <- subset_samples(phyloseq, `Host_disease`=="CPP")
CPP_endo_campylobacter <- subset_samples(phyloseq, `Host_disease`=="CPP Endo")
CPP_healthy_campylobacter <-subset_samples(phyloseq, `Host_disease`=="Control")

#filter by rectal samples 
CPP_only_rectal_campylobacter <- subset_samples(CPP_only_campylobacter,`collection_method` == "rectal_swab")
CPP_endo_rectal_campylobacter <- subset_samples(CPP_endo_campylobacter,`collection_method` == "rectal_swab")
CPP_healthy_rectal_campylobacter <-subset_samples(CPP_healthy_campylobacter, `collection_method` == "rectal_swab")

otu_table(CPP_only_rectal_campylobacter)
ASV_table_campylobacter <- as.data.frame(as.matrix(otu_table(CPP_only_rectal_campylobacter)))
ASV_table_flipped_campylobacter <- as.data.frame(t(ASV_table_campylobacter))
campylobacter <- "8365c56e7d1f57079c8821b426206d94"
samples_with_campylobacter_2 <- sum(ASV_table_flipped_campylobacter[, campylobacter] > 0)
total_samples_2 <- nrow(ASV_table_flipped_campylobacter)
percentage <- (samples_with_campylobacter_2 / total_samples_2) * 100
print(percentage)

otu_table(CPP_endo_rectal_campylobacter)
ASV_table_CPP_endo_rectal_campylobacter <- as.data.frame(as.matrix(otu_table(CPP_endo_rectal_campylobacter)))
ASV_table_flipped_CPP_endo_rectal_campylobacter <- as.data.frame(t(ASV_table_CPP_endo_rectal_campylobacter))
campylobacter <- "8365c56e7d1f57079c8821b426206d94"
samples_with_campylobacter_2 <- sum(ASV_table_flipped_CPP_endo_rectal_campylobacter[, campylobacter] > 0)
total_samples_2 <- nrow(ASV_table_flipped_CPP_endo_rectal_campylobacter)
percentage <- (samples_with_campylobacter_2 / total_samples_2) * 100
print(percentage)


otu_table(CPP_healthy_rectal_campylobacter)
ASV_table_CPP_healthy_rectal_campylobacter <- as.data.frame(as.matrix(otu_table(CPP_healthy_rectal_campylobacter)))
ASV_table_flipped_CPP_healthy_rectal_campylobacter  <- as.data.frame(t(ASV_table_CPP_healthy_rectal_campylobacter))
campylobacter <- "8365c56e7d1f57079c8821b426206d94"
samples_with_campylobacter_2 <- sum(ASV_table_flipped_CPP_healthy_rectal_campylobacter [, campylobacter] > 0)
total_samples_2 <- nrow(ASV_table_CPP_healthy_rectal_campylobacter )
percentage <- (samples_with_campylobacter_2 / total_samples_2) * 100
print(percentage)

#make vector 
percentages_camp <- c(50, 76.31, 0.429)
host_disease <- c("CPP Only", "CPP Endo", "Healthy")
dataframe <- data.frame(Percent_Campylobacter = percentages_camp, Disease_group = host_disease)

library(ggplot2)
plot_campylobacter <- ggplot(data = dataframe, aes(x = Disease_group, y = Percent_Campylobacter, fill = Disease_group))+ geom_bar(stat = "identity") + geom_col() + scale_fill_manual(values = c("CPP Endo" = "#225ea8", "CPP Only" = "#41b6c4", "Healthy" = "#c7e9b4")) +theme_classic() +
  xlab("Host Disease") +
  ylab ("Percentage of Samples with Campylobacter")+ 
  theme(legend.position = "none")
ggsave("results/aim2/04-core_microbiome/plot_campylobacter.png", plot_campylobacter , width = 7, height = 7)

### March 28 trying to find what is unqiue to healthy in rectal coremicrobiome
#get the taxonomic information for CPP_endo disease group and find the ASVs that are unique to the core of CPP-endo 
tax_mat_rectal_2 <- tax_table(CPP_endo_rectal)

core_taxonomy_rectal_2 <- as.data.frame(tax_mat_rectal_2[CPP_healthy_ASVs_rectal_0.001, ])
print(head(core_taxonomy_rectal_2)
      
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
