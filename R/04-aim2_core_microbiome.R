
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)

#load filtered phyloseq without rarification (expects uneven sequencing depths)
phyloseq <- readRDS('data/phyloseq_filtered.rds')

#view 
otu_table(phyloseq)
sample_data(phyloseq)

#transform OTU table into relative abundance from absolute counts 
CPP_RA <- transform_sample_counts(phyloseq, fun=function(x) x/sum(x))

#RUN WITHOUT STRATIFY BY BODY SITE 
#subset the dataset by our variables within the predictor of disease 
CPP_only <- subset_samples(phyloseq, `Host_disease`=="CPP")
CPP_endo <- subset_samples(phyloseq, `Host_disease`=="CPP Endo")
CPP_healthy <-subset_samples(phyloseq, `Host_disease`=="Control")

#set thresholds and prevalence for each subset, gives which ASVs meet the threshold 
CPP_only_ASVs <- core_members(CPP_only, detection=0, prevalence = 0.5)
CPP_endo_ASVs <- core_members(CPP_endo, detection=0, prevalence = 0.5)
CPP_healthy_ASVs <- core_members(CPP_healthy, detection=0, prevalence = 0.5)

#prune 
prune_taxa(CPP_only_ASVs,phyloseq) %>%
  tax_table()

prune_taxa(CPP_endo_ASVs,phyloseq) %>%
  tax_table()

prune_taxa(CPP_healthy_ASVs,phyloseq) %>%
  tax_table()

#plot at the genus level - gives the samples within each group that contain the ASVs only present within CPP only and abundance 
prune_taxa(CPP_only_ASVs,phyloseq) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`Host_disease`, scales ="free")

#plot at the genus level - gives the samples within each group that contain ASVs only present within CPP Endo and abundance 
prune_taxa(CPP_endo_ASVs,phyloseq) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`Host_disease`, scales ="free")

##plot at the genus level - gives the samples within each group that contain ASVs only present within CPP healthy and abundance 
prune_taxa(CPP_healthy_ASVs,phyloseq) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`Host_disease`, scales ="free")

#combine the CPP healthy and CPP only ASV vectors together to create a larger list 
CPP_list <- list(Healthy = CPP_healthy_ASVs, CPP_Only = CPP_only_ASVs)

# Create a Venn diagram to compare healthy and CPP 
venn_CPP_only_vs_Healthy <- ggVennDiagram(x = CPP_list)
venn_CPP_only_vs_Healthy
#there are 6 ASVs that are unique across all the samples in Healthy group 
#there are 8 ASVs that are shared across all the samples shared for all Disease states 
#there is only 1 ASVs that is shared across all the samples shared for the CPP group

ggsave("venn_CPP_only_vs_healthy.png", first_venn)

#Create a list of CPP and CPP endo ASVs and Venn diagram to compare CPP and CPP_endo
CPP_list_2 <- list(CPP_Only = CPP_only_ASVs, CPP_Endo = CPP_endo_ASVs)

venn_CPP_only_vs_CPP_endo <- ggVennDiagram(x = CPP_list_2)
venn_CPP_only_vs_CPP_endo
#there are 0 ASVs that are unique across all the samples in CPP Endo group 
#there are 7 ASVs that are shared across all the samples shared for all Disease states 
#there is only 2 ASVs that is shared across all the samples shared for the CPP Only group

#Create a list of CPP healthy and CPP endo ASVs and Venn diagram to compare CPP and CPP_endo
CPP_list_3 <- list(CPP_Healthy = CPP_healthy_ASVs, CPP_Endo = CPP_endo_ASVs)

venn_CPP_healthy_vs_CPP_endo <- ggVennDiagram(x = CPP_list_3)
venn_CPP_healthy_vs_CPP_endo
#there are 0 ASVs that are unique across all the samples in CPP Endo group 
#there are 7 ASVs that are shared across all the samples shared for all Disease states 
#there is only 7 ASVs that is shared across all the samples shared for the CPP Healthy group

#create a list of all three conditions 
CPP_list_full <- list(CPP_Healthy= CPP_healthy_ASVs, CPP_Endo = CPP_endo_ASVs, CPP_Only = CPP_only_ASVs)

venn_all_diseases <- ggVennDiagram(x = CPP_list_full)
venn_all_diseases
#there are 6 ASVs that are common to the healthy group 
#there are no ASVs that are shared to the CPP Endo group 
#there are no ASVs that are shared to the CPP Only group 

#try with a different threshold AND prevalence for each subset that is less stringent 
CPP_only_ASVs_0.4 <- core_members(CPP_only, detection=0.01, prevalence = 0.4)
CPP_endo_ASVs_0.4 <- core_members(CPP_endo, detection=0.01, prevalence = 0.4)
CPP_healthy_ASVs_0.4 <- core_members(CPP_healthy, detection=0.01, prevalence = 0.4)

CPP_list_0.4 <- list(Healthy = CPP_healthy_ASVs_0.4, CPP_Only = CPP_only_ASVs_0.4, CPP_Endo = CPP_endo_ASVs_0.4)

# Create a Venn diagram to compare healthy and CPP 
venn_all_diseases_0.4 <- ggVennDiagram(x = CPP_list_2)
venn_all_diseases_0.4

##STRATIFY BY BODY SITE 
#try stratifying by the rectal sample after stratify by disease 
CPP_only_2 <- subset_samples(phyloseq, `Host_disease`=="CPP")
CPP_endo_2 <- subset_samples(phyloseq, `Host_disease`=="CPP Endo")
CPP_healthy_2 <-subset_samples(phyloseq, `Host_disease`=="Control")

CPP_only_3 <- subset_samples(CPP_only_2,`collection_method` == "rectal_swab")
CPP_endo_3 <- subset_samples(CPP_endo_2,`collection_method` == "rectal_swab")
CPP_healthy_3 <-subset_samples(CPP_healthy_2, `collection_method` == "rectal_swab")

#then assign ASVs
CPP_only_ASVs_3 <- core_members(CPP_only_3, detection=0, prevalence = 0.5)
CPP_endo_ASVs_3 <- core_members(CPP_endo_3, detection=0, prevalence = 0.5)
CPP_healthy_ASVs_3 <- core_members(CPP_healthy_3, detection=0, prevalence = 0.5)

#get the names of the genus 
CPP_endo_ASVs_3_Names <- tax_table(phyloseq)[CPP_endo_ASVs_3, "Genus"]

#try finding what genus are unique to a particular group 
Unique_to_CPP_Endo <- setdiff(CPP_endo_ASVs_3,CPP_only_ASVs_3)
Unique_to_CPP_Endo_Genera <- tax_table(phyloseq)[Unique_to_CPP_Endo, "Genus"]

#then make a list, either with all of the diseases or just two 
CPP_list_3 <- list(Healthy = CPP_healthy_ASVs_3, CPP_Only = CPP_only_ASVs_3, CPP_Endo = CPP_endo_ASVs_3)

#then make venn diagram 
venn_all_diseases_3 <- ggVennDiagram(x = CPP_list_3)+ scale_fill_gradient(low ="810f7c", high = "#8c96c6")+labs(title = "Rectal Samples") + theme(plot.title = element_text(hjust = 0.5)) + theme(plot.title = element_text(size=20, face="bold"))
venn_all_diseases_3

#prune 
prune_taxa(CPP_only_ASVs_3,phyloseq) %>%
  tax_table()

prune_taxa(CPP_endo_ASVs_3,phyloseq) %>%
  tax_table()

prune_taxa(CPP_healthy_ASVs_3,phyloseq) %>%
  tax_table()

#plot at the genus level - gives the samples within each group that contain the ASVs only present within CPP only and abundance 
prune_taxa(CPP_only_ASVs_3,phyloseq) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`Host_disease`, scales ="free")

#plot at the genus level - gives the samples within each group that contain ASVs only present within CPP Endo and abundance 
prune_taxa(CPP_endo_ASVs_3,phyloseq) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`Host_disease`, scales ="free")

##plot at the genus level - gives the samples within each group that contain ASVs only present within CPP healthy and abundance 
prune_taxa(CPP_healthy_ASVs_3,phyloseq) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`Host_disease`, scales ="free")

#try stratifying by the vaginal sample and disease instead now  
CPP_only_3 <- subset_samples(phyloseq, `Host_disease`=="CPP")
CPP_endo_3 <- subset_samples(phyloseq, `Host_disease`=="CPP Endo")
CPP_healthy_3 <-subset_samples(phyloseq, `Host_disease`=="Control")

CPP_only_4 <- subset_samples(CPP_only_2,`collection_method` == "vaginal_swab")
CPP_endo_4 <- subset_samples(CPP_endo_2,`collection_method` == "vaginal_swab")
CPP_healthy_4 <-subset_samples(CPP_healthy_2, `collection_method` == "vaginal_swab")

#identify thresholds of prevalence/abundance 
CPP_only_ASVs_4 <- core_members(CPP_only_4, detection=0, prevalence = 0.5)
CPP_endo_ASVs_4 <- core_members(CPP_endo_4, detection=0, prevalence = 0.5)
CPP_healthy_ASVs_4 <- core_members(CPP_healthy_4, detection=0, prevalence = 0.5)

#create list 
CPP_list_4 <- list(Healthy = CPP_healthy_ASVs_4, CPP_Only = CPP_only_ASVs_4, CPP_Endo = CPP_endo_ASVs_4)

#create a venn diagram 
venn_all_diseases_4 <- ggVennDiagram(x = CPP_list_4)+ scale_fill_gradient(low ="810f7c", high = "#8c96c6")+labs(title = "Vaginal Samples") + theme(plot.title = element_text(hjust = 0.5)) + theme(plot.title = element_text(size=20, face="bold"))
venn_all_diseases_4
 
#prune 
prune_taxa(CPP_only_ASVs_4,phyloseq) %>%
  tax_table()

prune_taxa(CPP_endo_ASVs_4,phyloseq) %>%
  tax_table()

prune_taxa(CPP_healthy_ASVs_4,phyloseq) %>%
  tax_table()

#plot at the genus level - gives the samples within each group that contain the ASVs only present within CPP only and abundance 
prune_taxa(CPP_only_ASVs_4,phyloseq) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`Host_disease`, scales ="free")

#plot at the genus level - gives the samples within each group that contain ASVs only present within CPP Endo and abundance 
prune_taxa(CPP_endo_ASVs_4,phyloseq) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`Host_disease`, scales ="free")

##plot at the genus level - gives the samples within each group that contain ASVs only present within CPP healthy and abundance 
prune_taxa(CPP_healthy_ASVs_4,phyloseq) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`Host_disease`, scales ="free")
