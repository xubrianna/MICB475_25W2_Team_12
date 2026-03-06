
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)

##LOAD AND MODIFY PHYLOSEQ
#load filtered phyloseq without rarification (expects uneven sequencing depths)
phyloseq <- readRDS('data/phyloseq_filtered.rds')

#transform OTU table into relative abundance from absolute counts
CPP_RA <- transform_sample_counts(phyloseq, fun=function(x) x/sum(x))

#view 
otu_table(phyloseq)
sample_data(phyloseq)

## STRATIFY BY BODY SITE 
# RECTAL
#try stratifying by the rectal sample after stratify by disease 
CPP_only <- subset_samples(CPP_RA, `Host_disease`=="CPP")
CPP_endo <- subset_samples(CPP_RA, `Host_disease`=="CPP Endo")
CPP_healthy <-subset_samples(CPP_RA, `Host_disease`=="Control")

CPP_only_rectal <- subset_samples(CPP_only,`collection_method` == "rectal_swab")
CPP_endo_rectal <- subset_samples(CPP_endo,`collection_method` == "rectal_swab")
CPP_healthy_rectal <-subset_samples(CPP_healthy, `collection_method` == "rectal_swab")

#then assign ASVs - 0.001 filters out rare
CPP_only_ASVs_rectal_0.001 <- core_members(CPP_only_rectal, detection=0.001, prevalence = 0.5)
CPP_endo_ASVs_rectal_0.001 <- core_members(CPP_endo_rectal, detection=0.001, prevalence = 0.5)
CPP_healthy_ASVs_rectal_0.001 <- core_members(CPP_healthy_rectal, detection=0.001, prevalence = 0.5)

#then make a list
CPP_list_rectal_0.001 <- list(Healthy = CPP_healthy_ASVs_rectal_0.001, CPP_Only = CPP_only_ASVs_rectal_0.001, CPP_Endo = CPP_endo_ASVs_rectal_0.001)

#then make venn diagram 
venn_all_diseases_rectal_0.001 <- ggVennDiagram(x = CPP_list_rectal_0.001)+ scale_fill_gradient(low ="#c7e9b4", high = "#225ea8")+
  labs(title = "Rectal Samples") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(plot.title = element_text(size=20, face="bold"))
venn_all_diseases_rectal_0.001

ggsave("results/aim2/venn_all_diseases_rectal_0.001.png", venn_all_diseases_rectal_0.001)

#prune 
prune_taxa(CPP_only_ASVs_rectal_0.001,CPP_RA) %>%
  tax_table()

prune_taxa(CPP_endo_ASVs_rectal_0.001,CPP_RA) %>%
  tax_table()

prune_taxa(CPP_healthy_ASVs_rectal_0.001,CPP_RA) %>%
  tax_table()

#plot at the genus level - gives the samples within each group that contain the ASVs only present within CPP only and abundance 
prune_taxa(CPP_only_ASVs_rectal_0.001,CPP_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`Host_disease`, scales ="free")

#plot at the genus level - gives the samples within each group that contain ASVs only present within CPP Endo and abundance 
prune_taxa(CPP_endo_ASVs_rectal_0.001,CPP_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`Host_disease`, scales ="free")

##plot at the genus level - gives the samples within each group that contain ASVs only present within CPP healthy and abundance 
prune_taxa(CPP_healthy_ASVs_rectal_0.001,CPP_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`Host_disease`, scales ="free")

#try abundance of 0 - absence/presence 
CPP_only_ASVs_rectal_0 <- core_members(CPP_only_rectal, detection=0, prevalence = 0.5)
CPP_endo_ASVs_rectal_0 <- core_members(CPP_endo_rectal, detection=0, prevalence = 0.5)
CPP_healthy_ASVs_rectal_0 <- core_members(CPP_healthy_rectal, detection=0, prevalence = 0.5)

CPP_list_rectal_0 <- list(Healthy = CPP_healthy_ASVs_rectal_0, CPP_Only = CPP_only_ASVs_rectal_0, CPP_Endo = CPP_endo_ASVs_rectal_0)

#then make venn diagram
venn_all_diseases_rectal_0<- ggVennDiagram(x = CPP_list_rectal_0)+ scale_fill_gradient(low ="#c7e9b4", high = "#225ea8") +
  labs(title = "Rectal Samples") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(plot.title = element_text(size=20, face="bold"))
venn_all_diseases_rectal_0

ggsave("results/aim2/venn_all_diseases_rectal_0.png", venn_all_diseases_rectal_0)

#try abundance of 0.01 to keep only abundant 
CPP_only_ASVs_rectal_0.01 <- core_members(CPP_only_rectal, detection=0.01, prevalence = 0.5)
CPP_endo_ASVs_rectal_0.01 <- core_members(CPP_endo_rectal, detection=0.01, prevalence = 0.5)
CPP_healthy_ASVs_rectal_0.01 <- core_members(CPP_healthy_rectal, detection=0.01, prevalence = 0.5)

CPP_list_rectal_0.01 <- list(Healthy = CPP_healthy_ASVs_rectal_0.01, CPP_Only = CPP_only_ASVs_rectal_0.01, CPP_Endo = CPP_endo_ASVs_rectal_0.01)

#then make venn diagram
venn_all_diseases_rectal_0.01<- ggVennDiagram(x = CPP_list_rectal_0.01)+ scale_fill_gradient(low ="#c7e9b4", high = "#225ea8") +
  labs(title = "Rectal Samples") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(plot.title = element_text(size=20, face="bold"))
venn_all_diseases_rectal_0.01

ggsave("results/aim2/venn_all_diseases_rectal_0.01.png", venn_all_diseases_rectal_0.01)


#VAGINAL
#try stratifying by the vaginal sample instead now  
CPP_only_vaginal <- subset_samples(CPP_only,`collection_method` == "vaginal_swab")
CPP_endo_vaginal <- subset_samples(CPP_endo,`collection_method` == "vaginal_swab")
CPP_healthy_vaginal <-subset_samples(CPP_healthy, `collection_method` == "vaginal_swab")

#identify thresholds of prevalence/abundance 
CPP_only_ASVs_vaginal_0.001 <- core_members(CPP_only_vaginal, detection=0.001, prevalence = 0.5)
CPP_endo_ASVs_vaginal_0.001 <- core_members(CPP_endo_vaginal, detection=0.001, prevalence = 0.5)
CPP_healthy_ASVs_vaginal_0.001 <- core_members(CPP_healthy_vaginal, detection=0.001, prevalence = 0.5)

#create list 
CPP_list_vaginal_0.001 <- list(Healthy = CPP_healthy_ASVs_vaginal_0.001, CPP_Only = CPP_only_ASVs_vaginal_0.001, CPP_Endo = CPP_endo_ASVs_vaginal_0.001)

#create a venn diagram 
venn_all_diseases_vaginal_0.001 <- ggVennDiagram(x = CPP_list_vaginal_0.001)+ scale_fill_gradient(low ="#c7e9b4", high = "#225ea8")+
  labs(title = "Vaginal Samples") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(plot.title = element_text(size=20, face="bold"))
venn_all_diseases_vaginal_0.001

ggsave("results/aim2/venn_all_diseases_vaginal_0.001.png", venn_all_diseases_vaginal_0.001)
 
#prune 
prune_taxa(CPP_only_ASVs_vaginal_0.001,CPP_RA) %>%
  tax_table()

prune_taxa(CPP_endo_ASVs_vaginal_0.001,CPP_RA) %>%
  tax_table()

prune_taxa(CPP_healthy_ASVs_vaginal_0.001,CPP_RA) %>%
  tax_table()

#plot at the genus level - gives the samples within each group that contain the ASVs only present within CPP only and abundance 
prune_taxa(CPP_only_ASVs_vaginal_0.001,CPP_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`Host_disease`, scales ="free")

#plot at the genus level - gives the samples within each group that contain ASVs only present within CPP Endo and abundance 
prune_taxa(CPP_endo_ASVs_vaginal_001,CPP_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`Host_disease`, scales ="free")

##plot at the genus level - gives the samples within each group that contain ASVs only present within CPP healthy and abundance 
prune_taxa(CPP_healthy_ASVs_vaginal_0.001,CPP_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`Host_disease`, scales ="free")

#try abundance of 0 
CPP_only_ASVs_vaginal_0 <- core_members(CPP_only_vaginal, detection=0, prevalence = 0.5)
CPP_endo_ASVs_vaginal_0 <- core_members(CPP_endo_vaginal, detection=0, prevalence = 0.5)
CPP_healthy_ASVs_vaginal_0 <- core_members(CPP_healthy_vaginal, detection=0, prevalence = 0.5)

#create list 
CPP_list_vaginal_0 <- list(Healthy = CPP_healthy_ASVs_vaginal_0, CPP_Only = CPP_only_ASVs_vaginal_0, CPP_Endo = CPP_endo_ASVs_vaginal_0)

#create a venn diagram 
venn_all_diseases_vaginal_0 <- ggVennDiagram(x = CPP_list_vaginal_0)+ scale_fill_gradient(low ="#c7e9b4", high = "#225ea8")+
  labs(title = "Vaginal Samples") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(plot.title = element_text(size=20, face="bold"))
venn_all_diseases_vaginal_0

ggsave("results/aim2/venn_all_diseases_vaginal_0.png", venn_all_diseases_vaginal_0)

#try abundance of 0.01 
CPP_only_ASVs_vaginal_0.01 <- core_members(CPP_only_vaginal, detection=0.01, prevalence = 0.5)
CPP_endo_ASVs_vaginal_0.01 <- core_members(CPP_endo_vaginal, detection=0.01, prevalence = 0.5)
CPP_healthy_ASVs_vaginal_0.01 <- core_members(CPP_healthy_vaginal, detection=0.01, prevalence = 0.5)

#create list 
CPP_list_vaginal_0.01 <- list(Healthy = CPP_healthy_ASVs_vaginal_0.01, CPP_Only = CPP_only_ASVs_vaginal_0.01, CPP_Endo = CPP_endo_ASVs_vaginal_0.01)

#create a venn diagram 
venn_all_diseases_vaginal_0.01 <- ggVennDiagram(x = CPP_list_vaginal_0.01)+ scale_fill_gradient(low ="#c7e9b4", high = "#225ea8")+
  labs(title = "Vaginal Samples") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(plot.title = element_text(size=20, face="bold"))
venn_all_diseases_vaginal_0.01

ggsave("results/aim2/venn_all_diseases_vaginal_0.01.png", venn_all_diseases_vaginal_0.01)