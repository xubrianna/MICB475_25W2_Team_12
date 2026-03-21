library(tidyverse)
library(phyloseq)
library(indicspecies)

#### Load data ####
phyloseq <- readRDS('data/phyloseq_filtered.rds')

otu_table(phyloseq)
sample_data(phyloseq)

#### Indicator Species/Taxa Analysis ####

#Filter to the genus level 
CPP_genus <- tax_glom(phyloseq, "Genus", NArm = FALSE)

#### Indicator Species/Taxa Analysis Vaginal ####

#Filter only vaginal samples 
CPP_genus_vaginal <- subset_samples(CPP_genus, `collection_method`=="vaginal_swab")

#Convert to relative abundance 
CPP_genus_vaginal_RA <- transform_sample_counts(CPP_genus_vaginal, fun=function(x) x/sum(x))

#ISA using host disease as the predictor 
isa_CPP_vaginal <- multipatt(t(otu_table(CPP_genus_vaginal_RA)), cluster = sample_data(CPP_genus_vaginal_RA)$`Host_disease`)
summary(isa_CPP_vaginal)

#Extract taxa table from phloseq as a data frame, have the ASV as the row name 
taxtable_vaginal <- tax_table(phyloseq) %>% as.data.frame() %>% rownames_to_column(var="ASV")

#Join with taxatable based on ASV ID 
#Filter for anything with p value less than 0.05 and view

CPP_table_vaginal <- isa_CPP_vaginal$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable_vaginal) %>%
  filter(p.value<0.05) 

view(CPP_table_vaginal)

#Make a plot

indicator_plot_vaginal <- ggplot(CPP_table_vaginal, aes(x = reorder(Genus, stat), y = stat, fill = Genus)) +
  geom_col() +
  coord_flip() +
  labs(x = "Genus", y = "Indicator Value", title = "Significant Indicator Species") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 25),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none",
    strip.text = element_text(size = 14),
    plot.margin = margin(15, 15, 15, 15)
  )

ggsave("results/aim2/06-indicator_species/indicator_plot_vaginal.png", indicator_plot_vaginal, width = 10, height = 7)

#### Indicator Species/Taxa Analysis for Rectal ####

#Filter only rectal samples 
CPP_genus_rectal<- subset_samples(CPP_genus, `collection_method`=="rectal_swab")

#Convert to relative abundance 
CPP_genus_rectal_RA <- transform_sample_counts(CPP_genus_rectal, fun=function(x) x/sum(x))

#ISA using host disease as the predictor 
isa_CPP_rectal <- multipatt(t(otu_table(CPP_genus_rectal_RA)), cluster = sample_data(CPP_genus_rectal_RA)$`Host_disease`)
summary(isa_CPP_rectal)

#Extract taxa table from phloseq as a data frame, have the ASV as the row name 
taxtable_rectal <- tax_table(phyloseq) %>% as.data.frame() %>% rownames_to_column(var="ASV")

#Join with taxatable based on ASV ID 
#Filter for anything with p value less than 0.05 and view

CPP_table_rectal <- isa_CPP_rectal$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable_rectal) %>%
  filter(p.value<0.05) 

view(CPP_table_rectal)

#Make a plot

indicator_plot_rectal <- ggplot(CPP_table_rectal, aes(x = reorder(Genus, stat), y = stat, fill = Genus)) +
  geom_col() +
  coord_flip() +
  labs(x = "Genus", y = "Indicator Value", title = "Significant Indicator Species") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 25),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none",
    strip.text = element_text(size = 14),
    plot.margin = margin(15, 15, 15, 15)
  )

indicator_plot_rectal

ggsave("results/aim2/06-indicator_species/indicator_plot_rectal.png", indicator_plot_rectal, width = 10, height = 7)
