library(tidyverse)
library(phyloseq)
library(indicspecies)
set.seed(2026)
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
taxtable_vaginal <- tax_table(phyloseq) %>% 
  as.data.frame() %>% 
  rownames_to_column(var="ASV")%>% 
  mutate(
    Genus = Genus %>%
      str_replace_all("g__", "")   # remove g prefix
  )

#Join with taxatable based on ASV ID 
#Filter for anything with p value less than 0.05 and view
CPP_table_vaginal <- isa_CPP_vaginal$sign %>%
  rownames_to_column(var="ASV") %>%
  mutate(
    Host_disease = case_when(
      s.Control == 1 ~ "Control",
      s.CPP == 1 ~ "CPP",
      `s.CPP Endo` == 1 ~ "CPP Endo",
      TRUE ~ "Mixed"
    )
  ) %>%
  left_join(taxtable_vaginal) %>%
  filter(p.value<0.05) 

View(CPP_table_vaginal)

#Make a plot

indicator_plot_vaginal <- ggplot(CPP_table_vaginal, aes(x = reorder(Genus, stat), y = stat, fill = Host_disease)) +
  geom_col() +
  coord_flip() +
  labs(x = "Genus", y = "Indicator Value", fill = "Associated Condition") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 18),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 18),
    strip.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    plot.margin = margin(15, 15, 15, 15)
  )

indicator_plot_vaginal

ggsave("results/aim2/04b-indicator_species/indicator_plot_vaginal.png", indicator_plot_vaginal, width = 10, height = 7)

#### Indicator Species/Taxa Analysis for Rectal ####

#Filter only rectal samples 
CPP_genus_rectal<- subset_samples(CPP_genus, `collection_method`=="rectal_swab")

#Convert to relative abundance 
CPP_genus_rectal_RA <- transform_sample_counts(CPP_genus_rectal, fun=function(x) x/sum(x))

#ISA using host disease as the predictor 
isa_CPP_rectal <- multipatt(t(otu_table(CPP_genus_rectal_RA)), cluster = sample_data(CPP_genus_rectal_RA)$`Host_disease`)
summary(isa_CPP_rectal)

#Extract taxa table from phloseq as a data frame, have the ASV as the row name 
taxtable_rectal <- tax_table(phyloseq) %>% 
  as.data.frame() %>% 
  rownames_to_column(var="ASV")  %>% 
  mutate(
    Genus = Genus %>%
      str_replace_all("g__", "")   # remove g prefix
  )

#Join with taxatable based on ASV ID 
#Filter for anything with p value less than 0.05 and view

CPP_table_rectal <- isa_CPP_rectal$sign %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV") %>%
  mutate(
    Host_disease = apply(
      select(., s.Control, s.CPP, `s.CPP Endo`),
      1,
      function(x) paste(names(x)[x == 1], collapse = " + ")
    )
  ) %>%
  mutate(
    `Associated Conditions` = Host_disease %>%
      str_replace_all("s.", "")   # remove s. prefix
  ) %>% 
  left_join(taxtable_rectal, by = "ASV") %>%
  filter(p.value < 0.05)


view(CPP_table_rectal)

#Make a plot

indicator_plot_rectal <- ggplot(CPP_table_rectal, aes(x = reorder(Genus, stat), y = stat, fill = `Associated Conditions`)) +
  geom_col() +
  coord_flip() +
  labs(x = "Genus", y = "Indicator Value") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 18),
    strip.text = element_text(size = 18),
    strip.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    plot.margin = margin(15, 15, 15, 15),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )

indicator_plot_rectal

ggsave("results/aim2/04b-indicator_species/indicator_plot_rectal.png", indicator_plot_rectal, width = 10, height = 7)

