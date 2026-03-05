install.packages("devtools")
devtools::install_github("cafferychen777/ggpicrust2")
install.packages("MicrobiomeStat")
install.packages("GGally")
install.packages("phyloseq")

library(tidyverse)
library(phyloseq)
library(ggpicrust2)

ps <- readRDS("data/phyloseq_filtered.rds")

meta <- sample_data(ps) %>%
  data.frame() %>%
  rownames_to_column("sample_name")
ko <- read.delim("data/pred_metagenome_unstrat.tsv")
head(ko[,1:6])
