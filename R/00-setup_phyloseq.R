# MICB475 Team 12
# 00-setup_phyloseq.R
# Import data and construct phyloseq object


library(phyloseq)
library(tidyverse)
library(ape)
library(vegan)


meta <- read_delim("data/cpp_meta.tsv", delim = "\t")
otu <- read_delim("data/feature-table.txt", delim = "\t", skip = 1)
tax <- read_delim("data/taxonomy.tsv", delim = "\t")
phylotree <- read.tree("data/tree.nwk")


otu_mat <- as.matrix(otu[, -1])
rownames(otu_mat) <- otu$`#OTU ID`

OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)
class(OTU)


### Format sample metadata

samp_df <- as.data.frame(meta[,-1])
rownames(samp_df)<- meta$'sample_name'

SAMP <- sample_data(samp_df)
class(SAMP)


#### Formatting taxonomy ####

tax_mat <- tax %>% select(-Confidence)%>%
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix() 

tax_mat <- tax_mat[,-1]
rownames(tax_mat) <- tax$`Feature ID`
TAX <- tax_table(tax_mat)
class(TAX)

### Create phyloseq object
## mouse_phylo <- phyloseq(OTU, SAMP, TAX, phylotree)

