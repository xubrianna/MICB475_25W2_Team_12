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

#####load packages####
library(phyloseq)
library(vegan)
library(ggplot2)
library(ape)
library(tidyverse)
library(FSA)
library(picante)

####otu table#####
otu <- read.table("feature-table.txt",header=TRUE,sep="\t",skip=1,row.names=1)
otu <- as.matrix(otu)

####metadata####
meta <- read.table("cpp_metadata.tsv", header=TRUE, sep="\t", row.names=1)

####tree####
tree <- read.tree("tree.nwk")

####taxonomy####
tax <- read_delim("taxonomy.tsv", delim = "\t")
tax_mat <- tax %>% select(-Confidence) %>%
  separate(col=Taxon, sep=";",
           into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix()
rownames(tax_mat) <- tax$`Feature ID`
TAX <- tax_table(tax_mat)

####phyloseq#####
#correct OTU to match w/ meta
OTU <- otu_table(otu, taxa_are_rows=TRUE)
sample_names(OTU) <- gsub("^X", "", sample_names(OTU))
sample_names(OTU) <- gsub("\\.0$", "", sample_names(OTU))
sample_names(OTU) <- sample_names(META)

TREE <- phy_tree(tree)
ps <- phyloseq(OTU, sample_data(META), tax_table(TAX), phy_tree(tree))

####phyloseq qc####
#omit mitochondrial/chloroplast sequences
ps_filt <- subset_taxa(ps,
                       Domain == "d__Bacteria" &
                         !Class %in% c("c__Chloroplast") &
                         !Family %in% c("f__Mitochondria"))
#filter low abundance ASVs, less than 5
ps_filt_nolow <- filter_taxa(ps_filt, function(x) sum(x) > 5, prune = TRUE)

#filter for low reads
ps_final <- prune_samples(sample_sums(ps_filt_nolow) > 100, ps_filt_nolow)

#rarefy samples
set.seed(123)
ps_rarefied <- rarefy_even_depth(ps_final, sample.size = 43000, rngseed = 123,
                                 replace = FALSE, verbose = TRUE)

####diversity metrics####
# 1. Extract OTU matrix
otu_mat <- as(otu_table(ps_rarefied), "matrix")
if(taxa_are_rows(ps_rarefied)) {
  otu_mat <- t(otu_mat)
}

# 2. Calculate Faith's PD
alpha_faith <- pd(otu_mat, phy_tree(ps_rarefied), include.root = TRUE)

# 3. Add sample names
alpha_faith$sample <- rownames(alpha_faith)

# 4. Merge with metadata (use your existing META)
# Ensure sample names match
alpha_faith$sample <- trimws(alpha_faith$sample)
META$sample <- trimws(META$sample)
META$env_medium <- as.character(META$env_medium)  # ensure character

alpha_faith <- alpha_faith %>%
  left_join(as.data.frame(META) %>% rownames_to_column(var="sample"), by="sample") %>%
  filter(!is.na(Host_disease) & !is.na(env_medium))  # keep only samples with both

# 5. Check which samples remain
table(alpha_faith$env_medium, useNA="ifany")

# 6. Factor env_medium for consistent facet order
alpha_faith$env_medium <- factor(alpha_faith$env_medium,
                                 levels = c("Rectal","Vaginal"))

# 7. Plot Faith's PD by Host Disease, faceted by env_medium
ggplot(alpha_faith, aes(x=Host_disease, y=PD, fill=Host_disease)) +
  geom_boxplot(lwd=0.8) +
  geom_jitter(width=0.2, size=3) +
  facet_wrap(~env_medium) +
  theme_bw() +
  ylab("Faith's Phylogenetic Diversity") +
  xlab("Host Disease") +
  scale_fill_manual(values=c(
    "Healthy"="green",
    "CPP-only"="orange",
    "CPP-endo"="red"
  ))

# 8. Optional: Kruskal-Wallis per env_medium
kruskal.test(PD ~ Host_disease, data = subset(alpha_faith, env_medium=="Vaginal"))
kruskal.test(PD ~ Host_disease, data = subset(alpha_faith, env_medium=="Rectal"))

# 9. Optional: Dunn's test pairwise comparisons
dunnTest(PD ~ Host_disease, data = subset(alpha_faith, env_medium=="Vaginal"), method="bh")
dunnTest(PD ~ Host_disease, data = subset(alpha_faith, env_medium=="Rectal"), method="bh")

#beta diversity
unifrac_dm <- distance(ps_rarefied, method = "unifrac")  # unweighted by default

#pcoa
pcoa_uni <- ordinate(ps_rarefied, method = "PCoA", distance = unifrac_dm)

library(ggplot2)
gg_uni <- plot_ordination(ps_rarefied, pcoa_uni, type="samples",
                          color="days_post_transplant",
                          shape="cage_id") +
  geom_point(size=4) +
  scale_color_gradient(low="yellow", high="blue") +
  theme_bw()
gg_uni









