library(ggpicrust2)
library(DESeq2)
library(edgeR)


ps <- readRDS("data/phyloseq_filtered.rds")

meta <- sample_data(ps) %>%
  data.frame() %>%
  rownames_to_column("sample_name")

meta_cpp_cpp_endo <- meta |> filter(Host_disease != "Control")
meta_cpp_control <- meta |> filter(Host_disease != "CPP Endo")
meta_cpp_endo_control <- meta |> filter(Host_disease != "CPP")

kegg_abundance <- ko2kegg_abundance("data/picrust_out/KO_metagenome_out/pred_metagenome_unstrat.tsv")


## CPP vs CPP Endo
rectal_meta_cpp_cpp_endo <- meta_cpp_cpp_endo %>%
  filter(env_medium == "rectal")
vaginal_meta_cpp_cpp_endo <- meta_cpp_cpp_endo %>%
  filter(env_medium == "vaginal")

daa_results_df_cpp_cpp_endo_rect <- pathway_daa(abundance = kegg_abundance, metadata = rectal_meta_cpp_cpp_endo, group = "Host_disease", daa_method = "edgeR")
daa_results_df_cpp_cpp_endo_vag <- pathway_daa(abundance = kegg_abundance, metadata = vaginal_meta_cpp_cpp_endo, group = "Host_disease", daa_method = "edgeR")


