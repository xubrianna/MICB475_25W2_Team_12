# MICB475 Team 12
# 07-random_forest.R
# Random Forest classification of Host_disease (Control vs CPP vs CPP Endo)
# Separate models for rectal and vaginal samples

library(randomForest)
library(caret)
library(ranger)
library(pROC)
library(ggplot2)
library(phyloseq)
library(tidyverse)
library(cowplot)

set.seed(2026)

# =============================================================================
# Part 1: Load data
# =============================================================================

ps <- readRDS("data/phyloseq_filtered.rds") %>%
  tax_glom("Genus")

# =============================================================================
# Part 2: Helper function to prepare data for RF
# =============================================================================

avg_abundance = taxa_sums(ps)/sum(taxa_sums(ps)) 
# Sort high to low
avg_abundance = sort(avg_abundance, decreasing = T)

avg_abundance = names(avg_abundance)
# Take the top 10
top_10 = avg_abundance[1:10]
# Extract taxa names 
top_10 = names(top_10)

prepare_rf_data <- function(ps_subset, prev_threshold = 0.1) {
  # CLR transform on all retained taxa
  ps_clr <- ps_subset %>% microbiome::transform("clr")
  
  ps_clr <- prune_taxa(avg_abundance,ps_clr)
  
  # Melt and create unique feature names per OTU
  df <- psmelt(ps_clr)
  
  
  
  # Build a lookup: one clean name per OTU (unique across the dataset)
  otu_names <- df %>%
    distinct(OTU, Genus) %>%
    mutate(feature_name = ifelse(is.na(Genus) | Genus == "",
                                paste0("OTU_", OTU),
                                paste0(Genus, "_", OTU))) %>%
    mutate(feature_name = make.names(feature_name, unique = TRUE))
  
  df <- df %>%
    left_join(otu_names %>% select(OTU, feature_name), by = "OTU") %>%
    group_by(feature_name) %>%
    mutate(Abundance = as.numeric(scale(Abundance))) %>%
    ungroup()
  
  # Pivot to wide format: one column per feature
  df_pivot <- df %>%
    select(Sample, Host_disease, feature_name, Abundance) %>%
    pivot_wider(names_from = feature_name, values_from = Abundance)
  
  # # Pivot to wide format: one column per feature
  # df_pivot <- df %>%
  #   select(Sample, Host_disease, all_of(meta_vars), feature_name, Abundance) %>%
  #   pivot_wider(names_from = feature_name, values_from = Abundance)
  # 
  # # Ensure metadata columns are numeric
  # df_pivot <- df_pivot %>%
  #   mutate(across(all_of(meta_vars), as.numeric))
  
  # Remove NAs and sample ID column
  df_final <- df_pivot %>%
    na.omit() %>%
    select(-Sample)
  
  return(df_final)
}

# =============================================================================
# Part 3: Helper function to run RF with k-fold CV + hyperparameter tuning
# =============================================================================

run_rf <- function(df, outcome_col = "Host_disease", k = 5,
                   repeats = 10, seed = 421) {
  set.seed(seed)
  
  # Separate predictors and outcome
  predictors <- df %>% select(-all_of(outcome_col))
  outcome <- df %>% pull(all_of(outcome_col)) %>%
    factor(levels = c("Control", "CPP Only", "CPP Endo"))
  
  # Hyperparameter grid
  tune_grid <- expand.grid(
    mtry = c(2, 3, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15),
    splitrule = c("gini", "extratrees"),
    min.node.size = c(1, 2, 3, 5, 10)
  )
  
  # Cap mtry at number of predictors
  tune_grid <- tune_grid %>% filter(mtry <= ncol(predictors))
  
  # Repeated stratified k-fold CV for more stable estimates
  train_control <- trainControl(
    method = "repeatedcv",
    number = k,
    repeats = repeats,
    classProbs = TRUE,
    summaryFunction = multiClassSummary,
    savePredictions = "final"
  )
  
  # Make valid R names for factor levels (caret requirement)
  levels(outcome) <- make.names(levels(outcome))
  
  # Run RF
  rf_model <- train(
    x = predictors,
    y = outcome,
    method = "ranger",
    trControl = train_control,
    tuneGrid = tune_grid,
    metric = "AUC",
    importance = "impurity",
    num.trees = 1000
  )
  
  return(rf_model)
}

# =============================================================================
# Part 4: Subset by body site and prepare data
# =============================================================================

ps_rectal <- subset_samples(ps, env_medium == "rectal")
ps_vaginal <- subset_samples(ps, env_medium == "vaginal")

# Remove zero-abundance taxa after subsetting
ps_rectal <- prune_taxa(taxa_sums(ps_rectal) > 0, ps_rectal)
ps_vaginal <- prune_taxa(taxa_sums(ps_vaginal) > 0, ps_vaginal)

df_rectal <- prepare_rf_data(ps_rectal)
df_vaginal <- prepare_rf_data(ps_vaginal)

# =============================================================================
# Part 5: Run RF models
# =============================================================================

cat("=== Running RF for RECTAL samples ===\n")
rf_rectal <- run_rf(df_rectal, k = 5, seed = 123)
print(rf_rectal)

cat("\n=== Running RF for VAGINAL samples ===\n")
rf_vaginal <- run_rf(df_vaginal, k = 5, seed = 123)
print(rf_vaginal)

# =============================================================================
# Part 6: ROC curves (one-vs-rest for each class)
# =============================================================================

plot_multiclass_roc <- function(model, title) {
  preds <- model$pred
  
  best <- model$bestTune
  for (param in names(best)) {
    preds <- preds[preds[[param]] == best[[param]], ]
  }
  
  classes <- levels(preds$obs)
  roc_list <- list()
  auc_vals <- c()
  
  for (cls in classes) {
    binary_truth <- ifelse(preds$obs == cls, 1, 0)
    roc_obj <- pROC::roc(binary_truth, preds[[cls]], quiet = TRUE)
    roc_list[[cls]] <- roc_obj
    auc_vals[cls] <- as.numeric(pROC::auc(roc_obj))
  }
  
  roc_df <- bind_rows(lapply(classes, function(cls) {
    r <- roc_list[[cls]]
    data.frame(
      FPR = 1 - r$specificities,
      TPR = r$sensitivities,
      Class = cls
    )
  }))
  
  group_colors <- c(
    "Control" = "#c7e9b4",
    "CPP.Only" = "#41b6c4",
    "CPP.Endo" = "#225ea8"
  )
  
  roc_labels <- c(
    "Control" = paste0("Control (AUC=", round(auc_vals["Control"], 2), ")"),
    "CPP.Only" = paste0("CPP Only (AUC=", round(auc_vals["CPP.Only"], 2), ")"),
    "CPP.Endo" = paste0("CPP Endo (AUC=", round(auc_vals["CPP.Endo"], 2), ")")
  )
  
  print(unique(roc_df$Class))
  print(names(group_colors))
  
  p <- ggplot(roc_df, aes(x = FPR, y = TPR, colour = Class)) +
    geom_line(linewidth = 1) +
    geom_abline(
      slope = 1, intercept = 0,
      color = "gray", linetype = "dashed", linewidth = 0.8
    ) +
    scale_color_manual(
      values = group_colors,
      breaks = names(group_colors),
      labels = roc_labels
    ) +
    labs(
      x = "False Positive Rate",
      y = "True Positive Rate",
      title = title,
      color = "Host Disease"
    ) +
    theme_classic(base_size = 16) +
    theme(
      legend.position = "right",
      legend.direction = "vertical",
      axis.title = element_text(size=18),
    )
  
  return(list(plot = p, auc = auc_vals))
}


roc_rectal <- plot_multiclass_roc(rf_rectal, "")
roc_vaginal <- plot_multiclass_roc(rf_vaginal, "")

roc_rectal$plot
roc_vaginal$plot

ggsave("results/aim4/07-rf_roc_rectal.png", roc_rectal$plot,
       width = 8, height = 7, dpi = 300)
ggsave("results/aim4/07-rf_roc_vaginal.png", roc_vaginal$plot,
       width = 8, height = 7, dpi = 300)

# Print AUC values
cat("\n=== AUC values (Rectal) ===\n")
print(round(roc_rectal$auc, 3))
cat("\n=== AUC values (Vaginal) ===\n")
print(round(roc_vaginal$auc, 3))

# =============================================================================
# Part 7: Feature importance
# =============================================================================

plot_importance <- function(model, title) {
  imp <- varImp(model)$importance
  imp$Feature <- rownames(imp)
  
  imp <- imp %>%
    mutate(Importance = Overall) %>%
    mutate(Feature = gsub("g__", "", Feature)) %>%
    mutate(Feature = gsub("_.*", "", Feature)) %>%
    group_by(Feature) %>%                          # ← collapse duplicates
    summarise(Importance = mean(Importance), .groups = "drop") %>%
    arrange(desc(Importance)) %>%
    slice_max(Importance, n = 20) %>%
    mutate(Feature = factor(Feature, levels = Feature))  # ← now unique
  
  p <- ggplot(imp, aes(x = Feature, y = Importance, fill = Importance)) +
    geom_col() +
    theme_classic(base_size = 16) +
    theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)) +
    ylab("Importance (Gini)") + xlab(NULL) +
    ggtitle(title)
  
  return(list(plot = p, data = imp))
}

imp_rectal <- plot_importance(rf_rectal, "Feature Importance - Rectal")
imp_vaginal <- plot_importance(rf_vaginal, "Feature Importance - Vaginal")

imp_rectal$plot
imp_vaginal$plot


ggsave("results/aim4/07-rf_importance_rectal.png", imp_rectal$plot,
       width = 18, height = 10, dpi = 300)
ggsave("results/aim4/07-rf_importance_vaginal.png", imp_vaginal$plot,
       width = 18, height = 10, dpi = 300)

# =============================================================================
# Part 8: Summary table
# =============================================================================

results_table <- data.frame(
  Body_Site = c("Rectal", "Vaginal"),
  Best_mtry = c(rf_rectal$bestTune$mtry, rf_vaginal$bestTune$mtry),
  Best_splitrule = c(as.character(rf_rectal$bestTune$splitrule),
                     as.character(rf_vaginal$bestTune$splitrule)),
  Best_min_node = c(rf_rectal$bestTune$min.node.size,
                    rf_vaginal$bestTune$min.node.size),
  CV_Accuracy = c(max(rf_rectal$results$Accuracy),
                  max(rf_vaginal$results$Accuracy)),
  AUC_Control = c(roc_rectal$auc["Control"], roc_vaginal$auc["Control"]),
  AUC_CPP = c(roc_rectal$auc["CPP.Only"], roc_vaginal$auc["CPP.Only"]),
  AUC_CPP_Endo = c(roc_rectal$auc["CPP.Endo"], roc_vaginal$auc["CPP.Endo"])
)

print(results_table)
write_tsv(results_table, "results/aim4/07-rf_summary.tsv")


tax1 <- plot_grid(roc_rectal$plot + ggtitle("Rectal") + theme(plot.title = element_text(hjust = 0.5, size= 18),
                                                               plot.margin = ggplot2::margin(15, 15, 25, 15)), 
                   roc_vaginal$plot + ggtitle("Vaginal") + theme(plot.title = element_text(hjust = 0.5, size= 18)),
                   ncol = 1, label_size = 20, labels = c('A', 'B'))
tax2 <- plot_grid(imp_rectal$plot + ggtitle("Rectal") + theme(plot.title = element_text(hjust = 0.5, size= 18)), 
                   imp_vaginal$plot + ggtitle("Vaginal") + theme(plot.title = element_text(hjust = 0.5, size= 18)),
                   ncol = 1, label_size = 20, labels = c('C', 'D'))

tax_final <- plot_grid(tax1, tax2, ncol= 2, label_size=20)
ggsave("results/aim4/07-final_taxa_only.png", tax_final,
       width = 15, height = 10, dpi = 300)



###################################################################################################
###################################################################################################
###################################################################################################
#### METADATA ONLY ################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
# =============================================================================
# Part 1: Load data
# =============================================================================

ps <- readRDS("data/phyloseq_filtered.rds") %>%
  tax_glom("Genus")

# =============================================================================
# Part 2: Helper function to prepare data for RF
# =============================================================================

# Metadata columns to include as predictors
meta_vars <- c("abnormal_uterine_bleeding", "Heavy_menstrual_bleeding",
               "Irregular_menstrual_bleeding", "Confirmed_adenomyosis",
               "Confirmed_cysts", "Confirmed_fibroids")

avg_abundance = taxa_sums(ps)/sum(taxa_sums(ps)) 
# Sort high to low
avg_abundance = sort(avg_abundance, decreasing = T)
# Take the top 10
top_10 = avg_abundance[1:10]
# Extract taxa names 
top_10 = names(top_10)

prepare_rf_data_meta <- function(ps_subset, outcome_col = "Host_disease") {
  df <- data.frame(sample_data(ps_subset)) %>%
    rownames_to_column("Sample") %>%
    select(Sample, all_of(outcome_col), all_of(meta_vars)) %>%
    mutate(across(all_of(meta_vars), as.numeric)) %>%
    na.omit() %>%
    select(-Sample)
  
  return(df)
}

# =============================================================================
# Part 3: Helper function to run RF with k-fold CV + hyperparameter tuning
# =============================================================================
run_rf <- function(df, outcome_col = "Host_disease", k = 5,
                   repeats = 10, seed = 421) {
  set.seed(seed)
  
  predictors <- df %>% select(-all_of(outcome_col))
  outcome <- df %>%
    pull(all_of(outcome_col)) %>%
    factor(levels = c("Control", "CPP Only", "CPP Endo"))
  
  levels(outcome) <- make.names(levels(outcome))
  
  tune_grid <- expand.grid(
    mtry = c(1, 3, 5, 6, 10),
    splitrule = c("gini", "extratrees"),
    min.node.size = c(1, 2, 3, 5, 10)
  ) %>%
    filter(mtry <= ncol(predictors))
  
  train_control <- trainControl(
    method = "repeatedcv",
    number = k,
    repeats = repeats,
    classProbs = TRUE,
    summaryFunction = multiClassSummary,
    savePredictions = "final"
  )
  
  rf_model <- train(
    x = predictors,
    y = outcome,
    method = "ranger",
    trControl = train_control,
    tuneGrid = tune_grid,
    metric = "AUC",
    importance = "impurity",
    num.trees = 1000
  )
  
  return(rf_model)
}

# =============================================================================
# Part 4: Subset by body site and prepare data
# =============================================================================

ps_rectal <- subset_samples(ps, env_medium == "rectal")
ps_vaginal <- subset_samples(ps, env_medium == "vaginal")

# Remove zero-abundance taxa after subsetting
ps_rectal <- prune_taxa(taxa_sums(ps_rectal) > 0, ps_rectal)
ps_vaginal <- prune_taxa(taxa_sums(ps_vaginal) > 0, ps_vaginal)

df_rectal <- prepare_rf_data_meta(ps_rectal)
df_vaginal <- prepare_rf_data_meta(ps_vaginal)

# =============================================================================
# Part 5: Run RF models
# =============================================================================

cat("=== Running RF for RECTAL samples ===\n")
rf_rectal_metadata <- run_rf(df_rectal, k = 5, seed = 123)
print(rf_rectal)

cat("\n=== Running RF for VAGINAL samples ===\n")
rf_vaginal_metadata <- run_rf(df_vaginal, k = 5, seed = 123)
print(rf_vaginal)

# =============================================================================
# Part 6: ROC curves (one-vs-rest for each class)
# =============================================================================

plot_multiclass_roc <- function(model, title) {
  preds <- model$pred
  
  best <- model$bestTune
  for (param in names(best)) {
    preds <- preds[preds[[param]] == best[[param]], ]
  }
  
  classes <- levels(preds$obs)
  roc_list <- list()
  auc_vals <- c()
  
  for (cls in classes) {
    binary_truth <- ifelse(preds$obs == cls, 1, 0)
    roc_obj <- pROC::roc(binary_truth, preds[[cls]], quiet = TRUE)
    roc_list[[cls]] <- roc_obj
    auc_vals[cls] <- as.numeric(pROC::auc(roc_obj))
  }
  
  roc_df <- bind_rows(lapply(classes, function(cls) {
    r <- roc_list[[cls]]
    data.frame(
      FPR = 1 - r$specificities,
      TPR = r$sensitivities,
      Class = cls
    )
  }))
  
  group_colors <- c(
    "Control" = "#c7e9b4",
    "CPP.Only" = "#41b6c4",
    "CPP.Endo" = "#225ea8"
  )
  
  roc_labels <- c(
    "Control" = paste0("Control (AUC=", round(auc_vals["Control"], 2), ")"),
    "CPP.Only" = paste0("CPP Only (AUC=", round(auc_vals["CPP.Only"], 2), ")"),
    "CPP.Endo" = paste0("CPP Endo (AUC=", round(auc_vals["CPP.Endo"], 2), ")")
  )
  
  print(unique(roc_df$Class))
  print(names(group_colors))
  
  p <- ggplot(roc_df, aes(x = FPR, y = TPR, colour = Class)) +
    geom_line(linewidth = 1) +
    geom_abline(
      slope = 1, intercept = 0,
      color = "gray", linetype = "dashed", linewidth = 0.8
    ) +
    scale_color_manual(
      values = group_colors,
      breaks = names(group_colors),
      labels = roc_labels
    ) +
    labs(
      x = "False Positive Rate",
      y = "True Positive Rate",
      title = title,
      color = "Host Disease"
    ) +
    theme_classic(base_size = 16) +
    theme(
      legend.position = "right",
      legend.direction = "vertical",
      axis.title = element_text(size=18),
    )
  
  return(list(plot = p, auc = auc_vals))
}

roc_rectal <- plot_multiclass_roc(rf_rectal, "")
roc_vaginal <- plot_multiclass_roc(rf_vaginal, "")

roc_rectal$plot
roc_vaginal$plot

ggsave("results/aim4/07-rf_roc_rectal_metadata_only.png", roc_rectal$plot,
       width = 8, height = 7, dpi = 300)
ggsave("results/aim4/07-rf_roc_vaginal_metadata_only.png", roc_vaginal$plot,
       width = 8, height = 7, dpi = 300)

# Print AUC values
cat("\n=== AUC values (Rectal) ===\n")
print(round(roc_rectal$auc, 3))
cat("\n=== AUC values (Vaginal) ===\n")
print(round(roc_vaginal$auc, 3))

# =============================================================================
# Part 7: Feature importance
# =============================================================================

plot_importance <- function(model, title) {
  imp <- varImp(model)$importance
  imp$Feature <- rownames(imp)
  
  # Use Overall importance (average across classes)
  imp <- imp %>%
    mutate(Importance = Overall) %>%
    arrange(desc(Importance)) %>%
    mutate(Feature = factor(Feature, levels = Feature)) %>%
    slice_max(Importance, n = 20)
  
  
  p <- ggplot(imp, aes(x = Feature, y = Importance, fill = Importance)) +
    geom_col() +
    theme_classic(base_size = 16) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          legend.position = "right",
          legend.direction = "vertical",
          axis.title = element_text(size=18),
          axis.text = element_text(size=14)) +
    ylab("Importance (Gini)") + xlab(NULL) +
    ggtitle(title)
  
  return(list(plot = p, data = imp))
}

imp_rectal <- plot_importance(rf_rectal, "")
imp_vaginal <- plot_importance(rf_vaginal, "")

imp_rectal$plot
imp_vaginal$plot

ggsave("results/aim4/07-rf_importance_rectal_metadata_only.png", imp_rectal$plot,
       width = 18, height = 10, dpi = 300)
ggsave("results/aim4/07-rf_importance_vaginal_metadata_only.png", imp_vaginal$plot,
       width = 18, height = 10, dpi = 300)

# =============================================================================
# Part 8: Summary table
# =============================================================================

results_table <- data.frame(
  Body_Site = c("Rectal", "Vaginal"),
  Best_mtry = c(rf_rectal$bestTune$mtry, rf_vaginal$bestTune$mtry),
  Best_splitrule = c(as.character(rf_rectal$bestTune$splitrule),
                     as.character(rf_vaginal$bestTune$splitrule)),
  Best_min_node = c(rf_rectal$bestTune$min.node.size,
                    rf_vaginal$bestTune$min.node.size),
  CV_Accuracy = c(max(rf_rectal$results$Accuracy),
                  max(rf_vaginal$results$Accuracy)),
  AUC_Control = c(roc_rectal$auc["Control"], roc_vaginal$auc["Control"]),
  AUC_CPP = c(roc_rectal$auc["CPP.Only"], roc_vaginal$auc["CPP.Only"]),
  AUC_CPP_Endo = c(roc_rectal$auc["CPP.Endo"], roc_vaginal$auc["CPP.Endo"])
)

print(results_table)
write_tsv(results_table, "results/aim4/07-rf_summary_metadata.tsv")


meta1 <- plot_grid(roc_rectal$plot + ggtitle("Rectal") + theme(plot.title = element_text(hjust = 0.5, size= 18),
                                                               plot.margin = ggplot2::margin(15, 15, 25, 15)), 
          roc_vaginal$plot + ggtitle("Vaginal") + theme(plot.title = element_text(hjust = 0.5, size= 18)),
          ncol = 1, label_size = 20, labels = c('A', 'B'))
meta2 <- plot_grid(imp_rectal$plot + ggtitle("Rectal") + theme(plot.title = element_text(hjust = 0.5, size= 18)), 
          imp_vaginal$plot + ggtitle("Vaginal") + theme(plot.title = element_text(hjust = 0.5, size= 18)),
          ncol = 1, label_size = 20, labels = c('C', 'D'))

meta_final <- plot_grid(meta1, meta2, ncol= 2, label_size=20)
ggsave("results/aim4/07-final_metadata_only.png", meta_final,
       width = 15, height = 10, dpi = 300)




