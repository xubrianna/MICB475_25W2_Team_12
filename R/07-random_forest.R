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

set.seed(2026)

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

prepare_rf_data <- function(ps_subset, prev_threshold = 0.1) {
  # Prevalence filter: keep genera present in >= prev_threshold fraction of samples
  # This is less biased than top-N abundance selection (no data leakage)
  prev_cutoff <- ceiling(prev_threshold * nsamples(ps_subset))
  ps_prev <- prune_taxa(
    apply(otu_table(ps_subset), 1, function(x) sum(x > 0)) >= prev_cutoff,
    ps_subset
  )
  
  # CLR transform on all retained taxa
  ps_clr <- ps_prev %>% microbiome::transform("clr")
  
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
    select(Sample, Host_disease, all_of(meta_vars), feature_name, Abundance) %>%
    pivot_wider(names_from = feature_name, values_from = Abundance)
  
  # Ensure metadata columns are numeric
  df_pivot <- df_pivot %>%
    mutate(across(all_of(meta_vars), as.numeric))
  
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
    factor(levels = c("Control", "CPP", "CPP Endo"))
  
  # Hyperparameter grid
  tune_grid <- expand.grid(
    mtry = c(3, 6, 10),
    splitrule = c("gini", "extratrees"),
    min.node.size = c(2, 3, 4)
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
  # Filter to best tune using all tuning parameters generically
  best <- model$bestTune
  for (param in names(best)) {
    preds <- preds[preds[[param]] == best[[param]], ]
  }
  
  # Compute one-vs-rest ROC for each class
  classes <- levels(preds$obs)
  roc_list <- list()
  auc_vals <- c()
  
  for (cls in classes) {
    binary_truth <- ifelse(preds$obs == cls, 1, 0)
    prob_col <- cls
    roc_obj <- roc(binary_truth, preds[[prob_col]], quiet = TRUE)
    roc_list[[cls]] <- roc_obj
    auc_vals[cls] <- auc(roc_obj)
  }
  
  # Build data frame for plotting
  roc_df <- bind_rows(lapply(classes, function(cls) {
    r <- roc_list[[cls]]
    data.frame(
      FPR = 1 - r$specificities,
      TPR = r$sensitivities,
      Class = paste0(gsub("\\.", " ", cls), 
                     " (AUC=", round(auc_vals[cls], 2), ")")
    )
  }))
  
  p <- ggplot(roc_df, aes(x = FPR, y = TPR, color = Class)) +
    geom_line(size = 1) +
    geom_abline(slope = 1, intercept = 0, color = "gray", 
                linetype = "dashed", size = 0.8) +
    labs(x = "False Positive Rate", y = "True Positive Rate",
         title = title, color = "Class") +
    theme_minimal(base_size = 16) +
    theme(legend.position = "bottom",
          legend.direction = "vertical")
  
  return(list(plot = p, auc = auc_vals))
}

roc_rectal <- plot_multiclass_roc(rf_rectal, "ROC Curves - Rectal Samples")
roc_vaginal <- plot_multiclass_roc(rf_vaginal, "ROC Curves - Vaginal Samples")

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
  
  # Use Overall importance (average across classes)
  imp <- imp %>%
    mutate(Importance = Overall) %>%
    mutate(Feature = gsub("g__", "", Feature)) %>%
    arrange(desc(Importance)) %>%
    mutate(Feature = factor(Feature, levels = Feature)) %>%
    slice_max(Importance, n = 20)
    
  
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
  AUC_CPP = c(roc_rectal$auc["CPP"], roc_vaginal$auc["CPP"]),
  AUC_CPP_Endo = c(roc_rectal$auc["CPP.Endo"], roc_vaginal$auc["CPP.Endo"])
)

print(results_table)
write_tsv(results_table, "results/aim4/07-rf_summary.tsv")

# =============================================================================
# Part 9: Binary RF — CPP Endo vs CPP
# =============================================================================

# Helper to prepare binary data (filter to 2 groups only)
prepare_rf_data_binary <- function(ps_subset, group1, group2, prev_threshold = 0.1) {
  samp <- data.frame(sample_data(ps_subset))
  keep <- rownames(samp)[samp$Host_disease %in% c(group1, group2)]
  ps_sub <- prune_samples(keep, ps_subset)
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  
  # Prevalence filter instead of top-N (avoids data leakage)
  prev_cutoff <- ceiling(prev_threshold * nsamples(ps_sub))
  ps_prev <- prune_taxa(
    apply(otu_table(ps_sub), 1, function(x) sum(x > 0)) >= prev_cutoff,
    ps_sub
  )
  
  ps_clr <- ps_prev %>% microbiome::transform("clr")
  
  df <- psmelt(ps_clr)
  
  # Build a lookup: one clean name per OTU
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
  
  df_pivot <- df %>%
    select(Sample, Host_disease, all_of(meta_vars), feature_name, Abundance) %>%
    pivot_wider(names_from = feature_name, values_from = Abundance)
  
  # Ensure metadata columns are numeric
  df_pivot <- df_pivot %>%
    mutate(across(all_of(meta_vars), as.numeric))
  
  df_final <- df_pivot %>% na.omit() %>% select(-Sample)
  return(df_final)
}

# Helper to run binary RF
run_rf_binary <- function(df, ref_level, pos_level,
                          outcome_col = "Host_disease", k = 5,
                          repeats = 10, seed = 421) {
  set.seed(seed)
  
  predictors <- df %>% select(-all_of(outcome_col))
  outcome <- df %>% pull(all_of(outcome_col)) %>%
    factor(levels = c(ref_level, pos_level))
  
  # Make valid R names for caret
  levels(outcome) <- make.names(levels(outcome))
  
  tune_grid <- expand.grid(
    mtry = c(3, 6, 10),
    splitrule = c("gini", "extratrees"),
    min.node.size = c(2, 3, 4)
  )
  tune_grid <- tune_grid %>% filter(mtry <= ncol(predictors))
  
  # Repeated stratified k-fold CV
  train_control <- trainControl(
    method = "repeatedcv",
    number = k,
    repeats = repeats,
    classProbs = TRUE,
    summaryFunction = twoClassSummary,
    savePredictions = "final"
  )
  
  rf_model <- train(
    x = predictors,
    y = outcome,
    method = "ranger",
    trControl = train_control,
    tuneGrid = tune_grid,
    metric = "ROC",
    importance = "impurity",
    num.trees = 1000
  )
  
  return(rf_model)
}

# Helper to plot binary ROC
plot_binary_roc <- function(model, title) {
  preds <- model$pred
  best <- model$bestTune
  for (param in names(best)) {
    preds <- preds[preds[[param]] == best[[param]], ]
  }
  
  pos_class <- levels(preds$obs)[2]
  roc_obj <- roc(preds$obs, preds[[pos_class]],
                 levels = rev(levels(preds$obs)), quiet = TRUE)
  auc_val <- auc(roc_obj)
  ci_val <- ci.auc(roc_obj)
  
  roc_df <- data.frame(
    FPR = 1 - roc_obj$specificities,
    TPR = roc_obj$sensitivities
  )
  
  label <- sprintf("AUC = %.2f (%.2f-%.2f)",
                   auc_val, ci_val[1], ci_val[3])
  
  p <- ggplot(roc_df, aes(x = FPR, y = TPR)) +
    geom_line(size = 1, color = "black") +
    geom_abline(slope = 1, intercept = 0, color = "gray",
                linetype = "dashed", size = 0.8) +
    annotate("text", x = 0.65, y = 0.15, label = label, size = 5) +
    labs(x = "False Positive Rate", y = "True Positive Rate",
         title = title) +
    theme_minimal(base_size = 16)
  
  return(list(plot = p, auc = auc_val, ci = ci_val))
}

# --- Prepare binary data (CPP Endo vs CPP) ---
df_bin_rectal <- prepare_rf_data_binary(
  subset_samples(ps, env_medium == "rectal"),
  "CPP", "CPP Endo"
)

df_bin_vaginal <- prepare_rf_data_binary(
  subset_samples(ps, env_medium == "vaginal"),
  "CPP", "CPP Endo"
)

# --- Run binary RF ---
cat("\n=== Binary RF: CPP Endo vs CPP (Rectal) ===\n")
rf_bin_rectal <- run_rf_binary(df_bin_rectal, ref_level = "CPP",
                               pos_level = "CPP Endo", k = 5, seed = 421)
print(rf_bin_rectal)

cat("\n=== Binary RF: CPP Endo vs CPP (Vaginal) ===\n")
rf_bin_vaginal <- run_rf_binary(df_bin_vaginal, ref_level = "CPP",
                                pos_level = "CPP Endo", k = 5, seed = 421)
print(rf_bin_vaginal)

# --- ROC curves ---
roc_bin_rectal <- plot_binary_roc(rf_bin_rectal,
                                  "ROC: CPP Endo vs CPP - Rectal")
roc_bin_vaginal <- plot_binary_roc(rf_bin_vaginal,
                                   "ROC: CPP Endo vs CPP - Vaginal")

roc_bin_rectal$plot
roc_bin_vaginal$plot

ggsave("results/aim4/07-rf_roc_binary_rectal.png", roc_bin_rectal$plot,
       width = 8, height = 7, dpi = 300)
ggsave("results/aim4/07-rf_roc_binary_vaginal.png", roc_bin_vaginal$plot,
       width = 8, height = 7, dpi = 300)

cat("\n=== Binary AUC (Rectal): CPP Endo vs CPP ===\n")
cat(sprintf("AUC = %.3f (%.3f-%.3f)\n",
            roc_bin_rectal$auc, roc_bin_rectal$ci[1], roc_bin_rectal$ci[3]))
cat("\n=== Binary AUC (Vaginal): CPP Endo vs CPP ===\n")
cat(sprintf("AUC = %.3f (%.3f-%.3f)\n",
            roc_bin_vaginal$auc, roc_bin_vaginal$ci[1], roc_bin_vaginal$ci[3]))

# --- Feature importance (binary) ---
imp_bin_rectal <- plot_importance(rf_bin_rectal,
                                  "Feature Importance: CPP Endo vs CPP - Rectal")
imp_bin_vaginal <- plot_importance(rf_bin_vaginal,
                                   "Feature Importance: CPP Endo vs CPP - Vaginal")

imp_bin_rectal$plot
imp_bin_vaginal$plot

ggsave("results/aim4/07-rf_importance_binary_rectal.png", imp_bin_rectal$plot,
       width = 18, height = 10, dpi = 300)
ggsave("results/aim4/07-rf_importance_binary_vaginal.png", imp_bin_vaginal$plot,
       width = 18, height = 10, dpi = 300)

# =============================================================================
# Part 11: Model comparison summary
# =============================================================================

comparison_table <- data.frame(
  Model = c("RF", "RF", "SVM", "SVM", "GLMNet", "GLMNet"),
  Body_Site = rep(c("Rectal", "Vaginal"), 3),
  AUC_3class_Control = c(
    roc_rectal$auc["Control"], roc_vaginal$auc["Control"],
    roc_svm_rectal$auc["Control"], roc_svm_vaginal$auc["Control"],
    roc_glm_rectal$auc["Control"], roc_glm_vaginal$auc["Control"]
  ),
  AUC_3class_CPP = c(
    roc_rectal$auc["CPP"], roc_vaginal$auc["CPP"],
    roc_svm_rectal$auc["CPP"], roc_svm_vaginal$auc["CPP"],
    roc_glm_rectal$auc["CPP"], roc_glm_vaginal$auc["CPP"]
  ),
  AUC_3class_CPP_Endo = c(
    roc_rectal$auc["CPP.Endo"], roc_vaginal$auc["CPP.Endo"],
    roc_svm_rectal$auc["CPP.Endo"], roc_svm_vaginal$auc["CPP.Endo"],
    roc_glm_rectal$auc["CPP.Endo"], roc_glm_vaginal$auc["CPP.Endo"]
  ),
  AUC_binary = c(
    roc_bin_rectal$auc, roc_bin_vaginal$auc,
    roc_svm_bin_r$auc, roc_svm_bin_v$auc,
    roc_glm_bin_r$auc, roc_glm_bin_v$auc
  )
)

print(comparison_table)
write_tsv(comparison_table, "results/aim4/07-model_comparison.tsv")




