# Clear all variables in the environment
rm(list = ls())
print(ls())

# Source the custom functions
source("C:/Users/fenggu/University of Michigan Dropbox/Feng Gu/GitHub/MachineLearningOnHFData/RSFfunctions.R")

# Load and install required packages
required_packages <- c("randomForestSRC", "survival", "Hmisc", "readxl", "ggplot2", "caret", "dplyr", "VIM", "openxlsx", "survAUC", "survivalROC", "tidyr", "purrr")

for (package in required_packages) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package)
  }
  library(package, character.only = TRUE)
}

file_path <- "C:/Users/fenggu/University of Michigan Dropbox/Feng Gu/GitHub/MachineLearningOnHFData/PredictorsFinal.xlsx"  # training
# file_path <- "C:/Users/fenggu/University of Michigan Dropbox/Feng Gu/TriSeg/AM_data_Fengedited.xlsx"  # Testing
user_data <- load_data(file_path)
subsets <- create_subsets(user_data)

default_params <- list(
  ntree = 1000,         # Number of trees to grow
  nodesize = 15,        # Minimum size of terminal nodes
  nodedepth = NULL,     # Maximum depth to which a tree should be grown (NULL means no limit)
  splitrule = "logrank", # Splitting rule for survival forests
  nsplit = 10,          # Number of random splits to consider for each candidate variable
  samptype = "swor",    # Sampling type: "swor" for without replacement, "swr" for with replacement
  sampsize = function(n) min(1, 0.632 * n) # Sample size for each tree (typically around 63.2% of the original data size)
)

# Hyperparameter optimization
param_grid <- list(
  ntree = c(1000, 1500, 2000),
  mtry = c(4, 5, 6, 7),
  nodesize = c(10, 30, 50),
  nodedepth = c(1, 3, 5, 10, 15, 20, "NoLimit"),
  splitrule = c("logrank"),
  samptype = "swor",
  sampsize = function(n) max(1, 0.632 * n)
)

predict_times <- c(0.5, 1, 2, 3, 5) # training
# predict_times <- c(0.5, 1) # testing

# Choose an option to perform a specific action
# mode <- 1  # Set the option mode: 1, 2, or 3
for (mode in 1) {

if (mode == 1) {
  # Option 1: Train three RSF models (D, M, DM), compare with MAGGIC score, output relevant evaluation metrics and visualizations, and save models and evaluation metrics
  
  # Feature selection
  result_D <- rsf_vh_selection(subsets$D)
  result_M <- rsf_vh_selection(subsets$M)
  result_DM <- rsf_vh_selection(subsets$DM)
  
  # Process dataset D
  mandatory_features_D <- names(result_D$selected_vars_frequency)[result_D$selected_vars_frequency >= 20]
  valid_features_D <- names(result_D$selected_vars_frequency)[result_D$selected_vars_frequency >= 5]
  sorted_valid_features_D <- result_D$selected_vars_sorted[valid_features_D]
  Nresult_D <- find_best_N_for_dataset(1, default_params, subsets$D, sorted_valid_features_D, mandatory_features_D, "D")
  
  # Process dataset M
  mandatory_features_M <- names(result_M$selected_vars_frequency)[result_M$selected_vars_frequency >= 20]
  valid_features_M <- names(result_M$selected_vars_frequency)[result_M$selected_vars_frequency >= 5]
  sorted_valid_features_M <- result_M$selected_vars_sorted[valid_features_M]
  Nresult_M <- find_best_N_for_dataset(1, default_params, subsets$M, sorted_valid_features_M, mandatory_features_M, "M")
  
  # Process dataset DM
  mandatory_features_DM <- names(result_DM$selected_vars_frequency)[result_DM$selected_vars_frequency >= 20]
  valid_features_DM <- names(result_DM$selected_vars_frequency)[result_DM$selected_vars_frequency >= 5]
  sorted_valid_features_DM <- result_DM$selected_vars_sorted[valid_features_DM]
  Nresult_DM <- find_best_N_for_dataset(1, default_params, subsets$DM, sorted_valid_features_DM, mandatory_features_DM, "DM")
  
  # Get the best N and corresponding features
  best_N_D <- Nresult_D$best_N
  best_N_M <- Nresult_M$best_N
  best_N_DM <- Nresult_DM$best_N
  
  # Extract the feature set corresponding to the best N from cindex_by_N
  selected_features_D <- unlist(strsplit(Nresult_D$cindex_by_N[Nresult_D$cindex_by_N$N == best_N_D, "Features"][1], ","))
  selected_features_M <- unlist(strsplit(Nresult_M$cindex_by_N[Nresult_M$cindex_by_N$N == best_N_M, "Features"][1], ","))
  selected_features_DM <- unlist(strsplit(Nresult_DM$cindex_by_N[Nresult_DM$cindex_by_N$N == best_N_DM, "Features"][1], ","))
  
  # Initialize a list to store the final retained features
  final_features <- selected_features_DM
  
  # Iterate over features ending with _D
  for (feature in selected_features_DM) {
    if (grepl("_D$", feature)) {
      base_name <- gsub("_D$", "", feature)
      corresponding_s <- paste0(base_name, "_S")
      
      # If the corresponding _S feature exists, remove the _S feature
      if (corresponding_s %in% final_features) {
        final_features <- setdiff(final_features, corresponding_s)
      }
    }
    
    # Iterate over features ending with _S
    if (grepl("_S$", feature)) {
      base_name <- gsub("_S$", "", feature)
      corresponding_d <- paste0(base_name, "_D")
      
      # If the corresponding _D feature exists, remove the _S feature
      if (corresponding_d %in% final_features) {
        final_features <- setdiff(final_features, feature)
      }
    }
  }
  
  selected_features_DM <- final_features
  # Conditionally remove "R_PA" from selected_features_M if "PVR_S" is present
  if ("PVR_S" %in% selected_features_DM) {
    selected_features_DM <- selected_features_DM[selected_features_DM != "R_PA"]
  }
  
  
  if ("PVR_S" %in% selected_features_M) {
    selected_features_M <- selected_features_M[selected_features_M != "R_PA"]
  }
  
  opt_D <- optimize_rsf_params(subsets$D, param_grid, selected_features_D)
  opt_M <- optimize_rsf_params(subsets$M, param_grid, selected_features_M)
  opt_DM <- optimize_rsf_params(subsets$DM, param_grid, selected_features_DM)
  best_params_D <- opt_D$best_params
  best_params_M <- opt_M$best_params
  best_params_DM <- opt_DM$best_params

  result_D_save_path <- "rsf_models_D.RData"
  result_M_save_path <- "rsf_models_M.RData"
  result_DM_save_path <- "rsf_models_DM.RData"
  D_model_results <- train_and_evaluate_models(subsets$D, selected_features_D, best_params_D, predict_times, "D", result_D_save_path)
  M_model_results <- train_and_evaluate_models(subsets$M, selected_features_M, best_params_M, predict_times, "M", result_M_save_path)
  DM_model_results <- train_and_evaluate_models(subsets$DM, selected_features_DM, best_params_DM, predict_times, "DM", result_DM_save_path)
 
  # Calculate and plot ROC curve for MAGGIC
  MAGGIC_ROC_plots <- calculate_and_plot_MAGGIC_ROC(user_data, predict_times)
  
  MAGGIC_metrics <- calculate_MAGGIC_metrics(user_data)
  MAGGIC_cindex <- MAGGIC_metrics$cindex
  MAGGIC_iauc <- data.frame(Model = "MAGGIC", iAUC = MAGGIC_metrics$iauc)
  
  MAGGIC_auc <- data.frame(Time = seq(0.25, 2.75, by = 0.25), AUC = unlist(MAGGIC_metrics$auc))
  MAGGIC_auc$Model <- "MAGGIC"
  
  MAGGIC_brier <- data.frame(Time = seq(0.25, 2.75, by = 0.25), BrierScore = unlist(MAGGIC_metrics$brier))
  MAGGIC_brier$Model <- "MAGGIC"
  
  # Add MAGGIC score to each of the results
  
  # rbind to process cindex_results
  cindex_results <- rbind(D_model_results$cindex_results, M_model_results$cindex_results, DM_model_results$cindex_results)
  cindex_results <- rbind(cindex_results, data.frame(Model = "MAGGIC", CIndex = MAGGIC_cindex))
  
  # Process auc_results
  auc_results <- rbind(D_model_results$auc_results, M_model_results$auc_results, DM_model_results$auc_results)
  auc_results <- rbind(auc_results, MAGGIC_auc)
  
  # Process brier_results
  brier_results <- rbind(D_model_results$brier_results, M_model_results$brier_results, DM_model_results$brier_results)
  brier_results <- rbind(brier_results, MAGGIC_brier)
  
  # Calculate and print summary results
  summary_results <- aggregate(CIndex ~ Model, data = cindex_results, function(x) c(mean = mean(x), sd = sd(x)))
  print(summary_results)
  
  integrated_auc_values <- rbind(
    data.frame(Model = "D", iAUC = D_model_results$iauc_values$iAUC),
    data.frame(Model = "M", iAUC = M_model_results$iauc_values$iAUC),
    data.frame(Model = "DM", iAUC = DM_model_results$iauc_values$iAUC),
    MAGGIC_iauc
  )
  
  integrated_auc_values <- integrated_auc_values %>%
    group_by(Model) %>%
    summarize(iAUC.mean = mean(iAUC), iAUC.sd = sd(iAUC)) %>%
    ungroup()
  
  print("Integrated AUC values:")
  print(integrated_auc_values)
  
  # Summarize Brier Score and AUC data and plot
  brier_summary <- brier_results %>%
    group_by(Model, Time) %>%
    summarize(Mean_BrierScore = mean(BrierScore), SD_BrierScore = sd(BrierScore), .groups = 'drop')
  
  auc_summary <- auc_results %>%
    group_by(Model, Time) %>%
    summarize(Mean_AUC = mean(AUC), SD_AUC = sd(AUC), .groups = 'drop')
  
  ggplot(brier_summary, aes(x = Time, y = Mean_BrierScore, color = Model)) +
    geom_line() +
    geom_point() +
    geom_ribbon(aes(ymin = Mean_BrierScore - SD_BrierScore, ymax = Mean_BrierScore + SD_BrierScore), alpha = 0.2, show.legend = FALSE) +
    labs(title = "Brier Score over Time for Different Models",
         x = "Time", y = "Brier Score") +
    theme_minimal() +
    scale_color_manual(values = c("D" = "blue", "M" = "green", "DM" = "red", "MAGGIC" = "purple"))
  
  ggplot(auc_summary, aes(x = Time, y = Mean_AUC, color = Model)) +
    geom_line() +
    geom_point() +
    geom_ribbon(aes(ymin = Mean_AUC - SD_AUC, ymax = Mean_AUC + SD_AUC), alpha = 0.2, show.legend = FALSE) +
    labs(title = "AUC over Time for Different Models",
         x = "Time", y = "AUC") +
    theme_minimal() +
    scale_color_manual(values = c("D" = "blue", "M" = "green", "DM" = "red", "MAGGIC" = "purple"))
  
  ggplot(cindex_results, aes(x = Model, y = CIndex)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.5) +
    labs(title = "C-index Distribution for Different Models",
         x = "Model", y = "C-index") +
    theme_minimal()
  
  
} else if (mode == 2) {
  # Option 2: Perform MAGGIC-based feature selection using MAGGIC score, train a model based on M and MAGGIC, output relevant evaluation metrics and visualizations
  
  result_M_with_maggic <- rsf_vh_selection_with_maggic(user_data, subsets$M)
  M_with_maggic <- cbind(subsets$M, user_data$MAGGIC)
  colnames(M_with_maggic)[ncol(M_with_maggic)] <- "MAGGIC"
  
  # Variables with a frequency of 13 or more are considered mandatory features
  mandatory_features <- names(result_M_with_maggic$selected_vars_frequency)[result_M_with_maggic$selected_vars_frequency >= 20]
  
  # Filter out variables with frequency less than 3
  valid_features <- names(result_M_with_maggic$selected_vars_frequency)[result_M_with_maggic$selected_vars_frequency >= 5]
  sorted_valid_features <- result_M_with_maggic$selected_vars_sorted[valid_features]
  
  # Call function
  Nresult_M_with_maggic <- find_best_N_for_dataset(
    reps = 3,
    default_params = default_params,
    dataset = M_with_maggic,
    sorted_valid_features = sorted_valid_features,
    mandatory_features = mandatory_features,
    dataset_name = "MAGGIC_with_M"
  )
  
  best_M_with_maggic <- Nresult_M_with_maggic$best_N
  selected_features_M_with_maggic <- unlist(strsplit(Nresult_M_with_maggic$cindex_by_N[Nresult_M_with_maggic$cindex_by_N$N ==  best_M_with_maggic, "Features"][1], ","))
  
  if ("PVR_S" %in% selected_features_M_with_maggic) {
    selected_features_M_with_maggic <- selected_features_M_with_maggic[selected_features_M_with_maggic != "R_PA"]
  }
  
  opt_M_with_maggic <- optimize_rsf_params(M_with_maggic, param_grid, selected_features_M_with_maggic)
  best_params_M_with_maggic <- opt_M_with_maggic$best_params
  
  result_M_with_maggic_save_path <- "rsf_models_M_with_maggic.RData"
  results_M_with_maggic <- train_and_evaluate_models(M_with_maggic, selected_features_M_with_maggic, best_params_M_with_maggic, predict_times, "MAGGIC_with_M",result_M_with_maggic_save_path)
  
  # Save results to Excel
  
  MAGGIC_metrics <- calculate_MAGGIC_metrics(user_data)
  MAGGIC_cindex <- MAGGIC_metrics$cindex
  MAGGIC_iauc <- data.frame(Model = "MAGGIC", iAUC = MAGGIC_metrics$iauc)
  
  MAGGIC_auc <- data.frame(Time = seq(0.25, 4.75, by = 0.25), AUC = unlist(MAGGIC_metrics$auc))
  MAGGIC_auc$Model <- "MAGGIC"
  
  MAGGIC_brier <- data.frame(Time = seq(0.25, 4.75, by = 0.25), BrierScore = unlist(MAGGIC_metrics$brier))
  MAGGIC_brier$Model <- "MAGGIC"
  
  cindex_results <- rbind(results_M_with_maggic$cindex_results, data.frame(Model = "MAGGIC", CIndex = MAGGIC_cindex))
  auc_results <- rbind(results_M_with_maggic$auc_results, MAGGIC_auc)
  brier_results <- rbind(results_M_with_maggic$ brier_results, MAGGIC_brier)
  
  # Calculate and print summary results
  summary_results <- aggregate(CIndex ~ Model, data = cindex_results, function(x) c(mean = mean(x), sd = sd(x)))
  print(summary_results)
  integrated_auc_values <- rbind(
    results_M_with_maggic$iauc_values,
    MAGGIC_iauc
  )
  
  integrated_auc_values <- integrated_auc_values %>%
    group_by(Model) %>%
    summarize(iAUC.mean = mean(iAUC), iAUC.sd = sd(iAUC)) %>%
    ungroup()
  
  print("Integrated AUC values:")
  print(integrated_auc_values)
  
  # Summarize Brier Score and AUC data and plot
  brier_summary <- brier_results %>%
    group_by(Model, Time) %>%
    summarize(Mean_BrierScore = mean(BrierScore), SD_BrierScore = sd(BrierScore), .groups = 'drop')
  
  auc_summary <- auc_results %>%
    group_by(Model, Time) %>%
    summarize(Mean_AUC = mean(AUC), SD_AUC = sd(AUC), .groups = 'drop')
  
  # Plot Brier Score over time
  ggplot(brier_summary, aes(x = Time, y = Mean_BrierScore, color = Model)) +
    geom_line() +
    geom_ribbon(aes(ymin = Mean_BrierScore -SD_BrierScore, ymax = Mean_BrierScore + SD_BrierScore), alpha = 0.2) +
    labs(title = "Brier Score over Time for MAGGIC vs M_with_maggic", x = "Time", y = "Brier Score") +
    theme_minimal() +
    scale_color_manual(values = c("MAGGIC" = "blue", "MAGGIC_with_M" = "green"))
  
  # Plot AUC over time
  ggplot(auc_summary, aes(x = Time, y = Mean_AUC, color = Model)) +
    geom_line() +
    geom_ribbon(aes(ymin = Mean_AUC - SD_AUC, ymax = Mean_AUC + SD_AUC), alpha = 0.2) +
    labs(title = "AUC over Time for MAGGIC vs M_with_maggic", x = "Time", y = "AUC") +
    theme_minimal() +
    scale_color_manual(values = c("MAGGIC" = "blue", "MAGGIC_with_M" = "green"))
  
  save_results_to_excel(results_M_with_maggic$brier_results, results_M_with_maggic$auc_results, results_M_with_maggic$cindex_results, "MAGGIC_M")
  results_M_with_maggic$roc_results <- results_M_with_maggic$roc_results %>%
    rename(Time = PredictTime)
  save_roc_to_excel(results_M_with_maggic$roc_results, "MAGGIC_M")
} else if (mode == 3) {  # # New Option 3: Load and evaluate existing RSF models, through external validation or internal validation set
  
  # Load and evaluate model D
  load("rsf_models_D.RData")
  #results_D <- evaluate_existing_models(rsf_models_D, predict_times, subsets$D, "D")
  results_D <- evaluate_existing_models_external(rsf_models_D, predict_times, subsets$D, "D")
  rm(rsf_models_D)  
  gc()
  
  # Load and evaluate model M
  load("rsf_models_M.RData")
  #results_M <- evaluate_existing_models(rsf_models_M, predict_times, subsets$M , "M")
  results_M <- evaluate_existing_models_external(rsf_models_M, predict_times, subsets$M, "M")
  rm(rsf_models_M)  # Release memory
  gc()
  
  # Load and evaluate model DM
  load("rsf_models_DM.RData")
  #results_DM <- evaluate_existing_models(rsf_models_DM, predict_times, subsets$DM, "DM")
  results_DM <- evaluate_existing_models_external(rsf_models_DM, predict_times, subsets$DM, "DM")
  rm(rsf_models_DM)  # Release memory
  gc()

  # Initialize result DataFrames
  cindex_results <- rbind(results_D$cindex_results, results_M$cindex_results, results_DM$cindex_results)
  brier_results <- rbind(results_D$brier_results, results_M$brier_results, results_DM$brier_results)
  auc_results <- rbind(results_D$auc_results, results_M$auc_results, results_DM$auc_results)
  iauc_values <- rbind(results_D$iauc_values, results_M$iauc_values, results_DM$iauc_values)
  
  # Calculate and plot the ROC curve for MAGGIC
  MAGGIC_ROC_plots <- calculate_and_plot_MAGGIC_ROC(user_data, predict_times)
  MAGGIC_roc_data <-  MAGGIC_ROC_plots$roc_data
  MAGGIC_metrics <- calculate_MAGGIC_metrics(user_data)
  MAGGIC_cindex <- MAGGIC_metrics$cindex
  MAGGIC_iauc <- data.frame(Model = "MAGGIC", iAUC = MAGGIC_metrics$iauc)
  
  MAGGIC_auc <- data.frame(Time = seq(0.25, 4.75, by = 0.25), AUC = unlist(MAGGIC_metrics$auc))
  MAGGIC_auc$Model <- "MAGGIC"
  
  MAGGIC_brier <- data.frame(Time = seq(0.25, 4.75, by = 0.25), BrierScore = unlist(MAGGIC_metrics$brier))
  MAGGIC_brier$Model <- "MAGGIC"
  
  # Add MAGGIC to each result
  cindex_results <- rbind(cindex_results, data.frame(Model = "MAGGIC", CIndex = MAGGIC_cindex))
  auc_results <- rbind(auc_results, MAGGIC_auc)
  brier_results <- rbind(brier_results, MAGGIC_brier)
  
  # Calculate and print summary results
  summary_results <- aggregate(CIndex ~ Model, data = cindex_results, function(x) c(mean = mean(x), sd = sd(x)))
  print(summary_results)
  
  integrated_auc_values <- rbind(
    iauc_values
    # MAGGIC_iauc
  )
  
  integrated_auc_values <- integrated_auc_values %>%
    group_by(Model) %>%
    summarize(iAUC.mean = mean(iAUC), iAUC.sd = sd(iAUC)) %>%
    ungroup()
  
  print("Integrated AUC values:")
  print(integrated_auc_values)
  
  # Summarize Brier Score and AUC data and plot
  brier_summary <- brier_results %>%
    group_by(Model, Time) %>%
    summarize(Mean_BrierScore = mean(BrierScore), SD_BrierScore = sd(BrierScore), .groups = 'drop')
  
  auc_summary <- auc_results %>%
    group_by(Model, Time) %>%
    summarize(Mean_AUC = mean(AUC), SD_AUC = sd(AUC), .groups = 'drop')
  
  ggplot(brier_summary, aes(x = Time, y = Mean_BrierScore, color = Model)) +
    geom_line() +
    geom_point() +
    geom_ribbon(aes(ymin = Mean_BrierScore - SD_BrierScore, ymax = Mean_BrierScore + SD_BrierScore), alpha = 0.2, show.legend = FALSE) +
    labs(title = "Brier Score over Time for Different Models",
         x = "Time", y = "Brier Score") +
    theme_minimal() +
    scale_color_manual(values = c("D" = "blue", "M" = "green", "DM" = "red", "MAGGIC" = "purple"))
  
  ggplot(auc_summary, aes(x = Time, y = Mean_AUC, color = Model)) +
    geom_line() +
    geom_point() +
    geom_ribbon(aes(ymin = Mean_AUC - SD_AUC, ymax = Mean_AUC + SD_AUC), alpha = 0.2, show.legend = FALSE) +
    labs(title = "AUC over Time for Different Models",
         x = "Time", y = "AUC") +
    theme_minimal() +
    scale_color_manual(values = c("D" = "blue", "M" = "green", "DM" = "red", "MAGGIC" = "purple"))
  
  ggplot(cindex_results, aes(x = Model, y = CIndex)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.5) +
    labs(title = "C-index Distribution for Different Models",
         x = "Model", y = "C-index") +
    theme_minimal()
  
  # Save results to Excel
  save_results_to_excel(results_D$brier_results, results_D$auc_results, results_D$cindex_results, "D")
  save_results_to_excel(results_M$brier_results, results_M$auc_results, results_M$cindex_results, "M")
  save_results_to_excel(results_DM$brier_results, results_DM$auc_results, results_DM$cindex_results, "DM")
  save_roc_to_excel(results_D$roc_data, "D")
  save_roc_to_excel(results_M$roc_data, "M")
  save_roc_to_excel(results_DM$roc_data, "DM")
  save_results_to_excel(MAGGIC_brier, MAGGIC_auc, data.frame(Model = "MAGGIC", CIndex = MAGGIC_cindex), "MAGGIC")
  save_roc_to_excel(MAGGIC_roc_data, "MAGGIC")
}
}

