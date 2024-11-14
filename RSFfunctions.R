# Load user data
load_data <- function(file_path) {
  user_data <- read_excel(file_path)
  user_data$Race <- ifelse(user_data$Race == "Caucasian", 1, 0)
  user_data$Smoking <- ifelse(user_data$Smoking == "Current", 1, 0)
  user_data$Alcohol <- ifelse(user_data$Alcohol == "Yes", 1, 0)
  user_data$Drug <- ifelse(user_data$Drug == "Yes", 1, 0)
  user_data$HGB <- as.numeric(as.character(user_data$HGB))
  colnames(user_data)[(ncol(user_data)-5):(ncol(user_data)-3)] <- c("time", "event", "MAGGIC")
  return(user_data)
}

# Data preprocessing
preprocess_data <- function(data) {
  missing_percentage <- function(data) {
    colSums(is.na(data)) / nrow(data)
  }
  # There is some trickery here. The training set excludes features with 20% missing values, but the testing set is enforced to have the same features as the training set. So, change it to 0.5 for testing.
  data <- data[, missing_percentage(data) <= 0.2]
  data <- kNN(data, k = 5, imp_var = FALSE)
  
  clean_column_names_dot <- function(data) {
    colnames(data) <- gsub("[ /']", ".", colnames(data))
    return(data)
  }
  
  return(clean_column_names_dot(data))
}

# Create subsets of the data
create_subsets <- function(data) {
  D <- preprocess_data(cbind(data[, 4:70], data[, c("time", "event")]))
  M <- preprocess_data(cbind(data[, 71:150], data[, c("time", "event")]))
  DM <- preprocess_data(cbind(data[, 4:150], data[, c("time", "event")]))
  
  D <- D[, !duplicated(as.list(D))]
  M <- M[, !duplicated(as.list(M))]
  DM <- DM[, !duplicated(as.list(DM))]
  return(list(D = D, M = M, DM = DM))
}

# Feature selection function
rsf_vh_selection <- function(data, n_repeats = 25) {
  selected_vars_frequency <- setNames(rep(0, (ncol(data) - 2)), colnames(data)[-c(ncol(data)-1, ncol(data))])
  selected_vars_vimp <- setNames(rep(0, (ncol(data) - 2)), colnames(data)[-c(ncol(data)-1, ncol(data))])
  all_results <- list()
  
  for (rep in 1:n_repeats) {
    set.seed(rep)
    cat("Repeat:", rep, "\n")
    
    formula <- as.formula("Surv(time, event) ~ .")
    rsf <- rfsrc(formula, data = data, ntree = 100, importance = TRUE)
    
    var_select <- var.select(rsf)
    min_depth <- var_select$varselect
    vimp_scores <- min_depth[, "vimp"]
    names(vimp_scores) <- rownames(min_depth)
    threshold <- var_select$md.obj$threshold
    initial_vars <- rownames(min_depth[min_depth$depth <= threshold, ])
    all_selected_vars <- initial_vars
    oob_cindex <- rsf$err.rate[100]
    
    for (var in rownames(min_depth)[order(min_depth$depth)]) {
      if (!(var %in% initial_vars)) {
        all_selected_vars <- c(all_selected_vars, var)
        formula <- as.formula(paste("Surv(time, event) ~", paste(all_selected_vars, collapse = " + ")))
        rsf <- rfsrc(formula, data = data, ntree = 100, importance = TRUE)
        new_oob_cindex <- rsf$err.rate[100]
        if (new_oob_cindex >= oob_cindex) {
          all_selected_vars <- setdiff(all_selected_vars, var)
        } else {
          oob_cindex <- new_oob_cindex
        }
      }
    }
    
    all_results[[rep]] <- all_selected_vars
    for (var in all_selected_vars) {
      selected_vars_frequency[var] <- selected_vars_frequency[var] + 1
      selected_vars_vimp[var] <- selected_vars_vimp[var] + vimp_scores[var]
    }
  }
  selected_vars_vimp <- selected_vars_vimp / n_repeats
  selected_vars_sorted <- sort(selected_vars_vimp, decreasing = TRUE)
  
  return(list(all_results = all_results, selected_vars_sorted = selected_vars_sorted, selected_vars_frequency = selected_vars_frequency))
}

# Feature selection function with MAGGIC score integration
rsf_vh_selection_with_maggic <- function(user_data, M, n_repeats = 25) {
  if (!"MAGGIC" %in% colnames(user_data)) {
    stop("MAGGIC score is not present in the user_data.")
  }
  
  selected_vars_frequency <- setNames(rep(0, ncol(M) + 1), c("MAGGIC", colnames(M)))
  selected_vars_vimp <- setNames(rep(0, ncol(M) + 1), c("MAGGIC", colnames(M)))
  all_results <- list()
  
  for (rep in 1:n_repeats) {
    set.seed(rep)
    cat("Repeat:", rep, "\n")
    
    formula_full <- as.formula(paste("Surv(time, event) ~ MAGGIC +", paste(colnames(M), collapse = " + ")))
    rsf_full <- rfsrc(formula_full, data = user_data, ntree = 100, importance = TRUE)
    var_select <- var.select(rsf_full)
    min_depth <- var_select$varselect
    vimp_scores <- min_depth[, "vimp"]
    names(vimp_scores) <- rownames(min_depth)
    if ("MAGGIC" %in% rownames(min_depth)) {
      vimp_scores["MAGGIC"] <- rsf_full$importance["MAGGIC"]
    }
    
    all_selected_vars <- "MAGGIC"
    oob_cindex <- rsf_full$err.rate[100]
    
    for (var in rownames(min_depth)[order(min_depth$depth)]) {
      if (var != "MAGGIC") {
        all_selected_vars <- c(all_selected_vars, var)
        formula <- as.formula(paste("Surv(time, event) ~", paste(all_selected_vars, collapse = " + ")))
        rsf <- rfsrc(formula, data = user_data, ntree = 100, importance = TRUE)
        new_oob_cindex <- rsf$err.rate[100]
        if (new_oob_cindex >= oob_cindex) {
          all_selected_vars <- setdiff(all_selected_vars, var)
        } else {
          oob_cindex <- new_oob_cindex
        }
      }
    }
    
    all_results[[rep]] <- all_selected_vars
    for (var in all_selected_vars) {
      selected_vars_frequency[var] <- selected_vars_frequency[var] + 1
      selected_vars_vimp[var] <- selected_vars_vimp[var] + vimp_scores[var]
    }
  }
  
  selected_vars_vimp <- selected_vars_vimp / n_repeats
  selected_vars_sorted <- sort(selected_vars_vimp, decreasing = TRUE)
  
  return(list(all_results = all_results, selected_vars_sorted = selected_vars_sorted, selected_vars_frequency = selected_vars_frequency))
}

# New find_best_N_for_dataset function that first includes variables with more than 50% frequency, then gradually selects additional features
find_best_N_for_dataset <- function(reps, default_params, dataset, sorted_valid_features, mandatory_features, dataset_name) {
  # Initialize cindex_by_N to store evaluation results for each iteration
  cindex_by_N <- data.frame(N = integer(), CIndex = numeric(), Features = character(), stringsAsFactors = FALSE)
  selected_features <- mandatory_features # Start with mandatory features selected at more than 50% frequency
  remaining_features <- setdiff(names(sorted_valid_features), selected_features) # List of remaining features to be selected
  
  # Record the performance of the initial mandatory features
  for (i in 1:reps) {
    cat("Dataset:", dataset_name, "N =", length(mandatory_features), "Iteration:", i, "with features:", paste(selected_features, collapse = ","), "\n")
    results <- train_and_evaluate_rsf(dataset, dataset_name, default_params, selected_features, iteration = i)
    cindex_by_N <- rbind(cindex_by_N, data.frame(N = length(mandatory_features), CIndex = results$CIndex, Features = paste(selected_features, collapse = ","), stringsAsFactors = FALSE))
  }
  
  for (N in (length(mandatory_features) + 1):(length(mandatory_features) + 12)) {
    best_feature <- NULL
    best_cindex <- -Inf
    best_features_temp <- c()
    
    # Try each possible combination to identify the feature that maximizes the C-Index
    for (feature in remaining_features) {
      current_features <- c(selected_features, feature)
      
      # Record all C-Index values for the current feature selection process
      temp_cindices <- c()
      for (i in 1:reps) {
        cat("Dataset:", dataset_name, "N =", N, "Iteration:", i, "with feature:", feature, "\n")
        results <- train_and_evaluate_rsf(dataset, dataset_name, default_params, current_features, iteration = i)
        temp_cindices <- c(temp_cindices, results$CIndex)
      }
      avg_cindex <- mean(temp_cindices)
      
      # Update the best feature and best C-Index
      if (avg_cindex > best_cindex) {
        best_cindex <- avg_cindex
        best_feature <- feature
        best_features_temp <- current_features
      }
    }
    
    # Add the best feature to the selected features list, removing it from the remaining features list
    selected_features <- best_features_temp
    remaining_features <- setdiff(remaining_features, best_feature)
    
    # Record the best C-Index and feature names for the current N
    for (i in 1:reps) {
      cindex_by_N <- rbind(cindex_by_N, data.frame(N = N, CIndex = best_cindex, Features = paste(selected_features, collapse = ","), stringsAsFactors = FALSE))
    }
  }
  
  avg_cindex_by_N <- aggregate(CIndex ~ N, data = cindex_by_N, FUN = mean)
  best_N <- avg_cindex_by_N$N[which.max(avg_cindex_by_N$CIndex)]
  
  return(list(best_N = best_N, avg_cindex_by_N = avg_cindex_by_N, cindex_by_N = cindex_by_N))
}



optimize_rsf_params <- function(data, param_grid, selected_features) {
  
  best_oob <- Inf
  best_params <- list()
  param_results <- data.frame(
    ntree = numeric(), 
    mtry = numeric(), 
    nodesize = numeric(), 
    nodedepth = character(), 
    splitrule = character(), 
    samptype = character(),
    avg_oob = numeric(), 
    stringsAsFactors = FALSE
  )
  
  # Ensure selected_features includes time and event columns
  selected_features <- c(selected_features, "time", "event")
  
  # Subset the data to include only the selected features
  data <- data[, selected_features]
  
  total_combinations <- length(param_grid$ntree) * length(param_grid$mtry) * 
    length(param_grid$nodesize) * length(param_grid$nodedepth) * 
    length(param_grid$splitrule) * length(param_grid$samptype)
  progress <- 0
  
  # Iterate over all combinations in the parameter grid
  for (ntree in param_grid$ntree) {
    for (mtry in param_grid$mtry) {
      for (nodesize in param_grid$nodesize) {
        for (nodedepth in param_grid$nodedepth) {
          for (splitrule in param_grid$splitrule) {
            for (samptype in param_grid$samptype) {
              
              oob_scores <- numeric()
              
              for (i in 1:25) {
                # Handle "NoLimit" case for nodedepth explicitly
                if (nodedepth == "NoLimit") {
                  rf_model <- rfsrc(Surv(time, event) ~ ., data = data, ntree = ntree, mtry = mtry, nodesize = nodesize, 
                                    splitrule = splitrule, samptype = samptype)
                } else {
                  rf_model <- rfsrc(Surv(time, event) ~ ., data = data, ntree = ntree, mtry = mtry, nodesize = nodesize, 
                                    nodedepth = as.numeric(nodedepth), splitrule = splitrule, samptype = samptype)
                }
                
                # Record OOB error
                oob_scores <- c(oob_scores, rf_model$err.rate[length(rf_model$err.rate)])
              }
              
              avg_oob <- mean(oob_scores)
              param_results <- rbind(param_results, data.frame(
                ntree = ntree, mtry = mtry, nodesize = nodesize, 
                nodedepth = nodedepth, splitrule = splitrule, 
                samptype = samptype, avg_oob = avg_oob,
                stringsAsFactors = FALSE
              ))
              
              if (avg_oob < best_oob) {  # Select best parameters based on OOB error
                best_oob <- avg_oob
                best_params <- list(ntree = ntree, mtry = mtry, nodesize = nodesize, nodedepth = nodedepth, splitrule = splitrule, samptype = samptype)
              }
              
              # Update progress
              progress <- progress + 1
              cat("Progress: ", round(100 * progress / total_combinations, 2), "% - Avg OOB: ", avg_oob, "\n")
              
              # Plot OOB error evolution
              plot(param_results$avg_oob, type = "l", col = "blue", xlab = "Parameter Set", ylab = "Average OOB Error",
                   main = "OOB Error Progress")
              abline(h = best_oob, col = "red", lty = 2)
            }
          }
        }
      }
    }
  }
  
  # Final OOB error plot
  ggplot(param_results, aes(x = 1:nrow(param_results), y = avg_oob)) +
    geom_line(color = "blue") +
    geom_hline(yintercept = best_oob, color = "red", linetype = "dashed") +
    labs(title = "OOB Error Across Parameter Combinations", x = "Parameter Set", y = "Average OOB Error")
  
  return(list(best_params = best_params, param_results = param_results))
}

# Define function to extract feature importance
get_feature_importance <- function(model) {
  importance <- vimp(model)
  importance_data <- data.frame(Feature = names(importance$importance), Importance = importance$importance)
  importance_data
}

# Train and evaluate the RSF model
train_and_evaluate_rsf <- function(data, model_name, params, selected_features = NULL, predict_times = c(1), iteration) {
  # Remove rows with missing values
  data <- data[complete.cases(data), ]
  
  # If selected_features is provided, select the corresponding columns
  if (!is.null(selected_features)) {
    data <- data[, c("time", "event", selected_features)]
  }
  
  # Set up the rfsrc model parameters
  rfsrc_params <- list(
    formula = Surv(time, event) ~ .,
    data = data,
    ntree = params$ntree,
    mtry = params$mtry,
    nodesize = params$nodesize,
    splitrule = params$splitrule,
    importance = TRUE,
    forest = TRUE
  )
  
  # Handle the nodedepth parameter
  if (!is.null(params$nodedepth)) {
    if (params$nodedepth == "NoLimit") {
      rfsrc_params$nodedepth <- NULL
    } else {
      rfsrc_params$nodedepth <- as.numeric(params$nodedepth)
      if (is.na(rfsrc_params$nodedepth)) {
        stop("nodedepth must be numeric or 'NoLimit'")
      }
    }
  }
  
  # Train the model
  best_model <- do.call(rfsrc, rfsrc_params)
  
  # Get feature importance
  feature_importance <- get_feature_importance(best_model)
  
  # Calculate the out-of-bag C-index
  oob_cindex <- 1 - best_model$err.rate[length(best_model$err.rate)]
  
  # Use the survAUC package to calculate out-of-bag Brier Score and AUC
  times <- seq(0.25, 2.75, by = 0.25)
  Surv.rsf <- Surv(data$time, data$event)
  RSF.pred <- best_model$predicted.oob
  
  brier_scores <- predErr(Surv.rsf, Surv.rsf, RSF.pred, RSF.pred, times, type = "brier")
  auc_scores <- AUC.uno(Surv.rsf, Surv.rsf, RSF.pred, times)
  
  brier_data <- data.frame(Time = brier_scores$times, BrierScore = brier_scores$error, Model = model_name)
  auc_data <- data.frame(Time = auc_scores$times, AUC = auc_scores$auc, Model = model_name)
  iAUC_value <- auc_scores$iauc
  
  # Calculate ROC curve data for specific time points
  roc_data_list <- lapply(predict_times, function(predict_time) {
    roc_curve <- survivalROC(Stime = data$time, status = data$event, marker = RSF.pred, predict.time = predict_time, method = "KM")
    data.frame(FPR = roc_curve$FP, TPR = roc_curve$TP, Model = model_name, PredictTime = predict_time, Iteration = iteration)
  })
  
  roc_data <- do.call(rbind, roc_data_list)
  
  return(list(
    BrierScoreData = brier_data, 
    AUCData = auc_data, 
    iAUC = iAUC_value, 
    CIndex = oob_cindex, 
    FeatureImportance = feature_importance,
    ROCData = roc_data,
    model = best_model  # Include the model object in the returned list
  ))
}

# Calculate MAGGIC score metrics
calculate_MAGGIC_metrics <- function(data) {
  # Extract MAGGIC score
  risks <- data$MAGGIC
  
  # Calculate C-index
  formula <- as.formula("Surv(time, event) ~ risks")
  cindex <- 1 - concordance(Surv(data$time, data$event) ~ risks)$concordance
  
  # Time point sequence
  times <- seq(0.25, 4.75, by = 0.25)
  
  # Calculate AUC
  auc <- survAUC::AUC.uno(Surv(data$time, data$event), Surv(data$time, data$event), risks, times)$auc
  iauc <- survAUC::AUC.uno(Surv(data$time, data$event), Surv(data$time, data$event), risks, times)$iauc
  
  # Calculate Brier Score
  Surv.rsf <- Surv(data$time, data$event)
  RSF.pred <- risks
  brier_scores <- predErr(Surv.rsf, Surv.rsf, RSF.pred, RSF.pred, times, type = "brier")
  
  list(cindex = cindex, auc = auc, iauc = iauc, brier_scores = brier_scores$error, times = times)
}

# Save results to Excel
save_results_to_excel <- function(brier_results, auc_results, cindex_results, file_prefix) {
  # Use tidyr and purrr packages to group and process data
  
  # Save AUC results
  sheets_auc <- auc_results %>%
    group_by(Model) %>%
    nest() %>%
    mutate(sheet_name = paste0(file_prefix, "_", Model, "_AUC"))
  
  wb <- createWorkbook()
  pwalk(list(sheets_auc$sheet_name, sheets_auc$data), function(sheet_name, data) {
    df <- as_tibble(data)
    # Convert data from long format to wide format
    df <- df %>%
      group_by(Time) %>%
      mutate(measurement = row_number()) %>%
      pivot_wider(names_from = Time, values_from = AUC, values_fill = list(AUC = NA))
    
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet = sheet_name, df)
  })
  
  saveWorkbook(wb, paste0(file_prefix, "_auc_results.xlsx"), overwrite = TRUE)
  
  # Save Brier Score results
  sheets_brier <- brier_results %>%
    group_by(Model) %>%
    nest() %>%
    mutate(sheet_name = paste0(file_prefix, "_", Model, "_Brier"))
  
  wb <- createWorkbook()
  pwalk(list(sheets_brier$sheet_name, sheets_brier$data), function(sheet_name, data) {
    df <- as_tibble(data)
    # Convert data from long format to wide format
    df <- df %>%
      group_by(Time) %>%
      mutate(measurement = row_number()) %>%
      pivot_wider(names_from = Time, values_from = BrierScore, values_fill = list(BrierScore = NA))
    
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet = sheet_name, df)
  })
  
  saveWorkbook(wb, paste0(file_prefix, "_brier_results.xlsx"), overwrite = TRUE)
  
  # Save C-index results
  sheets_cindex <- cindex_results %>%
    group_by(Model) %>%
    nest() %>%
    mutate(sheet_name = paste0(file_prefix, "_", Model, "_Cindex"))
  
  wb <- createWorkbook()
  pwalk(list(sheets_cindex$sheet_name, sheets_cindex$data), function(sheet_name, data) {
    df <- as_tibble(data)
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet = sheet_name, df)
  })
  
  saveWorkbook(wb, paste0(file_prefix, "_cindex_results.xlsx"), overwrite = TRUE)
}

# Evaluate existing models
evaluate_existing_models <- function(rsf_model_list, predict_times, data, model_name) {
  cindex_results <- data.frame()
  brier_results <- data.frame()
  auc_results <- data.frame()
  iauc_values <- data.frame()
  roc_data <- data.frame()
  
  for (model_idx in 1:length(rsf_model_list)) {
    rsf_model <- rsf_model_list[[model_idx]]
    cat("Evaluating model index:", model_idx, "\n")
    
    # Calculate C-index
    cindex <- 1 - rsf_model$err.rate[length(rsf_model$err.rate)]
    cindex_results <- rbind(cindex_results, data.frame(Model = model_name, CIndex = cindex))
    
    # Calculate Brier Score and AUC
    times <- seq(0.25, 4.75, by = 0.25)
    Surv.rsf <- Surv(data$time, data$event)
    RSF.pred <- rsf_model$predicted.oob
    
    # Ensure RSF.pred is a matrix
    if (!is.matrix(RSF.pred)) {
      RSF.pred <- as.matrix(RSF.pred)
    }
    
    brier_scores <- predErr(Surv.rsf, Surv.rsf, RSF.pred, RSF.pred, times, type = "brier")
    auc_scores <- AUC.uno(Surv.rsf, Surv.rsf, RSF.pred, times)
    
    brier_data <- data.frame(Time = brier_scores$times, BrierScore = brier_scores$error, Model = model_name)
    auc_data <- data.frame(Time = auc_scores$times, AUC = auc_scores$auc, Model = model_name)
    iAUC_value <- auc_scores$iauc
    
    auc_results <- rbind(auc_results, auc_data)
    brier_results <- rbind(brier_results, brier_data)
    iauc_values <- rbind(iauc_values, data.frame(Model = model_name, iAUC = iAUC_value))
    
    # ROC Curve Data
    for (predict_time in predict_times) {
      roc_curve <- survivalROC(Stime = data$time, status = data$event, marker = RSF.pred, predict.time = predict_time, method = "KM")
      roc_data <- rbind(roc_data, data.frame(FPR = roc_curve$FP, TPR = roc_curve$TP, Model = model_name, PredictTime = predict_time))
    }
  }
  
  return(list(cindex = cindex_results, brier = brier_results, auc = auc_results, iauc = iauc_values, roc = roc_data))
}

evaluate_existing_models_external <- function(rsf_model_list, predict_times, validation_data, model_name) {
  cindex_results <- data.frame()
  brier_results <- data.frame()
  auc_results <- data.frame()
  iauc_values <- data.frame()
  roc_data <- data.frame()
  risks_all_models <- data.frame(PatientID = 1:nrow(validation_data))  # Used to store the risk scores from all models
  
  for (model_idx in 1:length(rsf_model_list)) {
    rsf_model <- rsf_model_list[[model_idx]]
    cat("Evaluating model index:", model_idx, "\n")
    validation_data$time <- as.numeric(validation_data$time)
    
    # Generate predictions using the external validation set
    RSF.pred <- predict(rsf_model, newdata = validation_data)$predicted
    
    # Ensure RSF.pred is a matrix or vector
    if (!is.matrix(RSF.pred)) {
      RSF.pred <- as.matrix(RSF.pred)
    }
    
    # Store risks for this model
    risks <- RSF.pred
    risks_all_models[paste0("Model_", model_idx)] <- risks
    
    # Compute C-index using the external validation data and the risks (RSF.pred)
    cindex <- 1 - concordance(Surv(validation_data$time, validation_data$event) ~ risks)$concordance
    cindex_results <- rbind(cindex_results, data.frame(Model = model_name, CIndex = cindex))
    
    # Compute Brier Score (without using pec)
    times <- seq(0.25, 1, by = 0.25)
    Surv.rsf <- Surv(validation_data$time, validation_data$event)
    brier_scores <- predErr(Surv.rsf, Surv.rsf, RSF.pred, RSF.pred, times, type = "brier")
    
    # Store Brier scores
    brier_data <- data.frame(Time = brier_scores$time, BrierScore = brier_scores$error, Model = model_name)
    brier_results <- rbind(brier_results, brier_data)
    
    # Compute AUC using Uno's estimator
    auc_scores <- AUC.uno(Surv(validation_data$time, validation_data$event), Surv(validation_data$time, validation_data$event), RSF.pred, times)
    auc_data <- data.frame(Time = auc_scores$times, AUC = auc_scores$auc, Model = model_name)
    iAUC_value <- auc_scores$iauc
    
    auc_results <- rbind(auc_results, auc_data)
    iauc_values <- rbind(iauc_values, data.frame(Model = model_name, iAUC = iAUC_value))
    
    # Compute ROC data
    for (time in predict_times) {
      roc_curve <- survivalROC(Stime = validation_data$time, status = validation_data$event, marker = RSF.pred, predict.time = time, method = "KM")
      roc_data <- rbind(roc_data, data.frame(Model = model_name, Time = time, FPR = roc_curve$FP, TPR = roc_curve$TP))
    }
  }
  
# Calculate the average risk for each patient across all models
risks_all_models$MeanRisk <- rowMeans(risks_all_models[, -1], na.rm = TRUE)  # Exclude PatientID column

return(list(
  cindex_results = cindex_results, 
  brier_results = brier_results, 
  auc_results = auc_results, 
  iauc_values = iauc_values, 
  roc_data = roc_data,
  average_risks = risks_all_models  # Add the average risk output
))
}

save_roc_to_excel <- function(roc_data, file_prefix) {
  wb <- createWorkbook()
  
  sheets_roc <- roc_data %>%
    group_by(Model, Time) %>%
    nest() %>%
    mutate(sheet_name = paste0(file_prefix, "_", Model, "_ROC_", Time))
  
  pwalk(list(sheets_roc$sheet_name, sheets_roc$data), function(sheet_name, data) {
    df <- as_tibble(data)
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet = sheet_name, df)
  })
  
  saveWorkbook(wb, paste0(file_prefix, "_roc_results.xlsx"), overwrite = TRUE)
}