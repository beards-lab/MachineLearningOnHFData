# Clear the environment variables
rm(list = ls())

# Ensure the environment is empty
print(ls())

# Install and load necessary packages
required_packages <- c("glmnet", "survival", "readxl", "openxlsx", "dplyr", "pec", "survAUC", "ggplot2", "caret", "VIM", "pROC")
for (package in required_packages) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package)
  }
}

library(ggplot2)
library(glmnet)
library(survival)
library(readxl)
library(openxlsx)
library(pec)
library(survAUC)
library(dplyr)
library(caret)
library(VIM)
library(pROC)

# Load user data
file_path <- "C:/Users/fenggu/University of Michigan Dropbox/Feng Gu/GitHub/MachineLearningOnHFData/PredictorsFinal.xlsx"  
user_data <- read_excel(file_path)

user_data$Race <- ifelse(user_data$Race == "Caucasian", 1, 0)
user_data$Smoking <- ifelse(user_data$Smoking == "Current", 1, 0)
user_data$Alcohol <- ifelse(user_data$Alcohol == "Yes", 1, 0)
user_data$Drug <- ifelse(user_data$Drug == "Yes", 1, 0)
user_data$HGB <- as.numeric(as.character(user_data$HGB))

# View data
head(user_data)

# Rename columns for time and event
colnames(user_data)[(ncol(user_data)-5):(ncol(user_data)-3)] <- c("time", "event", "MAGGIC")

# Create subsets D, M, DM
D <- cbind(user_data[, 4:70], user_data[, c("time", "event")])
M <- cbind(user_data[, c(71:150)], user_data[, c("time", "event")])
DM <- cbind(user_data[, c(4:150)], user_data[, c("time", "event")])

# Remove duplicate columns and retain one
D <- D[, !duplicated(as.list(D))]
M <- M[, !duplicated(as.list(M))]
DM <- DM[, !duplicated(as.list(DM))]

# Calculate missing values percentage
missing_percentage <- function(data) {
  colSums(is.na(data)) / nrow(data)
}

# Drop columns with more than 20% missing values
D <- D[, missing_percentage(D) <= 0.2]
M <- M[, missing_percentage(M) <= 0.2]
DM <- DM[, missing_percentage(DM) <= 0.2]

# Fill missing values (less than 20%) using KNN
D <- kNN(D, k = 5, imp_var = FALSE)
M <- kNN(M, k = 5, imp_var = FALSE)
DM <- kNN(DM, k = 5, imp_var = FALSE)

# Clean column names (replace special characters with '.')
clean_column_names_dot <- function(data) {
  colnames(data) <- gsub("[ /']", ".", colnames(data))
  return(data)
}

# Apply column name cleaning to D, M, DM
D <- clean_column_names_dot(D)
M <- clean_column_names_dot(M)
DM <- clean_column_names_dot(DM)

# Univariate Cox regression to filter significant variables
univariate_cox <- function(data) {
  significant_vars <- c()
  
  for (var in colnames(data)[-c(ncol(data)-1, ncol(data))]) {
    formula <- as.formula(paste("Surv(time, event) ~", var))
    
    # Handle variables causing non-convergence
    cox_model <- tryCatch({
      coxph(formula, data = data)
    }, warning = function(w) {
      message("Warning for variable ", var, ": ", w$message)
      return(NULL)
    }, error = function(e) {
      message("Error for variable ", var, ": ", e$message)
      return(NULL)
    })
    
    if (!is.null(cox_model)) {
      p_value <- summary(cox_model)$coefficients[,"Pr(>|z|)"]
      
      # Ensure p_value is a single value
      if (length(p_value) == 1 && p_value < 0.05) {
        significant_vars <- c(significant_vars, var)
      }
    }
  }
  
  return(significant_vars)
}

# Apply univariate Cox regression on D, M, DM
sig_vars_D <- univariate_cox(D)
sig_vars_M <- univariate_cox(M)
sig_vars_DM <- univariate_cox(DM)

# Lasso-Cox regression with 10-fold cross-validation
lasso_cox_cv <- function(data, significant_vars) {
  # Extract and check data matrix
  x <- data[, significant_vars]
  
  # Ensure all variables are numeric, and handle NA/NaN/Inf
  x <- lapply(x, function(col) as.numeric(as.character(col)))
  x <- as.matrix(do.call(cbind, x))
  
  # Check for NA/NaN/Inf and remove rows containing these values
  if (any(is.na(x) | is.infinite(x))) {
    warning("Removing rows with NA/NaN/Inf values.")
    x <- x[complete.cases(x), ]
  }
  
  # Create Surv object
  y <- Surv(data$time, data$event)
  
  # Perform Lasso-Cox regression with 10-fold cross-validation
  cvfit <- cv.glmnet(x, y, family = "cox", alpha = 1, nfolds = 10)
  best_lambda <- cvfit$lambda.min
  final_model <- glmnet(x, y, family = "cox", alpha = 1, lambda = best_lambda)
  
  # Extract selected variables
  selected_vars <- rownames(coef(final_model))[which(coef(final_model) != 0)]
  
  return(list(selected_vars = selected_vars, cvfit = cvfit))
}

# Apply Lasso-Cox regression and cross-validation on D, M, DM
lasso_D <- lasso_cox_cv(D, sig_vars_D)
lasso_M <- lasso_cox_cv(M, sig_vars_M)
lasso_DM <- lasso_cox_cv(DM, sig_vars_DM)

# Calculate C-index, Brier Score, and AUC
calculate_metrics <- function(data, selected_vars) {
  # Extract and build the data matrix
  x <- data[, selected_vars]
  x <- as.matrix(do.call(cbind, lapply(x, function(col) as.numeric(as.character(col)))))
  
  # Create a Surv object
  y <- Surv(data$time, data$event)
  
  # Fit a Cox regression model with x=TRUE to keep the design matrix
  cox_model <- coxph(Surv(time, event) ~ ., data = data[, c(selected_vars, "time", "event")], x = TRUE)
  
  # Calculate C-index
  cindex <- summary(cox_model)$concordance[1]
  
  # Calculate Brier Score
  lp <- predict(cox_model)  # Linear predictors
  times <- seq(0.25, 4.75, by = 0.25)  # Specify time points
  brier_scores <- predErr(Surv(data$time, data$event), Surv(data$time, data$event), lp, lp, times = times, type = "brier")
  
  # Calculate AUC
  auc <- AUC.uno(Surv.rsp = y, Surv.rsp.new = y, lpnew = lp, times = times)
  
  # Return C-index, Brier Score, and AUC
  return(list(model = cox_model, cindex = cindex, brier = brier_scores$error, auc = auc$auc, iauc = auc$iauc))
}

# Apply the selected variables to D, M, DM datasets, build Cox models, and calculate metrics
metrics_D <- calculate_metrics(D, lasso_D$selected_vars)
metrics_M <- calculate_metrics(M, lasso_M$selected_vars)
metrics_DM <- calculate_metrics(DM, lasso_DM$selected_vars)

# Load external validation dataset
file_path_validation <- "C:/Users/fenggu/University of Michigan Dropbox/Feng Gu/TriSeg/AM_data_Fengedited.xlsx" 
validation_data <- read_excel(file_path_validation)

# Preprocess categorical and numeric variables in the validation set to match the training set
validation_data$Race <- ifelse(validation_data$Race == "Caucasian", 1, 0)
validation_data$Smoking <- ifelse(validation_data$Smoking == "Current", 1, 0)
validation_data$Alcohol <- ifelse(validation_data$Alcohol == "Yes", 1, 0)
validation_data$Drug <- ifelse(validation_data$Drug == "Yes", 1, 0)
validation_data$HGB <- as.numeric(as.character(validation_data$HGB))

# Rename time and event columns for consistency
colnames(validation_data)[(ncol(validation_data)-5):(ncol(validation_data)-3)] <- c("time", "event", "MAGGIC")

# Create subsets for D, M, and DM
D_val <- cbind(validation_data[, 4:70], validation_data[, c("time", "event")])
M_val <- cbind(validation_data[, c(71:150)], validation_data[, c("time", "event")])
DM_val <- cbind(validation_data[, c(4:150)], validation_data[, c("time", "event")])

# Remove columns with all missing values
D_val <- D_val[, colSums(is.na(D_val)) < nrow(D_val)]
M_val <- M_val[, colSums(is.na(M_val)) < nrow(M_val)]
DM_val <- DM_val[, colSums(is.na(DM_val)) < nrow(DM_val)]

# Impute missing values using KNN
D_val <- kNN(D_val, k = 5, imp_var = FALSE)
M_val <- kNN(M_val, k = 5, imp_var = FALSE)
DM_val <- kNN(DM_val, k = 5, imp_var = FALSE)

# Function to clean column names (replace with dots)
clean_column_names_dot <- function(data) {
  colnames(data) <- gsub("[ /']", ".", colnames(data))
  return(data)
}

# Apply column name cleaning to D, M, and DM
D_val <- clean_column_names_dot(D_val)
M_val <- clean_column_names_dot(M_val)
DM_val <- clean_column_names_dot(DM_val)

# Predict and evaluate using the trained Cox models on new datasets
predict_and_evaluate <- function(cox_model, data, selected_vars) {
  # Predict using the trained model
  predictors <- predict(cox_model, newdata = data[, selected_vars], type = "risk")
  
  # Create a survival object
  surv_obj <- Surv(data$time, data$event)
  
  # Calculate C-index
  cindex <- 1 - concordance(surv_obj ~ predictors)$concordance
  
  # Calculate Brier Score
  times <- seq(0.25, 4.75, by = 0.25)  # Specify time points
  brier_scores <- predErr(Surv(data$time, data$event), Surv(data$time, data$event), predictors, predictors, times = times, type = "brier")
  
  # Calculate AUC
  auc <- AUC.uno(Surv.rsp = surv_obj, Surv.rsp.new = surv_obj, lpnew = predictors, times = times)
  
  # Return C-index, Brier Score, AUC, and iAUC
  return(list(cindex = cindex, brier = brier_scores$error, auc = auc$auc, iauc = auc$iauc))
}

# Evaluate on the validation dataset
cindex_D_val <- predict_and_evaluate(metrics_D$model, D_val, lasso_D$selected_vars)
cindex_M_val <- predict_and_evaluate(metrics_M$model, M_val, lasso_M$selected_vars)
cindex_DM_val <- predict_and_evaluate(metrics_DM$model, DM_val, lasso_DM$selected_vars)

# Print the C-index results
cat("Validation D model C-index:", cindex_D_val$cindex, "\n")
cat("Validation M model C-index:", cindex_M_val$cindex, "\n")
cat("Validation DM model C-index:", cindex_DM_val$cindex, "\n")

library(openxlsx)

# Function to write results to an Excel file
write_results_to_excel <- function(file_name, metrics_D, metrics_M, metrics_DM, cindex_D_val, cindex_M_val, cindex_DM_val) {
  
  # Create a new workbook
  wb <- createWorkbook()
  
  # Create worksheets
  addWorksheet(wb, "Metrics_D")
  addWorksheet(wb, "Metrics_M")
  addWorksheet(wb, "Metrics_DM")
  addWorksheet(wb, "Cindex_Validation")
  
  # Write Metrics_D results
  writeData(wb, "Metrics_D", data.frame(Cindex = metrics_D$cindex, 
                                        Brier_Score = metrics_D$brier, 
                                        AUC = metrics_D$auc, 
                                        iAUC = metrics_D$iauc))
  
  # Write Metrics_M results
  writeData(wb, "Metrics_M", data.frame(Cindex = metrics_M$cindex, 
                                        Brier_Score = metrics_M$brier, 
                                        AUC = metrics_M$auc, 
                                        iAUC = metrics_M$iauc))
  
  # Write Metrics_DM results
  writeData(wb, "Metrics_DM", data.frame(Cindex = metrics_DM$cindex, 
                                         Brier_Score = metrics_DM$brier, 
                                         AUC = metrics_DM$auc, 
                                         iAUC = metrics_DM$iauc))
  
  # Save the workbook
  saveWorkbook(wb, file_name, overwrite = TRUE)
}

# Compute and obtain metrics
metrics_D <- calculate_metrics(D, lasso_D$selected_vars)
metrics_M <- calculate_metrics(M, lasso_M$selected_vars)
metrics_DM <- calculate_metrics(DM, lasso_DM$selected_vars)

# Validation set C-index
metrics_D  <- predict_and_evaluate(metrics_D$model, D_val, lasso_D$selected_vars)
metrics_M <- predict_and_evaluate(metrics_M$model, M_val, lasso_M$selected_vars)
metrics_DM <- predict_and_evaluate(metrics_DM$model, DM_val, lasso_DM$selected_vars)

# Write results to an Excel file
write_results_to_excel("model_evaluation_results.xlsx", metrics_D, metrics_M, metrics_DM)