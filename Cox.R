# 清除环境变量
rm(list = ls())

# 确认环境中已经没有变量
print(ls())

# 安装并加载必要的包
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

# 加载用户数据
file_path <- "C:/Users/fenggu/University of Michigan Dropbox/Feng Gu/GitHub/MachineLearningOnHFData/PredictorsFinal.xlsx"  
user_data <- read_excel(file_path)

# 转换列
# user_data$HFtype <- ifelse(user_data$HFtype == "HFpEF", 1, 0)
user_data$Race <- ifelse(user_data$Race == "Caucasian", 1, 0)
user_data$Smoking <- ifelse(user_data$Smoking == "Current", 1, 0)
user_data$Alcohol <- ifelse(user_data$Alcohol == "Yes", 1, 0)
user_data$Drug <- ifelse(user_data$Drug == "Yes", 1, 0)
user_data$HGB <- as.numeric(as.character(user_data$HGB))

# 查看数据
head(user_data)



# 重命名时间和事件列
colnames(user_data)[(ncol(user_data)-5):(ncol(user_data)-3)] <- c("time", "event", "MAGGIC")

# 创建D, M, DM子集
D <- cbind(user_data[, 4:70], user_data[, c("time", "event")])
M <- cbind(user_data[, c(71:150)], user_data[, c("time", "event")])
DM <- cbind(user_data[, c(4:150)], user_data[, c("time", "event")])

# 删除重复的列，保留一个
D <- D[, !duplicated(as.list(D))]
M <- M[, !duplicated(as.list(M))]
DM <- DM[, !duplicated(as.list(DM))]

# 计算缺失值比例
missing_percentage <- function(data) {
  colSums(is.na(data)) / nrow(data)
}

# 删除缺失值超过20%的列
D <- D[, missing_percentage(D) <= 0.2]
M <- M[, missing_percentage(M) <= 0.2]
DM <- DM[, missing_percentage(DM) <= 0.2]

# 使用KNN填补缺失值小于20%的列
D <- kNN(D, k = 5, imp_var = FALSE)
M <- kNN(M, k = 5, imp_var = FALSE)
DM <- kNN(DM, k = 5, imp_var = FALSE)

# 清理列名函数（替换为.）
clean_column_names_dot <- function(data) {
  colnames(data) <- gsub("[ /']", ".", colnames(data))
  return(data)
}

# 对D, M, DM应用清理列名
D <- clean_column_names_dot(D)
M <- clean_column_names_dot(M)
DM <- clean_column_names_dot(DM)

# 获取变量的Z值
get_z_values <- function(data) {
  z_values <- sapply(names(data)[-c(ncol(data)-1, ncol(data))], function(var) {
    tryCatch({
      # 判断是否为二分类变量（0, 1）
      if (all(data[[var]] %in% c(0, 1))) {
        # 二分类变量不进行标准化
        model <- coxph(Surv(time, event) ~ data[[var]], data = data)
      } else {
        # 连续变量进行标准化
        standardized_var <- scale(data[[var]])
        model <- coxph(Surv(time, event) ~ standardized_var, data = data)
      }
      z_value <- summary(model)$coefficients[,"z"]
      return(z_value)
    }, error = function(e) {
      # 如果模型不收敛或有其他错误，返回NA
      return(NA)
    }, warning = function(w) {
      # 捕获警告信息并返回NA
      return(NA)
    })
  })
  
  # 删除返回NA的变量
  z_values <- z_values[!is.na(z_values)]
  
  # 创建一个数据框并按Z值绝对值排序
  z_df <- data.frame(Variable = names(z_values), Z = z_values)
  z_df <- z_df[order(-abs(z_df$Z)), ]
  return(z_df)
}

# 计算D, M, DM数据集的Z值并排序
D_z_sorted <- get_z_values(D)
M_z_sorted <- get_z_values(M)
DM_z_sorted <- get_z_values(DM)



# Stepwise Selection Function
stepwise_selection_multiple_cv <- function(data, z_sorted_vars, n_repeats = 100, k = 5) {
  selected_vars_frequency <- setNames(rep(0, length(z_sorted_vars$Variable)), z_sorted_vars$Variable)
  all_results <- list()
  
  for (rep in 1:n_repeats) {
    set.seed(NULL)  # 保证每次生成不同的随机数序列
    cat("重复次数:", rep, "\n")  # 显示目前是第几次重复
    
    folds <- createFolds(data$event, k = k, list = TRUE)
    best_cindex_test <- 0
    selected_vars <- character(0)
    final_train_cindex <- 0
    final_test_cindex <- 0
    final_brier_all_train <- NULL
    final_brier_all_test <- NULL
    final_auc_all_train <- NULL
    final_auc_all_test <- NULL
    
    for (i in 1:length(z_sorted_vars$Variable)) {
      candidate_vars <- setdiff(z_sorted_vars$Variable, selected_vars)
      best_var <- NULL
      best_var_cindex_test <- 0
      
      best_train_cindex_per_var <- c()
      best_test_cindex_per_var <- c()
      brier_train_per_var <- NULL
      brier_test_per_var <- NULL
      auc_train_per_var <- NULL
      auc_test_per_var <- NULL
      
      for (var in candidate_vars) {
        current_vars <- c(selected_vars, var)
        formula <- as.formula(paste("Surv(time, event) ~", paste(current_vars, collapse = " + ")))
        
        train_cindex_folds <- c()
        test_cindex_folds <- c()
        brier_train_folds <- matrix(0, ncol = length(seq(0.25, 4.75, by = 0.25)), nrow = k)
        brier_test_folds <- matrix(0, ncol = length(seq(0.25, 4.75, by = 0.25)), nrow = k)
        auc_train_folds <- matrix(0, ncol = length(seq(0.25, 4.75, by = 0.25)), nrow = k)
        auc_test_folds <- matrix(0, ncol = length(seq(0.25, 4.75, by = 0.25)), nrow = k)
        converged <- TRUE
        
        for (j in 1:k) {
          fold <- folds[[j]]
          train_data <- data[-fold, ]
          test_data <- data[fold, ]
          model <- tryCatch(
            coxph(formula, data = train_data),
            error = function(e) {
              converged = FALSE
              return(NULL)
            }
          )
          
          if (is.null(model) || !converged) {
            train_cindex = 0
            test_cindex = 0
          } else {
            surv_obj_train <- Surv(train_data$time, train_data$event)
            surv_obj_test <- Surv(test_data$time, test_data$event)
            predictors_train <- predict(model, newdata = train_data, type = "risk")
            predictors_test <- predict(model, newdata = test_data, type = "risk")
            train_cindex <- 1 - concordance(surv_obj_train ~ predictors_train)$concordance
            test_cindex <- 1 - concordance(surv_obj_test ~ predictors_test)$concordance
            
            brier_train <- predErr(surv_obj_train, surv_obj_train, predictors_train, predictors_train, seq(0.25, 4.75, by = 0.25), type = "brier")$error
            brier_test <- predErr(surv_obj_test, surv_obj_test, predictors_test, predictors_test, seq(0.25, 4.75, by = 0.25), type = "brier")$error
            auc_train <- survAUC::AUC.uno(surv_obj_train, surv_obj_train, predictors_train, seq(0.25, 4.75, by = 0.25))$auc
            auc_test <- survAUC::AUC.uno(surv_obj_test, surv_obj_test, predictors_test, seq(0.25, 4.75, by = 0.25))$auc
          }
          
          train_cindex_folds <- c(train_cindex_folds, train_cindex)
          test_cindex_folds <- c(test_cindex_folds, test_cindex)
          brier_train_folds[j, ] <- brier_train
          brier_test_folds[j, ] <- brier_test
          auc_train_folds[j, ] <- auc_train
          auc_test_folds[j, ] <- auc_test
        }
        
        mean_cindex_test <- mean(test_cindex_folds)
        if (mean_cindex_test > best_var_cindex_test) {
          best_var_cindex_test <- mean_cindex_test
          best_train_cindex_per_var <- mean(train_cindex_folds)
          best_test_cindex_per_var <- mean(test_cindex_folds)
          brier_train_per_var <- colMeans(brier_train_folds)
          brier_test_per_var <- colMeans(brier_test_folds)
          auc_train_per_var <- colMeans(auc_train_folds)
          auc_test_per_var <- colMeans(auc_test_folds)
          best_var <- var
        }
      }
      
      if (best_var_cindex_test > best_cindex_test) {
        best_cindex_test <- best_var_cindex_test
        final_train_cindex <- best_train_cindex_per_var
        final_test_cindex <- best_test_cindex_per_var
        final_brier_all_train <- brier_train_per_var
        final_brier_all_test <- brier_test_per_var
        final_auc_all_train <- auc_train_per_var
        final_auc_all_test <- auc_test_per_var
        selected_vars <- c(selected_vars, best_var)
        selected_vars_frequency[best_var] <- selected_vars_frequency[best_var] + 1
        cat("成功加入变量:", best_var, "当前最佳C-index (测试集):", best_cindex_test, "\n")  # 已加入变量
      } else {
        cat("C-index未增加，停止选择。\n")
        break
      }
    }
    
    all_results[[rep]] <- list(selected_vars = selected_vars, metrics = list(
      train_cindex = final_train_cindex,
      test_cindex = final_test_cindex,
      brier_values_train = final_brier_all_train,
      brier_values_test = final_brier_all_test,
      auc_values_train = final_auc_all_train,
      auc_values_test = final_auc_all_test
    ))
    cat(rep,"次运行结果: Train_Cindex=", final_train_cindex, "Test_Cindex=", final_test_cindex, "\n")  # 显示当前训练集和验证集的C-index
  }
  
  # Calculate frequency of selected variables
  selected_vars_frequency_sorted <- sort(selected_vars_frequency, decreasing = TRUE)
  
  # Aggregate结果
  metrics_summary <- data.frame(
    Train_Cindex_Mean = mean(sapply(all_results, function(x) x$metrics$train_cindex)),
    Train_Cindex_SE = sd(sapply(all_results, function(x) x$metrics$train_cindex)) / sqrt(n_repeats),
    Test_Cindex_Mean = mean(sapply(all_results, function(x) x$metrics$test_cindex)),
    Test_Cindex_SE = sd(sapply(all_results, function(x) x$metrics$test_cindex)) / sqrt(n_repeats),
    Brier_Mean_Train = rowMeans(do.call(cbind, lapply(all_results, function(x) x$metrics$brier_values_train))),
    Brier_SE_Train = apply(do.call(cbind, lapply(all_results, function(x) x$metrics$brier_values_train)), 1, sd) / sqrt(n_repeats),
    Brier_Mean_Test = rowMeans(do.call(cbind, lapply(all_results, function(x) x$metrics$brier_values_test))),
    Brier_SE_Test = apply(do.call(cbind, lapply(all_results, function(x) x$metrics$brier_values_test)), 1, sd) / sqrt(n_repeats),
    AUC_Mean_Train = rowMeans(do.call(cbind, lapply(all_results, function(x) x$metrics$auc_values_train))),
    AUC_SE_Train = apply(do.call(cbind, lapply(all_results, function(x) x$metrics$auc_values_train)), 1, sd) / sqrt(n_repeats),
    AUC_Mean_Test = rowMeans(do.call(cbind, lapply(all_results, function(x) x$metrics$auc_values_test))),
    AUC_SE_Test = apply(do.call(cbind, lapply(all_results, function(x) x$metrics$auc_values_test)), 1, sd) / sqrt(n_repeats)
  )
  
  times <- seq(0.25, 4.75, by = 0.25)
  all_brier_values_train <- lapply(all_results, function(x) x$metrics$brier_values_train)
  all_brier_values_test <- lapply(all_results, function(x) x$metrics$brier_values_test)
  all_auc_values_train <- lapply(all_results, function(x) x$metrics$auc_values_train)
  all_auc_values_test <- lapply(all_results, function(x) x$metrics$auc_values_test)
  
  names(all_brier_values_train) <- paste0("Iteration_", 1:n_repeats)
  names(all_brier_values_test) <- paste0("Iteration_", 1:n_repeats)
  names(all_auc_values_train) <- paste0("Iteration_", 1:n_repeats)
  names(all_auc_values_test) <- paste0("Iteration_", 1:n_repeats)
  
  all_train_cindex <- sapply(all_results, function(x) x$metrics$train_cindex)
  all_test_cindex <- sapply(all_results, function(x) x$metrics$test_cindex)
  
  list(
    metrics = metrics_summary,
    feature_frequency = selected_vars_frequency_sorted,
    train_cindex_values = all_train_cindex,
    test_cindex_values = all_test_cindex,
    brier_values_train = all_brier_values_train,
    brier_values_test = all_brier_values_test,
    auc_values_train = all_auc_values_train,
    auc_values_test = all_auc_values_test,
    times = times
  )
}

# 对D, M, DM数据集执行逐步变量选择，指定初始变量
result_D <- stepwise_selection_multiple_cv(D, D_z_sorted, n_repeats = 100, k = 5)
result_M <- stepwise_selection_multiple_cv(M, M_z_sorted, n_repeats = 100, k = 5)
result_DM <- stepwise_selection_multiple_cv(DM, DM_z_sorted, n_repeats = 100, k = 5)


# 计算 MAGGIC score 的 C-index 和 iAUC
calculate_MAGGIC_metrics <- function(data) {
  formula <- as.formula("Surv(time, event) ~ MAGGIC")
  model <- coxph(formula, data = data)
  cindex <- 1 - concordance(Surv(data$time, data$event) ~ predict(model, newdata = data, type = "risk"))$concordance
  
  times <- seq(0.25, 4.75, by = 0.25)
  auc <- survAUC::AUC.uno(Surv(data$time, data$event), Surv(data$time, data$event), predict(model, newdata = data, type = "risk"), times)$auc
  
  list(cindex = cindex, auc = auc, times = times)
}

MAGGIC_metrics <- calculate_MAGGIC_metrics(user_data)

MAGGIC_cindex <- MAGGIC_metrics$cindex
MAGGIC_auc <- MAGGIC_metrics$auc
times <- MAGGIC_metrics$times

# Violin Plot

# 创建 Violin Plot 所需的数据
violin_data <- data.frame(
  Cindex = c(result_D$test_cindex_values, result_M$test_cindex_values, result_DM$test_cindex_values),
  Model = rep(c("D", "M", "DM"), each = 100)
)

# 添加 MAGGIC score 的辅助线
ggplot(violin_data, aes(x = Model, y = Cindex, fill = Model)) + 
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white") +
  geom_hline(yintercept = MAGGIC_cindex, linetype = "dashed", color = "red") +  # 添加 MAGGIC score 的辅助线
  labs(title = "Test C-index Distribution for Each Model", y = "Test C-index", x = "Model") +
  theme_minimal()

# iAUC 曲线图

# 合并 iAUC 数据并创建数据框
merge_cv_metrics_iAUC <- function(times, auc_mean, auc_se, model_name) {
  data.frame(
    Time = times,
    Mean = auc_mean,
    SE = auc_se,
    Model = model_name
  )
}

df_D <- merge_cv_metrics_iAUC(times, result_D$metrics$AUC_Mean_Test, result_D$metrics$AUC_SE_Test, "Model D")
df_M <- merge_cv_metrics_iAUC(times, result_M$metrics$AUC_Mean_Test, result_M$metrics$AUC_SE_Test, "Model M")
df_DM <- merge_cv_metrics_iAUC(times, result_DM$metrics$AUC_Mean_Test, result_DM$metrics$AUC_SE_Test, "Model DM")
df_MAGGIC <- data.frame(Time = times, Mean = MAGGIC_auc, SE = rep(0, length(times)), Model = "MAGGIC")

# 合并所有数据
df_all <- rbind(df_D, df_M, df_DM, df_MAGGIC)

# 可视化 iAUC 随时间变化的情况
ggplot(df_all, aes(x = Time, y = Mean, color = Model, fill = Model)) +
  geom_line(size = 1) +
  geom_ribbon(data = subset(df_all, Model != "MAGGIC"), aes(ymin = Mean - SE, ymax = Mean + SE, fill = Model), alpha = 0.2) +  # 排除 MAGGIC 的误差带
  labs(title = "iAUC Over Time", y = "iAUC", x = "Time") +
  theme_minimal() +
  theme(legend.position = "bottom")

# 从 feature_frequency 中挑选前 N 个特征
N <- 24  # 你可以根据需要调整前 N 个特征的数量

get_top_n_features <- function(feature_frequency, n) {
  names(feature_frequency)[1:n]
}

# 挑选 D, M, DM 数据集的前 N 个特征
top_n_features_D <- get_top_n_features(result_D$feature_frequency, N)
top_n_features_M <- get_top_n_features(result_M$feature_frequency, N)
top_n_features_DM <- get_top_n_features(result_DM$feature_frequency, N)




# 函数用于在全部的数据集上做 Cox 回归并计算 C-index 和 iAUC
evaluate_cox_model <- function(data, selected_vars) {
  formula <- as.formula(paste("Surv(time, event) ~", paste(selected_vars, collapse = " + ")))
  model <- coxph(formula, data = data)
  
  cindex <- 1 - concordance(Surv(data$time, data$event) ~ predict(model, newdata = data, type = "risk"))$concordance
  
  times <- seq(0.25, 4.75, by = 0.25)
  auc <- survAUC::AUC.uno(Surv(data$time, data$event), Surv(data$time, data$event), predict(model, newdata = data, type = "risk"), times)$auc
  
  list(cindex = cindex, auc = auc, times = times)
}

# 计算 D, M, DM 数据集前 N 个特征的 C-index 和 iAUC
results_top_n_D <- evaluate_cox_model(D, top_n_features_D)
results_top_n_M <- evaluate_cox_model(M, top_n_features_M)
results_top_n_DM <- evaluate_cox_model(DM, top_n_features_DM)

# 打印结果
cat("D dataset - C-index:", results_top_n_D$cindex, "\n")
cat("M dataset - C-index:", results_top_n_M$cindex, "\n")
cat("DM dataset - C-index:", results_top_n_DM$cindex, "\n")

# iAUC 曲线图的合并和绘制
merge_cv_metrics_iAUC <- function(times, auc_mean, model_name) {
  data.frame(
    Time = times,
    Mean = auc_mean,
    Model = model_name
  )
}

df_top_n_D <- merge_cv_metrics_iAUC(results_top_n_D$times, results_top_n_D$auc, "Top N Model D")
df_top_n_M <- merge_cv_metrics_iAUC(results_top_n_M$times, results_top_n_M$auc, "Top N Model M")
df_top_n_DM <- merge_cv_metrics_iAUC(results_top_n_DM$times, results_top_n_DM$auc, "Top N Model DM")
df_MAGGIC <- data.frame(Time = times, Mean = MAGGIC_auc, Model = "MAGGIC")

# 合并所有数据
df_all_iAUC <- rbind(df_top_n_D, df_top_n_M, df_top_n_DM, df_MAGGIC)

# 可视化 iAUC 随时间变化的情况
ggplot(df_all_iAUC, aes(x = Time, y = Mean, color = Model, fill = Model)) +
  geom_line(size = 1) +
  labs(title = "iAUC Over Time for Top N Features", y = "iAUC", x = "Time") +
  theme_minimal() +
  theme(legend.position = "bottom")

# 函数用于在全部的数据集上做 Cox 回归
train_cox_model <- function(data, selected_vars) {
  formula <- as.formula(paste("Surv(time, event) ~", paste(selected_vars, collapse = " + ")))
  coxph(formula, data = data)
}

# 训练 Cox 回归模型
cox_model_D <- train_cox_model(D, top_n_features_D)
cox_model_M <- train_cox_model(M, top_n_features_M)
cox_model_DM <- train_cox_model(DM, top_n_features_DM)


saveRDS(cox_model_D, "cox_model_D.rds")
saveRDS(cox_model_M, "cox_model_M.rds")
saveRDS(cox_model_DM, "cox_model_DM.rds")
# 创建 Excel 文件
wb <- createWorkbook()

# 保存 C-index 数据
cindex_data <- data.frame(
  D = result_D$test_cindex_values,
  M = result_M$test_cindex_values,
  DM = result_DM$test_cindex_values
)

addWorksheet(wb, "C-index")
writeData(wb, "C-index", cindex_data, rowNames = TRUE)

# 保存 iAUC 数据
save_iAUC_data <- function(result, sheet_name) {
  iAUC_data <- do.call(rbind, result$auc_values_test)
  rownames(iAUC_data) <- paste0("Iteration_", 1:100)
  colnames(iAUC_data) <- paste0("Time_", result$times)
  
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, iAUC_data, rowNames = TRUE)
}

save_iAUC_data(result_D, "iAUC_D")
save_iAUC_data(result_M, "iAUC_M")
save_iAUC_data(result_DM, "iAUC_DM")

# 保存 Excel 文件
saveWorkbook(wb, "results.xlsx", overwrite = TRUE)

