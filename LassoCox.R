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

# 单因素Cox回归筛选显著变量
univariate_cox <- function(data) {
  significant_vars <- c()
  
  for (var in colnames(data)[-c(ncol(data)-1, ncol(data))]) {
    formula <- as.formula(paste("Surv(time, event) ~", var))
    
    # 捕获不收敛的变量
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
      
      # 确保p_value是单个值
      if (length(p_value) == 1 && p_value < 0.05) {
        significant_vars <- c(significant_vars, var)
      }
    }
  }
  
  return(significant_vars)
}


# 对D, M, DM进行单因素Cox回归筛选
sig_vars_D <- univariate_cox(D)
sig_vars_M <- univariate_cox(M)
sig_vars_DM <- univariate_cox(DM)

# Lasso-Cox回归和10折交叉验证
lasso_cox_cv <- function(data, significant_vars) {
  # 提取并检查数据矩阵
  x <- data[, significant_vars]
  
  # 确保所有变量是数值型，并处理NA/NaN/Inf
  x <- lapply(x, function(col) as.numeric(as.character(col)))
  x <- as.matrix(do.call(cbind, x))
  
  # 检查NA/NaN/Inf并删除包含这些值的行
  if (any(is.na(x) | is.infinite(x))) {
    warning("Removing rows with NA/NaN/Inf values.")
    x <- x[complete.cases(x), ]
  }
  
  # 创建Surv对象
  y <- Surv(data$time, data$event)
  
  # 执行Lasso-Cox回归和10折交叉验证
  cvfit <- cv.glmnet(x, y, family = "cox", alpha = 1, nfolds = 10)
  best_lambda <- cvfit$lambda.min
  final_model <- glmnet(x, y, family = "cox", alpha = 1, lambda = best_lambda)
  
  # 提取选择的变量
  selected_vars <- rownames(coef(final_model))[which(coef(final_model) != 0)]
  
  return(list(selected_vars = selected_vars, cvfit = cvfit))
}

# 对D, M, DM进行Lasso-Cox回归和10折交叉验证
lasso_D <- lasso_cox_cv(D, sig_vars_D)
lasso_M <- lasso_cox_cv(M, sig_vars_M)
lasso_DM <- lasso_cox_cv(DM, sig_vars_DM)


# 计算C-index, Brier Score和AUC
calculate_metrics <- function(data, selected_vars) {
  # 提取并构建数据矩阵
  x <- data[, selected_vars]
  x <- as.matrix(do.call(cbind, lapply(x, function(col) as.numeric(as.character(col)))))
  
  # 创建Surv对象
  y <- Surv(data$time, data$event)
  
  # 构建Cox回归模型，设置x=TRUE确保设计矩阵保存
  cox_model <- coxph(Surv(time, event) ~ ., data = data[, c(selected_vars, "time", "event")], x = TRUE)
  
  # 计算C-index
  cindex <- summary(cox_model)$concordance[1]
  
  # 计算Brier Score
  lp <- predict(cox_model)  # 线性预测
  times <- seq(0.25, 4.75, by = 0.25)  # 指定时间点
  brier_scores <- predErr(Surv(data$time, data$event), Surv(data$time, data$event), lp, lp, times = times, type = "brier")
  
  # 计算AUC
  auc <- AUC.uno(Surv.rsp = y, Surv.rsp.new = y, lpnew = lp, times = times)
  
  # 返回C-index, Brier Score, AUC
  return(list(model = cox_model,cindex = cindex, brier = brier_scores$error, auc = auc$auc, iauc = auc$iauc))
}

# 对D, M, DM数据集应用所选变量，构建Cox回归模型并计算指标
metrics_D <- calculate_metrics(D, lasso_D$selected_vars)
metrics_M <- calculate_metrics(M, lasso_M$selected_vars)
metrics_DM <- calculate_metrics(DM, lasso_DM$selected_vars)


# 加载外源验证集数据
file_path_validation <- "C:/Users/fenggu/University of Michigan Dropbox/Feng Gu/TriSeg/AM_data_Fengedited.xlsx" 
validation_data <- read_excel(file_path_validation)

# 预处理验证集中的分类变量和数值变量，保持与训练集一致
validation_data$Race <- ifelse(validation_data$Race == "Caucasian", 1, 0)
validation_data$Smoking <- ifelse(validation_data$Smoking == "Current", 1, 0)
validation_data$Alcohol <- ifelse(validation_data$Alcohol == "Yes", 1, 0)
validation_data$Drug <- ifelse(validation_data$Drug == "Yes", 1, 0)
validation_data$HGB <- as.numeric(as.character(validation_data$HGB))

# 重命名时间和事件列，保持一致
colnames(validation_data)[(ncol(validation_data)-5):(ncol(validation_data)-3)] <- c("time", "event", "MAGGIC")

# 创建D, M, DM子集
D_val <- cbind(validation_data[, 4:70], validation_data[, c("time", "event")])
M_val <- cbind(validation_data[, c(71:150)], validation_data[, c("time", "event")])
DM_val <- cbind(validation_data[, c(4:150)], validation_data[, c("time", "event")])

# 去除全空的列
D_val <- D_val[, colSums(is.na(D_val)) < nrow(D_val)]
M_val <- M_val[, colSums(is.na(M_val)) < nrow(M_val)]
DM_val <- DM_val[, colSums(is.na(DM_val)) < nrow(DM_val)]

# 使用KNN填补缺失值
D_val <- kNN(D_val, k = 5, imp_var = FALSE)
M_val <- kNN(M_val, k = 5, imp_var = FALSE)
DM_val <- kNN(DM_val, k = 5, imp_var = FALSE)


# 清理列名函数（替换为.）
clean_column_names_dot <- function(data) {
  colnames(data) <- gsub("[ /']", ".", colnames(data))
  return(data)
}

# 对D, M, DM应用清理列名
D_val <- clean_column_names_dot(D_val)
M_val <- clean_column_names_dot(M_val)
DM_val <- clean_column_names_dot(DM_val)

# 使用训练好的Cox模型对新数据进行预测并计算C-index
predict_and_evaluate <- function(cox_model, data, selected_vars) {
  # 用已经训练好的模型对新数据集进行预测
  predictors <- predict(cox_model, newdata = data[, selected_vars], type = "risk")
  
  # 构建生存对象
  surv_obj <- Surv(data$time, data$event)
  
  # 计算C-index
  cindex <- 1 - concordance(surv_obj ~ predictors)$concordance
  
  # 计算Brier Score
  times <- seq(0.25, 4.75, by = 0.25)  # 指定时间点
  brier_scores <- predErr(Surv(data$time, data$event), Surv(data$time, data$event), predictors, predictors, times = times, type = "brier")
  
  # 计算AUC
  auc <- AUC.uno(Surv.rsp = surv_obj, Surv.rsp.new = surv_obj, lpnew = predictors, times = times)
  
  # 返回C-index, Brier Score, AUC 和 iAUC
  return(list(cindex = cindex, brier = brier_scores$error, auc = auc$auc, iauc = auc$iauc))
}




# 在外源验证集上评估
cindex_D_val <- predict_and_evaluate(metrics_D$model, D_val, lasso_D$selected_vars)
cindex_M_val <- predict_and_evaluate(metrics_M$model, M_val, lasso_M$selected_vars)
cindex_DM_val <- predict_and_evaluate(metrics_DM$model, DM_val, lasso_DM$selected_vars)

# 打印C-index结果
cat("Validation D model C-index:", cindex_D_val$cindex, "\n")
cat("Validation M model C-index:", cindex_M_val$cindex, "\n")
cat("Validation DM model C-index:", cindex_DM_val$cindex, "\n")


library(openxlsx)

# 将结果写入Excel的函数
write_results_to_excel <- function(file_name, metrics_D, metrics_M, metrics_DM, cindex_D_val, cindex_M_val, cindex_DM_val) {
  
  # 创建一个新的工作簿
  wb <- createWorkbook()
  
  # 创建工作表
  addWorksheet(wb, "Metrics_D")
  addWorksheet(wb, "Metrics_M")
  addWorksheet(wb, "Metrics_DM")
  addWorksheet(wb, "Cindex_Validation")
  
  # 写入 Metrics_D 的结果
  writeData(wb, "Metrics_D", data.frame(Cindex = metrics_D$cindex, 
                                        Brier_Score = metrics_D$brier, 
                                        AUC = metrics_D$auc, 
                                        iAUC = metrics_D$iauc))
  
  # 写入 Metrics_M 的结果
  writeData(wb, "Metrics_M", data.frame(Cindex = metrics_M$cindex, 
                                        Brier_Score = metrics_M$brier, 
                                        AUC = metrics_M$auc, 
                                        iAUC = metrics_M$iauc))
  
  # 写入 Metrics_DM 的结果
  writeData(wb, "Metrics_DM", data.frame(Cindex = metrics_DM$cindex, 
                                         Brier_Score = metrics_DM$brier, 
                                         AUC = metrics_DM$auc, 
                                         iAUC = metrics_DM$iauc))
  

  
  # 保存工作簿
  saveWorkbook(wb, file_name, overwrite = TRUE)
}

# 计算并获取metrics
metrics_D <- calculate_metrics(D, lasso_D$selected_vars)
metrics_M <- calculate_metrics(M, lasso_M$selected_vars)
metrics_DM <- calculate_metrics(DM, lasso_DM$selected_vars)

# 验证集上的C-index
metrics_D  <- predict_and_evaluate(metrics_D$model, D_val, lasso_D$selected_vars)
metrics_M <- predict_and_evaluate(metrics_M$model, M_val, lasso_M$selected_vars)
metrics_DM <- predict_and_evaluate(metrics_DM$model, DM_val, lasso_DM$selected_vars)

# 将结果写入Excel文件
write_results_to_excel("model_evaluation_results.xlsx", metrics_D, metrics_M, metrics_DM)


# 提取Cox模型的总结
summary_cox <- summary(metrics_DM$model)

# 提取Z值
z_values <- summary_cox$coefficients[, "z"]

# 创建包含变量名和Z值的表
z_table <- data.frame(
  variable = rownames(summary_cox$coefficients),
  z_value = z_values
)

# 按Z值的绝对值排序
z_table_sorted <- z_table[order(abs(z_table$z_value), decreasing = TRUE), ]

# 查看排序后的表
print(z_table_sorted)
