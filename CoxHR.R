# 清除环境中的所有变量
rm(list = ls())

# 确认环境中已经没有变量
print(ls())

# 安装并加载必要的包
required_packages <- c("randomForestSRC", "survival", "Hmisc", "readxl", "ggplot2", "caret", "dplyr", "VIM", "openxlsx","survAUC","survivalROC")
for (package in required_packages) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package)
  }
}

library(randomForestSRC)
library(survival)
library(Hmisc)
library(readxl)
library(ggplot2)
library(caret)
library(dplyr)
library(VIM)
library(openxlsx)
library(survAUC)
library(survivalROC)

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
rsf_vh_selection <- function(data, n_repeats = 25) {
  
  selected_vars_frequency <- setNames(rep(0, (ncol(data) - 2)), colnames(data)[-c(ncol(data)-1, ncol(data))])
  selected_vars_vimp <- setNames(rep(0, (ncol(data) - 2)), colnames(data)[-c(ncol(data)-1, ncol(data))])
  all_results <- list()
  
  for (rep in 1:n_repeats) {
    set.seed(rep)  # 保证每次生成不同的随机数序列
    cat("重复次数:", rep, "\n")
    
    formula <- as.formula("Surv(time, event) ~ .")
    rsf <- rfsrc(formula, data = data, ntree = 100, importance = TRUE)
    
    var_select <- var.select(rsf)
    min_depth <- var_select$varselect
    
    vimp_scores <- min_depth[, "vimp"]  # 获取VIMP得分
    names(vimp_scores) <- rownames(min_depth)  # 为VIMP得分向量命名
    
    # 使用合适的阈值进行选择
    threshold <- var_select$md.obj$threshold
    
    # 使用最小深度阈值选择初始变量
    initial_vars <- rownames(min_depth[min_depth$depth <= threshold, ])
    
    # 逐一添加变量，直至误差不再减少
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
    
    all_results[[rep]] <- all_selected_vars  # 记录每次的选择结果
    
    for (var in all_selected_vars) {
      selected_vars_frequency[var] <- selected_vars_frequency[var] + 1
      selected_vars_vimp[var] <- selected_vars_vimp[var] + vimp_scores[var]
    }
  }
  
  # 计算VIMP的平均值
  selected_vars_vimp <- selected_vars_vimp / n_repeats
  
  # 对选中的变量按VIMP得分排序
  selected_vars_sorted <- sort(selected_vars_vimp, decreasing = TRUE)
  
  return(list(all_results = all_results, selected_vars_sorted = selected_vars_sorted, selected_vars_frequency = selected_vars_frequency))
}



# 计算D, M, DM数据集的变量重要性评分并排序
result_D <- rsf_vh_selection(D)
result_M <- rsf_vh_selection(M)
result_DM <- rsf_vh_selection(DM)



default_params <- list(
  ntree = 1000,         # Number of trees to grow
  nodesize = 15,        # Minimum size of terminal nodes
  nodedepth = NULL,     # Maximum depth to which a tree should be grown (NULL means no limit)
  splitrule = "logrank", # Splitting rule for survival forests
  nsplit = 10,          # Number of random splits to consider for each candidate variable
  samptype = "swor",    # Sampling type: "swor" for without replacement, "swr" for with replacement
  sampsize = function(n) min(1, 0.632 * n) # Sample size for each tree (typically around 63.2% of the original data size)
)

# 加载并调用函数定义文件
source("C:/Users/fenggu/University of Michigan Dropbox/Feng Gu/GitHub/MachineLearningOnHFData/train_and_evaluate_RSF.R")  # 请替换为实际路径
source("C:/Users/fenggu/University of Michigan Dropbox/Feng Gu/GitHub/MachineLearningOnHFData/Coxregression/optimize_rsf_params.R")  # 请替换为实际路径

# 设置重复次数
reps <- 25

# 定义函数来找出每个数据集的最佳特征数量N
find_best_N_for_dataset <- function(reps, default_params, dataset, selected_vars_sorted, dataset_name) {
  cindex_by_N <- data.frame(N = integer(), CIndex = numeric())
  
  # 对每个可能的N进行前向选择
  for (N in 1:50) {
    selected_features <- names(selected_vars_sorted)[1:N]
    
    # 重复多次训练和评估
    for (i in 1:reps) {
      cat("Dataset:", dataset_name, "N =", N, "Iteration:", i, "\n")
      
      # 传递 iteration 参数给 train_and_evaluate_rsf 函数
      results <- train_and_evaluate_rsf(dataset, dataset_name, default_params, selected_features, iteration = i)
      
      cindex_results <- data.frame(N = N, CIndex = results$CIndex)
      cindex_by_N <- rbind(cindex_by_N, cindex_results)
    }
  }
  
  # 计算每个N的平均C-index
  avg_cindex_by_N <- aggregate(CIndex ~ N, data = cindex_by_N, FUN = mean)
  
  # 找到C-index最大的N
  best_N <- avg_cindex_by_N$N[which.max(avg_cindex_by_N$CIndex)]
  return(list(best_N = best_N, avg_cindex_by_N = avg_cindex_by_N))
}


# 运行函数来找出每个数据集的最佳N
Nresult_D <- find_best_N_for_dataset(reps, default_params, D, result_D$selected_vars_sorted, "D")
Nresult_M <- find_best_N_for_dataset(reps, default_params, M, result_M$selected_vars_sorted, "M")
Nresult_DM <- find_best_N_for_dataset(reps, default_params, DM, result_DM$selected_vars_sorted, "DM")

# 最佳的N值
best_N_D <- Nresult_D$best_N
best_N_M <- Nresult_M$best_N
best_N_DM <- Nresult_DM$best_N

cat("Best N for D:", best_N_D, "\n")
cat("Best N for M:", best_N_M, "\n")
cat("Best N for DM:", best_N_DM, "\n")

# 绘制C-index随N变化的图
par(mfrow = c(1, 3)) # 三个子图并排显示

plot(Nresult_D$avg_cindex_by_N$N, Nresult_D$avg_cindex_by_N$CIndex, type = "b", col = "blue", xlab = "Number of Features (N)", ylab = "Average C-Index", main = "C-Index vs N (D)")
plot(Nresult_M$avg_cindex_by_N$N, Nresult_M$avg_cindex_by_N$CIndex, type = "b", col = "red", xlab = "Number of Features (N)", ylab = "Average C-Index", main = "C-Index vs N (M)")
plot(Nresult_DM$avg_cindex_by_N$N, Nresult_DM$avg_cindex_by_N$CIndex, type = "b", col = "green", xlab = "Number of Features (N)", ylab = "Average C-Index", main = "C-Index vs N (DM)")



# 提取前N个最常被选择的特征
get_top_features <- function(result, top_n = NULL) {
  selected_vars_sorted <- result$selected_vars_sorted
  top_selected_vars <- names(selected_vars_sorted)[1:min(top_n, length(selected_vars_sorted))]
  return(top_selected_vars)
}

selected_features_D <- get_top_features(result_D,top_n =best_N_D)
selected_features_M <- get_top_features(result_M,top_n = best_N_M)
selected_features_DM <- get_top_features(result_DM,top_n = best_N_DM)

# 初始化一个列表来存储最终保留的特征
final_features <- selected_features_DM

# 迭代以 _D 结尾的特征
for (feature in selected_features_DM) {
  if (grepl("_D$", feature)) {
    base_name <- gsub("_D$", "", feature)
    corresponding_s <- paste0(base_name, "_S")
    
    # 如果存在对应的 _S 特征，删除 _S 特征
    if (corresponding_s %in% final_features) {
      final_features <- setdiff(final_features, corresponding_s)
    }
  }
  
  # 迭代以 _S 结尾的特征
  if (grepl("_S$", feature)) {
    base_name <- gsub("_S$", "", feature)
    corresponding_d <- paste0(base_name, "_D")
    
    # 如果存在对应的 _D 特征，删除 _S 特征
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


# 输出结果
selected_features_DM
# 设置更细致的超参数网格
param_grid <- list(
  ntree = c(1000,1500,2000),
  mtry = c(4,5,6,7),
  nodesize = c(10,30,50),
  nodedepth = c(1,3,5,10,15,20,"NoLimit"),
  splitrule = c("logrank"),
  samptype = "swor",    # Sampling type: "swor" for without replacement, "swr" for with replacement
  sampsize = function(n) max(1, 0.632 * n) # Sample size for each tree 
)

# 调用优化函数
opt_D <- optimize_rsf_params(D, param_grid,selected_features_D)
opt_M <- optimize_rsf_params(M, param_grid,selected_features_M)
opt_DM <- optimize_rsf_params(DM, param_grid,selected_features_DM)




best_params_D <- opt_D$best_params
best_params_M <- opt_M$best_params
best_params_DM <- opt_DM$best_params
save(best_params_D, best_params_M, best_params_DM,
     selected_features_D, selected_features_M, selected_features_DM,
     file = "saved_variables.RData")
load("saved_variables.RData")
# 运行多次并记录结果
reps <- 25
cindex_results <- data.frame()
brier_results <- data.frame()
auc_results <- data.frame()
iauc_values <- data.frame()
roc_results <- data.frame()

predict_times <- c(0.5,1, 2, 3,5)  # 可以根据需要设置多个预测时间点

for (i in 1:reps) {
  cat("Iteration:", i, "\n")
  
  results_D <- train_and_evaluate_rsf(D, "D", best_params_D, selected_features_D, predict_times, iteration = i)
  cindex_results <- rbind(cindex_results, data.frame(Model = "D", CIndex = results_D$CIndex))
  brier_results <- rbind(brier_results, results_D$BrierScoreData)
  auc_results <- rbind(auc_results, results_D$AUCData)
  iauc_values <- rbind(iauc_values, data.frame(Model = "D", iAUC = results_D$iAUC))
  roc_results <- rbind(roc_results, results_D$ROCData)
  
  results_M <- train_and_evaluate_rsf(M, "M", best_params_M, selected_features_M, predict_times, iteration = i)
  cindex_results <- rbind(cindex_results, data.frame(Model = "M", CIndex = results_M$CIndex))
  brier_results <- rbind(brier_results, results_M$BrierScoreData)
  auc_results <- rbind(auc_results, results_M$AUCData)
  iauc_values <- rbind(iauc_values, data.frame(Model = "M", iAUC = results_M$iAUC))
  roc_results <- rbind(roc_results, results_M$ROCData)
  
  results_DM <- train_and_evaluate_rsf(DM, "DM", best_params_DM, selected_features_DM, predict_times, iteration = i)
  cindex_results <- rbind(cindex_results, data.frame(Model = "DM", CIndex = results_DM$CIndex))
  brier_results <- rbind(brier_results, results_DM$BrierScoreData)
  auc_results <- rbind(auc_results, results_DM$AUCData)
  iauc_values <- rbind(iauc_values, data.frame(Model = "DM", iAUC = results_DM$iAUC))
  roc_results <- rbind(roc_results, results_DM$ROCData)
}

# 函数用于绘制MAGGIC score的ROC曲线并计算AUC
calculate_and_plot_MAGGIC_ROC <- function(data, predict_times) {
  Surv.rsf <- Surv(data$time, data$event)
  risks <- data$MAGGIC
  plot_list <- list()
  
  for (time in predict_times) {
    # 计算指定时间点的ROC曲线和AUC
    roc_curve <- survivalROC(Stime = data$time, status = data$event, marker = risks, predict.time = time, method = "KM")
    auc <- roc_curve$AUC
    
    # 绘制ROC曲线
    p <- ggplot(data.frame(FPR = roc_curve$FP, TPR = roc_curve$TP), aes(x = FPR, y = TPR)) +
      geom_line(size = 1, color = "blue") +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
      labs(title = paste("ROC Curve for MAGGIC Score at", time, "Years"),
           subtitle = paste("AUC =", round(auc, 3)),
           x = "False Positive Rate (1 - Specificity)", 
           y = "True Positive Rate (Sensitivity)") +
      theme_minimal()
    
    plot_list[[as.character(time)]] <- p
    print(p)  # 在R环境中显示图
  }
  
  return(plot_list)
}

# 调用函数
MAGGIC_ROC_plots <- calculate_and_plot_MAGGIC_ROC(user_data, predict_times)

# 分别绘制每个时间点的ROC曲线
for (time in predict_times) {
  roc_time_data <- subset(roc_results, PredictTime == time)
  
  # 计算均值和标准误
  roc_summary <- roc_time_data %>%
    group_by(Model, FPR) %>%
    summarize(
      Mean_TPR = mean(TPR),
      SE_TPR = sd(TPR) / sqrt(n()),
      .groups = 'drop'
    )
  
  # 绘制ROC曲线
  p <- ggplot(roc_summary, aes(x = FPR, y = Mean_TPR, color = Model)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = Mean_TPR - SE_TPR, ymax = Mean_TPR + SE_TPR), alpha = 0.2) +
    labs(title = paste("ROC Curve for", time, "Year Death Prediction"), 
         x = "False Positive Rate (1 - Specificity)", 
         y = "True Positive Rate (Sensitivity)") +
    theme_minimal() +
    scale_color_manual(values = c("D" = "blue", "M" = "green", "DM" = "red"))
  
  print(p)
}
# 计算Cindex
summary_results <- aggregate(CIndex ~ Model, data = cindex_results, function(x) c(mean = mean(x), sd = sd(x)))
print(summary_results)

# 计算 Integrative AUC
integrated_auc_values <- aggregate(iAUC ~ Model, data = iauc_values, function(x) c(mean = mean(x), sd = sd(x)))
print("Integrated AUC values:")
print(integrated_auc_values)

# 汇总 Brier Score 和 AUC 数据
brier_summary <- brier_results %>%
  group_by(Model, Time) %>%
  summarize(Mean_BrierScore = mean(BrierScore), SD_BrierScore = sd(BrierScore), .groups = 'drop')

auc_summary <- auc_results %>%
  group_by(Model, Time) %>%
  summarize(Mean_AUC = mean(AUC), SD_AUC = sd(AUC), .groups = 'drop')

# 绘制 Brier Score 随时间变化的图
ggplot(brier_summary, aes(x = Time, y = Mean_BrierScore, color = Model)) +
  geom_line() +
  geom_ribbon(aes(ymin = Mean_BrierScore - SD_BrierScore, ymax = Mean_BrierScore + SD_BrierScore), 
              alpha = 0.2) +
  labs(title = "Brier Score over Time for Different Models", x = "Time", y = "Brier Score") +
  theme_minimal() +
  scale_color_manual(values = c("D" = "blue", "M" = "green", "DM" = "red"))

# 绘制 AUC 随时间变化图
ggplot(auc_summary, aes(x = Time, y = Mean_AUC, color = Model)) +
  geom_line() +
  geom_ribbon(aes(ymin = Mean_AUC - SD_AUC, ymax = Mean_AUC + SD_AUC), 
              alpha = 0.2) +
  labs(title = "AUC over Time for Different Models (survAUC)", x = "Time", y = "AUC") +
  theme_minimal() +
  scale_color_manual(values = c("D" = "blue", "M" = "green", "DM" = "red"))

# 绘制C-index分布图
ggplot(cindex_results, aes(x = Model, y = CIndex)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(title = "C-index Distribution for Different Models", x = "Model", y = "C-index") +
  theme_minimal()

# 准备数据和特征
datasets <- list(D = D, M = M, DM = DM)
selected_features <- list(D = selected_features_D, M = selected_features_M, DM = selected_features_DM)
best_params <- list(D = best_params_D, M = best_params_M, DM = best_params_DM)

# 用于存储模型结果
rsf_models <- list(D = list(), M = list(), DM = list())

# 训练25个RSF模型
n_repeats <- 25

for (dataset_name in names(datasets)) {
  data <- datasets[[dataset_name]]
  features <- selected_features[[dataset_name]]
  params <- best_params[[dataset_name]]
  
  for (i in 1:n_repeats) {
    cat("训练", dataset_name, "模型", i, "\n")
    
    # 构建公式
    formula <- as.formula(paste("Surv(time, event) ~", paste(features, collapse = " + ")))
    
    # 训练RSF模型
    rsf_models[[dataset_name]][[i]] <- rfsrc(
      formula,
      data = data[, c("time", "event", features)],
      ntree = params$ntree,
      mtry = params$mtry,
      nodesize = params$nodesize,
      importance = "none"
    )
  }
}

# 保存模型到文件
rsf_models_D <- rsf_models$D
rsf_models_M <- rsf_models$M
rsf_models_DM <- rsf_models$DM

save(rsf_models_D, file = "rsf_models_D.RData")
save(rsf_models_M, file = "rsf_models_M.RData")
save(rsf_models_DM, file = "rsf_models_DM.RData")

library(dplyr)
library(tidyr)
library(purrr)
library(openxlsx)

# 假设 auc_results 是你的数据框
sheets <- auc_results %>%
  group_by(Model) %>%
  nest()

# 创建一个Excel文件，并为每个 Model 创建一个子表
wb <- createWorkbook()

# 使用purrr::pwalk来遍历sheets的每一行，并对每个 Model 处理数据并写入Excel
pwalk(list(sheets$Model, sheets$data), function(Model, data) {
  
  # 确保data是一个数据框
  df <- as_tibble(data)
  
  df <- df %>% 
    group_by(Time) %>%  # 先按时间列重排
    mutate(measurement = row_number()) %>%  # 添加测量次数列
    ungroup() %>%
    pivot_wider(names_from = Time, values_from = AUC, values_fill = list(AUC = NA))  # 转换为宽格式
  
  addWorksheet(wb, Model)
  writeData(wb, sheet = Model, df)  # 把数据写入对应表单
})

# 保存Excel文件
saveWorkbook(wb, "auc_results.xlsx", overwrite = TRUE)


# 修改后的函数，直接使用 MAGGIC score 计算 C-index, AUC 和 Brier Score
calculate_MAGGIC_metrics <- function(data) {
  # 提取 MAGGIC score
  risks <- data$MAGGIC
  
  # 计算 C-index
  formula <- as.formula("Surv(time, event) ~ risks")
  cindex <- 1 - concordance(Surv(data$time, data$event) ~ risks)$concordance
  
  # 时间点序列
  times <- seq(0.25, 4.75, by = 0.25)
  
  # 计算 AUC
  auc <- survAUC::AUC.uno(Surv(data$time, data$event), Surv(data$time, data$event), risks, times)$auc
  iauc <- survAUC::AUC.uno(Surv(data$time, data$event), Surv(data$time, data$event), risks, times)$iauc
  
  # 计算 Brier Score
  Surv.rsf <- Surv(data$time, data$event)
  RSF.pred <- risks
  brier_scores <- predErr(Surv.rsf, Surv.rsf, RSF.pred, RSF.pred, times, type = "brier")
  
  list(cindex = cindex, auc = auc,iauc = iauc, brier_scores = brier_scores$error, times = times)
}

# 使用函数计算MAGGIC score的结果
MAGGIC_metrics <- calculate_MAGGIC_metrics(user_data)

# 查看结果
print(MAGGIC_metrics)


MAGGIC_cindex <- MAGGIC_metrics$cindex
MAGGIC_auc <- MAGGIC_metrics$auc
MAGGIC_Brier <- MAGGIC_metrics$brier_scores
times <- MAGGIC_metrics$times
# 假设MAGGIC_metrics[["auc"]] 是你的向量
auc_values <- MAGGIC_metrics[["brier_scores"]]

# 创建一个包含向量的单行数据框
df <- as.data.frame(t(auc_values))
colnames(df) <- paste0("AUC_", seq_along(auc_values))  # 可选：命名列

# 如果您将其添加到Excel文件中，使用现有的Workbook对象
wb <- createWorkbook()

# 添加一个新的工作表
addWorksheet(wb, "MAGGIC_metrics")

# 将数据写入Excel工作表
writeData(wb, sheet = "MAGGIC_metrics", df)

# 保存Excel文件
saveWorkbook(wb, "MAGGIC_metrics.xlsx", overwrite = TRUE)



# 假设 beir_results 是你的数据框
sheets <- brier_results %>%
  group_by(Model) %>%
  nest()

# 创建一个Excel文件，并为每个 Model 创建一个子表
wb <- createWorkbook()

# 使用purrr::pwalk来遍历sheets的每一行，并对每个 Model 处理数据并写入Excel
pwalk(list(sheets$Model, sheets$data), function(Model, data) {
  
  # 确保data是一个数据框
  df <- as_tibble(data)
  
  df <- df %>% 
    group_by(Time) %>%  # 先按时间列重排
    mutate(measurement = row_number()) %>%  # 添加测量次数列
    ungroup() %>%
    pivot_wider(names_from = Time, values_from = BrierScore, values_fill = list(BeirScore = NA))  # 转换为宽格式
  
  addWorksheet(wb, Model)
  writeData(wb, sheet = Model, df)  # 把数据写入对应表单
})

# 保存Excel文件
saveWorkbook(wb, "brier_results.xlsx", overwrite = TRUE)
# 定义基于MAGGIC的特征选择函数
rsf_vh_selection_with_maggic <- function(user_data, M, n_repeats = 25) {
  
  # 确保MAGGIC列存在于user_data中
  if (!"MAGGIC" %in% colnames(user_data)) {
    stop("MAGGIC score is not present in the user_data.")
  }
  
  # 初始特征为MAGGIC
  selected_vars_frequency <- setNames(rep(0, ncol(M) + 1), c("MAGGIC", colnames(M)))
  selected_vars_vimp <- setNames(rep(0, ncol(M) + 1), c("MAGGIC", colnames(M)))
  all_results <- list()
  
  for (rep in 1:n_repeats) {
    set.seed(rep)  # 保证每次生成不同的随机数序列
    cat("重复次数:", rep, "\n")
    
    # 在M集特征基础上手动加入MAGGIC列
    formula_full <- as.formula(paste("Surv(time, event) ~ MAGGIC +", paste(colnames(M), collapse = " + ")))
    
    # 使用包含MAGGIC和M集特征的随机森林模型
    rsf_full <- rfsrc(formula_full, data = user_data, ntree = 100, importance = TRUE)
    
    # 对完整模型进行变量选择，使用var.select
    var_select <- var.select(rsf_full)
    min_depth <- var_select$varselect
    
    vimp_scores <- min_depth[, "vimp"]  # 获取VIMP得分
    names(vimp_scores) <- rownames(min_depth)  # 为VIMP得分向量命名
    
    # 确保MAGGIC的VIMP得分被计算
    if ("MAGGIC" %in% rownames(min_depth)) {
      vimp_scores["MAGGIC"] <- rsf_full$importance["MAGGIC"]
    }
    
    # 逐一添加M集中未被选择的变量，初始从MAGGIC出发
    all_selected_vars <- "MAGGIC"
    oob_cindex <- rsf_full$err.rate[100]  # 基于完整模型的初始误差
    
    for (var in rownames(min_depth)[order(min_depth$depth)]) {
      if (var != "MAGGIC") {
        # 添加新特征到现有特征集
        all_selected_vars <- c(all_selected_vars, var)
        formula <- as.formula(paste("Surv(time, event) ~", paste(all_selected_vars, collapse = " + ")))
        rsf <- rfsrc(formula, data = user_data, ntree = 100, importance = TRUE)
        
        new_oob_cindex <- rsf$err.rate[100]
        
        # 如果新特征导致误差增加，移除该特征
        if (new_oob_cindex >= oob_cindex) {
          all_selected_vars <- setdiff(all_selected_vars, var)
        } else {
          oob_cindex <- new_oob_cindex
        }
      }
    }
    
    all_results[[rep]] <- all_selected_vars  # 记录每次的选择结果
    
    for (var in all_selected_vars) {
      selected_vars_frequency[var] <- selected_vars_frequency[var] + 1
      selected_vars_vimp[var] <- selected_vars_vimp[var] + vimp_scores[var]
    }
  }
  
  # 计算VIMP的平均值
  selected_vars_vimp <- selected_vars_vimp / n_repeats
  
  # 对选中的变量按VIMP得分排序
  selected_vars_sorted <- sort(selected_vars_vimp, decreasing = TRUE)
  
  return(list(all_results = all_results, selected_vars_sorted = selected_vars_sorted, selected_vars_frequency = selected_vars_frequency))
}

# 调用新函数，基于MAGGIC和M子集选择特征
result_M_with_maggic <- rsf_vh_selection_with_maggic(user_data, M)
# 将MAGGIC从user_data中加入到M集
M_with_maggic <- cbind(M, user_data$MAGGIC)
colnames(M_with_maggic)[ncol(M_with_maggic)] <- "MAGGIC"  # 重命名最后一列为MAGGIC

# 检查M_with_maggic是否包含MAGGIC
head(M_with_maggic)
Nresult_M_with_maggic <- find_best_N_for_dataset(reps, default_params, M_with_maggic, result_M_with_maggic$selected_vars_sorted,"MAGGIC_with_M")
# 最佳的N值
best_M_with_maggic <- Nresult_M_with_maggic$best_N
cat("Best N for MAGGIC WITH M:", best_M_with_maggic, "\n")


# 绘制C-index随N变化的图
par(mfrow = c(1, 1)) # 三个子图并排显示

plot(Nresult_M_with_maggic$avg_cindex_by_N$N, Nresult_M_with_maggic$avg_cindex_by_N$CIndex, type = "b", col = "blue", xlab = "Number of Features (N)", ylab = "Average C-Index", main = "C-Index vs N (D)")


# Assuming you have a function to get top features and optimize the model
selected_features_M_with_maggic <- get_top_features(result_M_with_maggic, top_n = 6)
opt_M_with_maggic <- optimize_rsf_params(M_with_maggic, param_grid, selected_features_M_with_maggic)
best_params_M_with_maggic <- opt_M_with_maggic$best_params

# Initialize result dataframes for C-index, Brier score, AUC, etc.
cindex_results <- data.frame()
brier_results <- data.frame()
auc_results <- data.frame()
iauc_values <- data.frame()

# Loop for model evaluation
for (i in 1:reps) {
  cat("Iteration:", i, "\n")
  
  # Train and evaluate the model with MAGGIC
  results_M_with_maggic <- train_and_evaluate_rsf(M_with_maggic, "M_with_maggic", best_params_M_with_maggic, 
                                                  selected_features_M_with_maggic, predict_times, iteration = i)
  
  # Collect results for C-Index, Brier Score, AUC, iAUC
  cindex_results <- rbind(cindex_results, data.frame(Model = "M_with_maggic", CIndex = results_M_with_maggic$CIndex))
  brier_results <- rbind(brier_results, results_M_with_maggic$BrierScoreData)
  auc_results <- rbind(auc_results, results_M_with_maggic$AUCData)
  iauc_values <- rbind(iauc_values, data.frame(Model = "M_with_maggic", iAUC = results_M_with_maggic$iAUC))
}

# Assuming `calculate_MAGGIC_metrics` has already been run for MAGGIC metrics
MAGGIC_metrics <- calculate_MAGGIC_metrics(user_data)
MAGGIC_cindex <- MAGGIC_metrics$cindex
MAGGIC_auc <- MAGGIC_metrics$auc
MAGGIC_Brier <- MAGGIC_metrics$brier_scores
times <- MAGGIC_metrics$times

# Collect Brier Score data
brier_summary <- data.frame(
  Time = rep(times, 2),
  BrierScore = c(MAGGIC_Brier, brier_results$BrierScore),
  Model = factor(rep(c("MAGGIC Only", "M_with_maggic"), each = length(times)))
)

# Plot Brier Score over time
ggplot(brier_summary, aes(x = Time, y = BrierScore, color = Model)) +
  geom_line() +
  geom_ribbon(aes(ymin = BrierScore - sd(BrierScore), ymax = BrierScore + sd(BrierScore)), alpha = 0.2) +
  labs(title = "Brier Score over Time for MAGGIC vs M_with_maggic", x = "Time", y = "Brier Score") +
  theme_minimal() +
  scale_color_manual(values = c("MAGGIC Only" = "blue", "M_with_maggic" = "green"))

# Collect AUC data
auc_summary <- data.frame(
  Time = rep(times, 2),
  AUC = c(MAGGIC_auc, auc_results$AUC),
  Model = factor(rep(c("MAGGIC Only", "M_with_maggic"), each = length(times)))
)

# Plot AUC over time
ggplot(auc_summary, aes(x = Time, y = AUC, color = Model)) +
  geom_line() +
  geom_ribbon(aes(ymin = AUC - sd(AUC), ymax = AUC + sd(AUC)), alpha = 0.2) +
  labs(title = "AUC over Time for MAGGIC vs M_with_maggic", x = "Time", y = "AUC") +
  theme_minimal() +
  scale_color_manual(values = c("MAGGIC Only" = "blue", "M_with_maggic" = "green"))

# Summarize C-index for M_with_maggic
summary_results <- aggregate(CIndex ~ Model, data = subset(cindex_results, Model == "M_with_maggic"), 
                             function(x) c(mean = mean(x), sd = sd(x)))
print(summary_results)

# Summarize integrative AUC for M_with_maggic
integrated_auc_values <- aggregate(iAUC ~ Model, data = subset(iauc_values, Model == "M_with_maggic"), 
                                   function(x) c(mean = mean(x), sd = sd(x)))
print("Integrated AUC values:")
print(integrated_auc_values)
