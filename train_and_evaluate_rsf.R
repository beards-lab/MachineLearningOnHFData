# 定义函数以提取变量重要性
get_feature_importance <- function(model) {
  importance <- vimp(model)
  importance_data <- data.frame(Feature = names(importance$importance), Importance = importance$importance)
  importance_data
}
train_and_evaluate_rsf <- function(data, model_name, params, selected_features = NULL, predict_times = c(1), iteration) {
  # 移除包含缺失值的行
  data <- data[complete.cases(data), ]
  
  # 如果提供了selected_features，则选择相应的列
  if (!is.null(selected_features)) {
    data <- data[, c("time", "event", selected_features)]
  }
  
  # 设置rfsrc模型参数的列表
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
  
  # 处理nodedepth参数
  if (!is.null(params$nodedepth)) {
    if (params$nodedepth == "NoLimit") {
      rfsrc_params$nodedepth <- NULL
    } else {
      # 尝试将nodedepth转换为数值型
      rfsrc_params$nodedepth <- as.numeric(params$nodedepth)
      
      # 检查是否转换成功
      if (is.na(rfsrc_params$nodedepth)) {
        stop("nodedepth must be numeric or 'NoLimit'")
      }
    }
  }
  
  # 训练模型
  best_model <- do.call(rfsrc, rfsrc_params)
  
  # 获取特征重要性
  feature_importance <- get_feature_importance(best_model)
  
  # 计算袋外C-index
  oob_cindex <- 1 - best_model$err.rate[length(best_model$err.rate)]
  
  # 使用survAUC包计算袋外Brier Score和AUC
  times <- seq(0.25, 4.75, by = 0.25)
  Surv.rsf <- Surv(data$time, data$event)
  RSF.pred <- best_model$predicted.oob
  
  brier_scores <- predErr(Surv.rsf, Surv.rsf, RSF.pred, RSF.pred, times, type = "brier")
  auc_scores <- AUC.uno(Surv.rsf, Surv.rsf, RSF.pred, times)
  
  brier_data <- data.frame(Time = brier_scores$times, BrierScore = brier_scores$error, Model = model_name)
  auc_data <- data.frame(Time = auc_scores$times, AUC = auc_scores$auc, Model = model_name)
  iAUC_value <- auc_scores$iauc  # Integrative AUC 值
  
  # 计算指定时间点的ROC曲线数据
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
    ROCData = roc_data
  ))
}


# 定义函数以选择重要特征
select_top_features <- function(data, params, reps = 100, threshold = 0.001, percentage = 0.7, top_n = NULL) {
  feature_importance_list <- list()
  
  for (i in 1:reps) {
    cat("Iteration:", i, "\n")
    results <- train_and_evaluate_rsf(data, "temp", params)
    feature_importance_list <- append(feature_importance_list, list(results$FeatureImportance))
  }
  
  all_importance_data <- do.call(rbind, feature_importance_list)
  
  feature_counts <- all_importance_data %>%
    group_by(Feature) %>%
    summarise(Count = sum(Importance > threshold), Total = n(), Percentage = Count / Total) %>%
    filter(Percentage >= percentage) %>%
    arrange(desc(Percentage)) # 按照Percentage降序排列
  
  if (!is.null(top_n)) {
    feature_counts <- feature_counts %>% head(top_n)
  }
  
  selected_features <- feature_counts$Feature
  
  return(selected_features)
}
