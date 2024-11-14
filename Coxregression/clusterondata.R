# 清除环境中的所有变量
rm(list = ls())

# 安装并加载必要的软件包
necessary_packages <- c("clustMixType", "FactoMineR", "factoextra", "cluster", "openxlsx", "fpc","dplyr","clusterSim")

installed_packages <- rownames(installed.packages())

for (pkg in necessary_packages) {
  if (!(pkg %in% installed_packages)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# 读取数据
file_path <- "C:/Users/fenggu/University of Michigan Dropbox/Feng Gu/GitHub/MachineLearningOnHFData/test.xlsx"  # 请替换为实际路径
user_data <- readxl::read_excel(file_path)

# 选择第三列到最后一列
user_data_selected <- user_data[, 3:ncol(user_data)]
user_data_processed <- user_data_selected %>%
  mutate(across(where(is.character), as.factor))
set.seed(42)
# FAMD分析
famd_result <- FAMD(user_data_processed, ncp = 10, graph = FALSE)
# 提取主成分结果
principal_components <- as.data.frame(famd_result$ind$coord)

fviz_screeplot(famd_result)
fviz_famd_ind(famd_result)
fviz_famd_var(famd_result)


# 计算Gower距离
gower_dist <- daisy(user_data_processed, metric = "gower")
gower_dist_matrix <- as.dist(gower_dist)

evaluate_clusters <- function(data, dist_matrix, max_clusters, B = 50) {
  sil_width <- numeric(max_clusters)
  ch_index <- numeric(max_clusters)
  gap_stat <- numeric(max_clusters)
  dbi_index <- numeric(max_clusters)
  
  for (k in 2:max_clusters) {
    # 进行层次聚类
    hc <- hclust(dist_matrix, method = "ward.D2")
    clusters <- cutree(hc, k = k)
    
    # 计算轮廓系数
    sil <- silhouette(clusters, dist_matrix)
    sil_width[k] <- mean(sil[, 3])
    
    # 计算DBI指数
    numeric_data <- data
    for (col in seq_along(data)) {
      if (is.factor(data[[col]])) {
        numeric_data[[col]] <- as.numeric(data[[col]])
      }
    }
    dbi_index[k] <- index.DB(as.matrix(numeric_data), clusters, d = dist_matrix, centrotypes = "centroids")$DB
    
    # 计算Gap统计量
    gap_stat_res <- clusGap(as.matrix(numeric_data), FUN = function(x, k) {
      list(cluster = cutree(hclust(dist_matrix, method = "ward.D2"), k = k))
    }, K.max = max_clusters, B = B)
    gap_stat[k] <- gap_stat_res$Tab[k, "gap"]
  }
  
  return(list(sil = sil_width, dbi = dbi_index, gap = gap_stat))
}

# 评估层次聚类的最佳聚类数目
max_clusters <- 10  # 最大聚类数，可以根据需要调整
hc_results <- evaluate_clusters(user_data_processed, gower_dist_matrix, max_clusters)

# 可视化评估结果
par(mfrow = c(1, 3))

# 轮廓系数
plot(2:max_clusters, hc_results$sil[2:max_clusters], type = "b", xlab = "Number of Clusters", ylab = "Silhouette Width", main = "Silhouette Coefficient")

# DBI指数
plot(2:max_clusters, hc_results$dbi[2:max_clusters], type = "b", xlab = "Number of Clusters", ylab = "DBI Index", main = "DBI Index")

# Gap统计量
plot(2:max_clusters, hc_results$gap[2:max_clusters], type = "b", xlab = "Number of Clusters", ylab = "Gap Statistic", main = "Gap Statistic")

par(mfrow = c(1, 1))

# 根据轮廓系数或Calinski-Harabasz指数确定最佳聚类数
optimal_clusters_hc <- 2

# 进行最终的层次聚类
hc_final <- hclust(gower_dist_matrix, method = "ward.D2")
final_clusters_hc <- cutree(hc_final, k = optimal_clusters_hc)
par(mfrow = c(1, 1))
plot(hc_final, labels = FALSE, main = "Hierarchical Clustering Dendrogram",xlab = "", sub = "", cex = 0.6)
rect.hclust(hc_final, k = optimal_clusters_hc, border = "red")
# K原型聚类 (k-prototypes clustering)
# 定义评估函数
optimal_clusters_kproto <- function(data, max_clusters, B = 50) {
  sil_width <- numeric(max_clusters)
  dbi_index <- numeric(max_clusters)
  gap_stat <- numeric(max_clusters)
  wss <- numeric(max_clusters)
  numeric_data <- data
  for (col in seq_along(data)) {
    if (is.factor(data[[col]])) {
      numeric_data[[col]] <- as.numeric(data[[col]])
    }
  }
  # 计算Gap统计量
  gap_stat_res <- clusGap(as.matrix(numeric_data), FUN = function(x, k) {
    kproto_res <- kproto(data, k)
    list(cluster = kproto_res$cluster)
  }, K.max = max_clusters, B = B)
  for (k in 2:max_clusters) {
    # K原型聚类
    kproto_result <- kproto(data, k)
    clusters <- kproto_result$cluster
    
    # 计算轮廓系数
    # 注意，这里仍然用Gower距离来计算轮廓系数
    gower_dist <- daisy(data, metric = "gower")
    sil <- silhouette(clusters, as.dist(gower_dist))
    sil_width[k] <- mean(sil[, 3])
    
    # 计算DBI指数
  
    dbi_index[k] <- index.DB(as.matrix(numeric_data), clusters, d = as.dist(gower_dist), centrotypes = "centroids")$DB
    

    gap_stat[k] <- gap_stat_res$Tab[k, "gap"]
    
    # 计算WSS（肘部法则）
    wss[k] <- sum(kproto_result$tot.withinss)
  }
  
  return(list(sil = sil_width, dbi = dbi_index, gap = gap_stat, wss = wss))
}

# 评估K原型聚类的最佳聚类数目
kproto_results <- optimal_clusters_kproto(user_data_processed, max_clusters)
# 可视化评估结果
par(mfrow = c(2, 2))

# 轮廓系数
plot(2:max_clusters, kproto_results$sil[2:max_clusters], type = "b", xlab = "Number of Clusters", ylab = "Silhouette Width", main = "Silhouette Coefficient (k-prototypes)")

# DBI指数
plot(2:max_clusters, kproto_results$dbi[2:max_clusters], type = "b", xlab = "Number of Clusters", ylab = "DBI Index", main = "DBI Index (k-prototypes)")

# Gap统计量
plot(2:max_clusters, kproto_results$gap[2:max_clusters], type = "b", xlab = "Number of Clusters", ylab = "Gap Statistic", main = "Gap Statistic (k-prototypes)")

# 肘部法则
plot(2:max_clusters, kproto_results$wss[2:max_clusters], type = "b", xlab = "Number of Clusters", ylab = "Within-Cluster Sum of Squares", main = "Elbow Method (k-prototypes)")

par(mfrow = c(1, 1))

# 根据轮廓系数或Calinski-Harabasz指数确定最佳聚类数
optimal_clusters_kproto <- 2

# 进行最终的K原型聚类
kproto_final <- kproto(user_data_processed, optimal_clusters_kproto)
final_clusters_kproto <- kproto_final$cluster

# 创建最终的数据框只包含主成分结果和聚类结果
final_result <- principal_components
final_result$HC_Cluster <- final_clusters_hc
final_result$KProto_Cluster <- final_clusters_kproto

# 创建一个新的Excel工作簿
wb <- createWorkbook()

# 添加一个工作表
addWorksheet(wb, "Clustered Data")

# 写入数据到工作表
writeData(wb, sheet = "Clustered Data", x = final_result)

# 保存Excel文件
saveWorkbook(wb, file = "C:/Users/fenggu/University of Michigan Dropbox/Feng Gu/GitHub/MachineLearningOnHFData/Final_Clustered_Data.xlsx", overwrite = TRUE)

# 打印完成信息
cat("数据处理及聚类分析已完成，结果已保存至 Final_Clustered_Data.xlsx\n")
