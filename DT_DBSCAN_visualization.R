library(dbscan)
library(gridExtra)
library(ggplot2)
library(MASS)

# 数据生成参数
n <- 250  # 样本数量
p <- 100  # 总特征数
p1 <- 20  # 前p1个特征有信号
mu_1 <- 0
mu_2 <- 5
sigma_1 <- 1
sigma_2 <- 1
rho <- 0  # 相关性参数

# 生成数据函数
generate_data <- function(n, p, p1, mu_1, mu_2, sigma_1, sigma_2, rho) {
  # 生成协方差矩阵
  Sigma <- rho^abs(matrix(rep(1:p, each=p), p, p, byrow = F) - matrix(rep(1:p, each=p), p, p, byrow = T))
  Sigma <- sqrt(diag(c(rep(sigma_1, p1), rep(sigma_2, p-p1)))) %*% Sigma %*% sqrt(diag(c(rep(sigma_1, p1), rep(sigma_2, p-p1))))
  
  # 生成真实标签
  npor <- 0.6
  n1 <- ceiling(n * npor)  # 第一类样本数
  n0 <- n - n1            # 第二类样本数
  
  # 生成中心点
  loc <- 1:p1
  center0 <- rep(mu_1, p)
  center1 <- rep(mu_1, p)
  center1[loc] <- mu_2
  
  # 生成数据
  X1 <- t(t(matrix(rnorm(n1*p), n1, p) %*% chol(Sigma)) + center1)
  X0 <- t(t(matrix(rnorm(n0*p), n0, p) %*% chol(Sigma)) + center0)
  X <- rbind(X1, X0)
  
  # 真实标签
  true_labels <- c(rep(1, n1), rep(0, n0))
  
  return(list(X = X, true_labels = true_labels))
}

# DT方法的数据分割函数
dt_data_split <- function(X, p1, epi = 0.5, eps = 15, minPts = 3) {
  # 使用DBSCAN进行聚类
  dbscan_result <- dbscan::dbscan(X, eps = eps, minPts = minPts)
  clu <- dbscan_result$cluster
  
  # 检查是否找到了聚类
  unique_clusters <- unique(clu[clu != 0])  # 排除噪声点
  if(length(unique_clusters) < 2) {
    stop("DBSCAN未能找到足够的聚类")
  }
  
  # 计算每个聚类的样本数量
  cluster_sizes <- sapply(unique_clusters, function(c) sum(clu == c))
  
  # 按样本数量排序，选择前两个最大的聚类
  sorted_indices <- order(cluster_sizes, decreasing = TRUE)
  cluster1 <- unique_clusters[sorted_indices[1]]
  cluster2 <- unique_clusters[sorted_indices[2]]
  
  # 获取两个聚类的数据点
  clu1_indices <- which(clu == cluster1)
  clu2_indices <- which(clu == cluster2)
  
  if(length(clu1_indices) == 0 || length(clu2_indices) == 0) {
    stop("聚类结果无效")
  }
  
  # 计算协方差矩阵
  sigma <- (var(X[clu1_indices,]) + var(X[clu2_indices,])) / 2
  
  # 生成测试集和训练集
  X_test <- mvrnorm(n = nrow(X), mu = rep(0, ncol(X)), Sigma = epi * (1 - epi) * sigma) + epi * X
  X_train <- X - X_test
  
  return(list(X_train = X_train, X_test = X_test, 
              original_clusters = clu, 
              cluster1 = cluster1, cluster2 = cluster2))
}

# 生成数据
data_result <- generate_data(n, p, p1, mu_1, mu_2, sigma_1, sigma_2, rho)
X <- data_result$X
true_labels <- data_result$true_labels

# DBSCAN参数设置
# 原始数据集的聚类参数
eps_original <- 13
minPts_original <- 5

# 训练集的聚类参数（可以不同）
eps_train <- 10
minPts_train <- 5

epi <- 0.5  # DT方法的参数

# 执行DT数据分割
dt_result <- dt_data_split(X, p1, epi = epi, eps = eps_original, minPts = minPts_original)
X_train <- dt_result$X_train
X_test <- dt_result$X_test
original_clusters <- dt_result$original_clusters

# 在训练集上进行DBSCAN聚类（使用不同的参数）
dbscan_train <- dbscan::dbscan(X_train, eps = eps_train, minPts = minPts_train)
train_clusters <- dbscan_train$cluster

# 测试集不直接进行聚类，而是直接使用训练集的聚类结果
# 在DT方法中，训练集和测试集是一一对应的，训练集的聚类标签直接对应到测试集
test_clusters <- train_clusters

# 创建数据框用于绘图
plot_data <- data.frame(
  x = X[, 1],  # 第1个特征作为x坐标
  y = X[, 2],  # 第2个特征作为y坐标
  true_label = as.factor(true_labels),
  original_cluster = as.factor(original_clusters),
  train_cluster = as.factor(train_clusters),
  test_cluster = as.factor(test_clusters),
  # 添加数据来源标识
  data_source = "all"  # 默认为全部数据
)

# 创建训练集数据框
train_data <- data.frame(
  x = X_train[, 1],  # 训练集第1个特征作为x坐标
  y = X_train[, 2],  # 训练集第2个特征作为y坐标
  true_label = as.factor(true_labels),
  original_cluster = as.factor(original_clusters),
  train_cluster = as.factor(train_clusters),
  test_cluster = as.factor(test_clusters),
  data_source = "train"
)

# 创建测试集数据框
test_data <- data.frame(
  x = X_test[, 1],  # 测试集第1个特征作为x坐标
  y = X_test[, 2],  # 测试集第2个特征作为y坐标
  true_label = as.factor(true_labels),
  original_cluster = as.factor(original_clusters),
  train_cluster = as.factor(train_clusters),
  test_cluster = as.factor(test_clusters),
  data_source = "test"
)

# 创建真实标签的散点图
plot_true <- ggplot(plot_data, aes(x = x, y = y, color = true_label)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(
    title = "真实标签分布",
    x = "Signal 1",
    y = "Signal 2",
    color = "真实标签"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  scale_color_manual(values = c("0" = "#FF6B6B", "1" = "#4ECDC4"))

# 创建原始数据DBSCAN聚类标签的散点图
plot_original <- ggplot(plot_data, aes(x = x, y = y, color = original_cluster)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(
    title = paste0("原始数据DBSCAN聚类 (eps=", eps_original, ", minPts=", minPts_original, ")"),
    x = "Signal 1",
    y = "Signal 2",
    color = "聚类标签"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  scale_color_manual(values = c("0" = "#999999", "1" = "#FF6B6B", "2" = "#4ECDC4", "3" = "#45B7D1", "4" = "#96CEB4"))

# 创建训练集DBSCAN聚类标签的散点图
plot_train <- ggplot(train_data, aes(x = x, y = y, color = train_cluster)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(
    title = paste0("训练集DBSCAN聚类 (eps=", eps_train, ", minPts=", minPts_train, ")"),
    x = "Signal 1",
    y = "Signal 2",
    color = "聚类标签"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  scale_color_manual(values = c("0" = "#999999", "1" = "#FF6B6B", "2" = "#4ECDC4", "3" = "#45B7D1", "4" = "#96CEB4"))

# 创建测试集DBSCAN聚类标签的散点图
plot_test <- ggplot(test_data, aes(x = x, y = y, color = test_cluster)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(
    title = "测试集聚类结果（与训练集一一对应）",
    x = "Signal 1",
    y = "Signal 2",
    color = "聚类标签"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  scale_color_manual(values = c("0" = "#999999", "1" = "#FF6B6B", "2" = "#4ECDC4", "3" = "#45B7D1", "4" = "#96CEB4"))

# 组合四个图
combined_plot <- grid.arrange(plot_true, plot_original, plot_train, plot_test, ncol = 2, nrow = 2)

# 保存图像
ggsave("DT_DBSCAN_clustering_comparison.png", combined_plot, width = 16, height = 12, dpi = 300)

# 打印聚类结果统计
cat("=== DT方法DBSCAN聚类结果统计 ===\n")
cat("数据生成参数:\n")
cat("样本数量 n =", n, "\n")
cat("总特征数 p =", p, "\n")
cat("信号特征数 p1 =", p1, "\n")
cat("DT参数 epi =", epi, "\n")
cat("原始数据DBSCAN参数: eps =", eps_original, ", minPts =", minPts_original, "\n")
cat("训练集DBSCAN参数: eps =", eps_train, ", minPts =", minPts_train, "\n\n")

cat("真实标签分布:\n")
print(table(true_labels))

cat("\n数据集分割信息:\n")
cat("训练集样本数量:", nrow(X_train), "\n")
cat("测试集样本数量:", nrow(X_test), "\n")
cat("总样本数量:", nrow(X), "\n")

cat("\n原始数据DBSCAN聚类结果:\n")
print(table(original_clusters))

cat("\n训练集DBSCAN聚类结果:\n")
print(table(train_clusters))

cat("\n测试集聚类结果（与训练集一一对应）:\n")
print(table(test_clusters))

# 计算聚类质量指标
if (require(cluster, quietly = TRUE)) {
  cat("\n=== 聚类质量评估 ===\n")
  
  # 原始数据聚类质量
  if (length(unique(original_clusters)) > 1) {
    silhouette_original <- silhouette(original_clusters, dist(X[, 1:p1]))
    cat("原始数据平均轮廓系数:", mean(silhouette_original[, 3]), "\n")
  }
  
  # 训练集聚类质量
  if (length(unique(train_clusters)) > 1) {
    silhouette_train <- silhouette(train_clusters, dist(X_train[, 1:p1]))
    cat("训练集平均轮廓系数:", mean(silhouette_train[, 3]), "\n")
  }
  
  # 测试集聚类质量
  if (length(unique(test_clusters)) > 1) {
    silhouette_test <- silhouette(test_clusters, dist(X_test[, 1:p1]))
    cat("测试集平均轮廓系数:", mean(silhouette_test[, 3]), "\n")
  }
}

# 计算聚类一致性
cat("\n=== 聚类一致性分析 ===\n")
# 比较原始聚类和训练集聚类的一致性
original_vs_train <- table(original_clusters, train_clusters)
cat("原始聚类 vs 训练集聚类交叉表:\n")
print(original_vs_train)

# 比较原始聚类和测试集聚类的一致性
original_vs_test <- table(original_clusters, test_clusters)
cat("\n原始聚类 vs 测试集聚类交叉表:\n")
print(original_vs_test)

# 比较训练集和测试集聚类的一致性
train_vs_test <- table(train_clusters, test_clusters)
cat("\n训练集聚类 vs 测试集聚类交叉表:\n")
print(train_vs_test)

# 显示图像
print(combined_plot)

cat("\n图像已保存为 DT_DBSCAN_clustering_comparison.png\n")
