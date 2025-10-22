library(dbscan)
library(gridExtra)
library(ggplot2)


# 数据生成参数
n <- 250  # 样本数量
p <- 100  # 总特征数
p1 <- 20  # 前p1个特征有信号
mu_1 <- 0
mu_2 <- 3

# 生成数据
generate_data <- function(n, p, p1, mu_1, mu_2) {
  # 生成真实标签
  true_labels <- sample(c(1, 2), n, replace = TRUE)
  
  # 初始化数据矩阵
  X <- matrix(0, nrow = n, ncol = p)
  
  # 为前p1个特征生成有信号的数据
  for (i in 1:p1) {
    X[true_labels == 1, i] <- rnorm(sum(true_labels == 1), mean = mu_1, sd = 1)
    X[true_labels == 2, i] <- rnorm(sum(true_labels == 2), mean = mu_2, sd = 1)
  }
  
  # 为第p1+1和p1+2个特征生成噪声（无信号）
  # 为剩余特征生成噪声
  for (i in (p1 + 1):p) {
    X[, i] <- rnorm(n, mean = 0, sd = 1)
  }
  
  return(list(X = X, true_labels = true_labels))
}

# 生成数据
data_result <- generate_data(n, p, p1, mu_1, mu_2)
X <- data_result$X
true_labels <- data_result$true_labels

# 使用DBSCAN进行聚类
# 选择eps和minPts参数
eps <- 13
minPts <- 5

# 对所有特征进行DBSCAN聚类
dbscan_result <- dbscan(X, eps = eps, minPts = minPts)
dbscan_labels <- dbscan_result$cluster

# 创建数据框用于绘图
plot_data <- data.frame(
  x = X[, 1],  # 第1个特征作为x坐标
  y = X[, 2],  # 第2个特征作为y坐标
  true_label = as.factor(true_labels),
  dbscan_label = as.factor(dbscan_labels)
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
  scale_color_manual(values = c("1" = "#FF6B6B", "2" = "#4ECDC4"))

# 创建DBSCAN聚类标签的散点图
plot_dbscan <- ggplot(plot_data, aes(x = x, y = y, color = dbscan_label)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(
    title = "DBSCAN聚类结果",
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

# 组合两个图
combined_plot <- grid.arrange(plot_true, plot_dbscan, ncol = 2)

# 保存图像
ggsave("DBSCAN_clustering_comparison.png", combined_plot, width = 12, height = 6, dpi = 300)

# 打印聚类结果统计
cat("数据生成参数:\n")
cat("样本数量 n =", n, "\n")
cat("总特征数 p =", p, "\n")
cat("信号特征数 p1 =", p1, "\n")
cat("真实标签分布:\n")
print(table(true_labels))
cat("\nDBSCAN聚类结果:\n")
print(table(dbscan_labels))
cat("\nDBSCAN参数: eps =", eps, ", minPts =", minPts, "\n")

# 计算聚类质量指标
# 计算轮廓系数（需要cluster包）
if (require(cluster, quietly = TRUE)) {
  if (length(unique(dbscan_labels)) > 1) {
    silhouette_score <- silhouette(dbscan_labels, dist(X[, 1:p1]))
    cat("\n平均轮廓系数:", mean(silhouette_score[, 3]), "\n")
  }
}

# 显示图像
print(combined_plot)
