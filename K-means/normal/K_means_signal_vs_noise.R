library(ClusterR)
library(gridExtra)
library(ggplot2)
library(cluster)

# 数据生成参数
n <- 40  # 样本数量
p <- 5  # 总特征数
p1 <- 2  # 前p1个特征有信号
mu_1 <- 0
mu_2 <- 1

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

# 使用K-means进行聚类
k <- 2  # 聚类数量
set.seed(123)
kmeans_result <- ClusterR::KMeans_arma(X, clusters = k, n_iter = 100, seed_mode = "random_subset", 
                                       verbose = FALSE, CENTROIDS = NULL, seed = 123)
kmeans_labels <- ClusterR::predict_KMeans(X, kmeans_result)

# 创建数据框用于绘图
plot_data <- data.frame(
  x_signal = X[, 1],  # 第1个特征（有信号）
  y_signal = X[, 2],  # 第2个特征（有信号）
  x_noise = X[, 1],   # 第1个特征（用于噪声对比）
  y_noise = X[, p1 + 1],  # 第p1+1个特征（噪声）
  true_label = as.factor(true_labels),
  kmeans_label = as.factor(kmeans_labels)
)

# 创建真实标签的散点图（有信号特征）
plot_true_signal <- ggplot(plot_data, aes(x = x_signal, y = y_signal, color = true_label)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(
    title = "真实标签分布（有信号特征）",
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

# 创建K-means聚类标签的散点图（有信号特征）
plot_kmeans_signal <- ggplot(plot_data, aes(x = x_signal, y = y_signal, color = kmeans_label)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(
    title = "K-means聚类结果（有信号特征）",
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
  scale_color_manual(values = c("1" = "#FF6B6B", "2" = "#4ECDC4", "3" = "#45B7D1", "4" = "#96CEB4"))

# 创建真实标签的散点图（噪声特征）
plot_true_noise <- ggplot(plot_data, aes(x = x_noise, y = y_noise, color = true_label)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(
    title = "真实标签分布（噪声特征）",
    x = "Signal 1",
    y = "Noise 1",
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

# 创建K-means聚类标签的散点图（噪声特征）
plot_kmeans_noise <- ggplot(plot_data, aes(x = x_noise, y = y_noise, color = kmeans_label)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(
    title = "K-means聚类结果（噪声特征）",
    x = "Signal 1",
    y = "Noise 1",
    color = "聚类标签"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  scale_color_manual(values = c("1" = "#FF6B6B", "2" = "#4ECDC4", "3" = "#45B7D1", "4" = "#96CEB4"))

# 组合四个图
combined_plot <- grid.arrange(plot_true_signal, plot_kmeans_signal, 
                             plot_true_noise, plot_kmeans_noise, 
                             ncol = 2, nrow = 2)

# 保存图像
ggsave("K_means_signal_vs_noise_comparison.png", combined_plot, width = 12, height = 12, dpi = 300)

# 打印聚类结果统计
cat("数据生成参数:\n")
cat("样本数量 n =", n, "\n")
cat("总特征数 p =", p, "\n")
cat("信号特征数 p1 =", p1, "\n")
cat("真实标签分布:\n")
print(table(true_labels))
cat("\nK-means聚类结果:\n")
print(table(kmeans_labels))

# 计算聚类质量指标
if (length(unique(kmeans_labels)) > 1) {
  silhouette_score <- silhouette(kmeans_labels, dist(X[, 1:p1]))
  cat("\n平均轮廓系数:", mean(silhouette_score[, 3]), "\n")
}

