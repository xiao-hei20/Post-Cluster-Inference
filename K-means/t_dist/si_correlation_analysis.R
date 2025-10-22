# SI方法相关系数分析脚本
# 固定设置，只使用SI方法，记录所有null特征的p值并计算相关系数矩阵

# 加载必要的包
library(ggplot2)
library(gridExtra)
library(CADET)
library(ClusterR)
library(MASS)
library(gamlss.dist)

# 设置随机种子
set.seed(123)

# 参数设置（与t_dist_normtest.R保持一致）
n <- 100   # 样本数
p <- 200   # 总特征数
p1 <- 40   # 信号特征数
K <- 2     # 聚类数
n_repetitions <- 100  # 重复次数

# t分布参数
df_t <- 3  # t分布自由度参数
mu_1 <- 0  # 第一类均值
mu_2 <- 2.9  # 第二类均值
sigma_1 <- 1  # 第一类标准差
sigma_2 <- 1  # 第二类标准差
rho <- 0.2  # 相关系数

# SI方法函数（使用t_dist_SI_multi）
t_dist_SI_multi <- function(X, j){
  # 计算第j维的p值
  feat <- j
  result <- CADET::kmeans_inference_1f(X, k = 2, 1, 2, feat = feat, iso = F, covMat = cov(X), seed = 2021, iter.max = 30)
  return(result$pval)
}

# 构建协方差矩阵（与t_dist_normtest.R相同）
Sigma <- rho^abs(matrix(rep(1:p, each=p), p, p, byrow = F) - matrix(rep(1:p, each=p), p, p, byrow = T))
Sigma <- sqrt(diag(c(rep(sigma_1, p1), rep(sigma_2, p-p1)))) %*% Sigma %*% sqrt(diag(c(rep(sigma_1, p1), rep(sigma_2, p-p1))))

# 维度选择（在循环外固定）
# 除了第p1+1维之外的所有维度
available_dims <- setdiff(1:p, p1+1)
# 从可用维度中随机选择p1个维度作为信号维度
signal_dims <- sample(available_dims, size = p1, replace = FALSE)
# 确保第p1+1维是噪声（不在信号维度中）
noise_dims <- setdiff(1:p, signal_dims)

cat("信号维度:", signal_dims, "\n")
cat("噪声维度:", noise_dims, "\n\n")

# 主循环：记录所有null特征的p值
cat("开始SI方法相关系数分析...\n")
cat("参数设置: n =", n, ", p =", p, ", p1 =", p1, ", K =", K, ", 重复次数 =", n_repetitions, "\n")
cat("t分布参数: df =", df_t, ", mu1 =", mu_1, ", mu2 =", mu_2, ", rho =", rho, "\n\n")

# 存储所有null特征的p值矩阵
null_pvalues_matrix <- matrix(0, nrow = n_repetitions, ncol = length(noise_dims))

# 存储Power和FDR结果
power_results <- numeric(n_repetitions)
fdr_results <- numeric(n_repetitions)

for (rep in 1:n_repetitions) {
  if (rep %% 10 == 0) {
    cat("完成重复次数:", rep, "/", n_repetitions, "\n")
  }
  
  # 使用与t_dist_normtest.R相同的数据生成方式
  set.seed(rep + 1000)  # 确保每次重复使用不同的随机种子
  
  # 计算协方差矩阵的1/2次方（Cholesky分解）
  Sigma_sqrt <- chol(Sigma)
  
  # 设置样本分配
  npor <- 0.2
  alpha <- 0.2
  n1 <- ceiling(n * npor)  # 第一类样本数
  n0 <- n - n1            # 第二类样本数
  
  # 设置中心点
  center0 <- rep(mu_1, p)
  center1 <- rep(mu_1, p)
  center1[signal_dims] <- mu_2  # 只在选定的信号维度上添加信号
  
  # 生成独立的t分布数据
  X1 <- matrix(rt(n1 * p, df = df_t), n1, p)
  X0 <- matrix(rt(n0 * p, df = df_t), n0, p)
  
  # 应用协方差结构：乘以协方差矩阵的1/2次方
  X1 <- X1 %*% t(Sigma_sqrt)
  X0 <- X0 %*% t(Sigma_sqrt)
  
  # 添加均值偏移（信号构造）
  X1 <- t(t(X1) + center1)
  X0 <- t(t(X0) + center0)
  
  # 组合数据
  X <- rbind(X1, X0)
  
  # 对所有特征计算p值
  all_pvalues <- numeric(p)
  for (j in 1:p) {
    all_pvalues[j] <- t_dist_SI_multi(X, j)
  }
  
  # 使用BH算法进行多重检验校正
  alpha <- 0.05
  bh_adjusted <- p.adjust(all_pvalues, method = "BH")
  rejected_features <- which(bh_adjusted < alpha)
  
  # 计算Power和FDR
  # Power = |拒绝集合 ∩ 信号维度| / |信号维度|
  power_results[rep] <- length(intersect(rejected_features, signal_dims)) / length(signal_dims)
  
  # FDR = |拒绝集合 ∩ 噪声维度| / |拒绝集合|
  if (length(rejected_features) > 0) {
    fdr_results[rep] <- length(intersect(rejected_features, noise_dims)) / length(rejected_features)
  } else {
    fdr_results[rep] <- 0  # 如果没有拒绝任何特征，FDR为0
  }
  
  # 对每个null特征（噪声特征）计算p值（用于相关性分析）
  null_pvalues <- all_pvalues[noise_dims]
  null_pvalues_matrix[rep, ] <- null_pvalues
}

cat("完成所有重复次数!\n\n")

# 输出Power和FDR结果
cat("=== Power和FDR分析结果 ===\n")
cat("平均Power:", round(mean(power_results), 4), "\n")
cat("Power标准差:", round(sd(power_results), 4), "\n")
cat("Power范围:", round(range(power_results), 4), "\n")
cat("平均FDR:", round(mean(fdr_results), 4), "\n")
cat("FDR标准差:", round(sd(fdr_results), 4), "\n")
cat("FDR范围:", round(range(fdr_results), 4), "\n\n")

# 计算相关系数矩阵
cat("计算相关系数矩阵...\n")
correlation_matrix <- cor(null_pvalues_matrix)

# 输出相关系数矩阵的基本统计信息
cat("相关系数矩阵维度:", dim(correlation_matrix), "\n")
cat("相关系数范围:", round(range(correlation_matrix), 4), "\n")
cat("平均相关系数:", round(mean(correlation_matrix[upper.tri(correlation_matrix)]), 4), "\n")
cat("相关系数标准差:", round(sd(correlation_matrix[upper.tri(correlation_matrix)]), 4), "\n\n")

# 可视化相关系数矩阵
library(reshape2)
library(viridis)

# 准备数据
cor_data <- melt(correlation_matrix)
colnames(cor_data) <- c("Var1", "Var2", "Correlation")

# 创建热图
p_cor_heatmap <- ggplot(cor_data, aes(x = Var1, y = Var2, fill = Correlation)) +
  geom_tile() +
  scale_fill_viridis(name = "Correlation", limits = c(-1, 1)) +
  labs(title = "SI Method: Correlation Matrix of Null Feature P-values",
       x = "Null Feature Index", 
       y = "Null Feature Index") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))

# 保存热图
ggsave("si_correlation_heatmap.png", p_cor_heatmap, width = 10, height = 8, dpi = 300)
cat("相关系数热图已保存为 si_correlation_heatmap.png\n")

# 创建相关系数分布直方图
cor_values <- correlation_matrix[upper.tri(correlation_matrix)]
cor_df <- data.frame(Correlation = cor_values)

p_cor_hist <- ggplot(cor_df, aes(x = Correlation)) +
  geom_histogram(aes(y = ..density..), bins = 30, 
                 fill = "lightblue", color = "black", alpha = 0.7) +
  geom_density(color = "red", size = 1) +
  labs(title = "SI Method: Distribution of Correlation Coefficients",
       x = "Correlation Coefficient", 
       y = "Density") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# 保存直方图
ggsave("si_correlation_distribution.png", p_cor_hist, width = 8, height = 6, dpi = 300)
cat("相关系数分布图已保存为 si_correlation_distribution.png\n")

# 输出详细的统计信息
cat("\n=== SI方法相关系数分析结果 ===\n")
cat("相关系数统计摘要:\n")
cat("  最小值:", round(min(cor_values), 4), "\n")
cat("  25%分位数:", round(quantile(cor_values, 0.25), 4), "\n")
cat("  中位数:", round(median(cor_values), 4), "\n")
cat("  75%分位数:", round(quantile(cor_values, 0.75), 4), "\n")
cat("  最大值:", round(max(cor_values), 4), "\n")
cat("  均值:", round(mean(cor_values), 4), "\n")
cat("  标准差:", round(sd(cor_values), 4), "\n\n")

# 检查高相关性的特征对
high_cor_threshold <- 0.5
high_cor_pairs <- which(abs(correlation_matrix) > high_cor_threshold & 
                       upper.tri(correlation_matrix), arr.ind = TRUE)

if (nrow(high_cor_pairs) > 0) {
  cat("高相关性特征对 (|相关系数| >", high_cor_threshold, "):\n")
  for (i in 1:nrow(high_cor_pairs)) {
    row_idx <- high_cor_pairs[i, 1]
    col_idx <- high_cor_pairs[i, 2]
    cor_val <- correlation_matrix[row_idx, col_idx]
    cat("  特征", row_idx, "与特征", col_idx, ": 相关系数 =", round(cor_val, 4), "\n")
  }
} else {
  cat("没有发现高相关性的特征对 (|相关系数| >", high_cor_threshold, ")\n")
}

# 保存相关系数矩阵到文件
write.csv(correlation_matrix, "si_correlation_matrix.csv", row.names = FALSE)
cat("\n相关系数矩阵已保存为 si_correlation_matrix.csv\n")

# 创建Power和FDR的可视化
png("si_power_fdr_analysis.png", width = 1200, height = 600)
par(mfrow = c(1, 2))

# Power分布
hist(power_results, breaks = 20, main = "SI方法Power分布", 
     xlab = "Power", ylab = "频数", col = "lightblue", border = "black")
abline(v = mean(power_results), col = "red", lwd = 2, lty = 2)
legend("topright", legend = paste("均值 =", round(mean(power_results), 4)), 
       col = "red", lty = 2, lwd = 2)

# FDR分布
hist(fdr_results, breaks = 20, main = "SI方法FDR分布", 
     xlab = "FDR", ylab = "频数", col = "lightcoral", border = "black")
abline(v = mean(fdr_results), col = "red", lwd = 2, lty = 2)
legend("topright", legend = paste("均值 =", round(mean(fdr_results), 4)), 
       col = "red", lty = 2, lwd = 2)

dev.off()
cat("Power和FDR分析图已保存为 si_power_fdr_analysis.png\n")

# 保存Power和FDR结果到文件
results_summary <- data.frame(
  Repetition = 1:n_repetitions,
  Power = power_results,
  FDR = fdr_results
)
write.csv(results_summary, "si_power_fdr_results.csv", row.names = FALSE)
cat("Power和FDR详细结果已保存为 si_power_fdr_results.csv\n")

# 保存p值矩阵到文件
write.csv(null_pvalues_matrix, "si_null_pvalues_matrix.csv", row.names = FALSE)
cat("null特征p值矩阵已保存为 si_null_pvalues_matrix.csv\n")

cat("\n分析完成!\n")
