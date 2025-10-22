# t_dist_SI_multi与IVS_cluster的Power和FDR比较分析

# 加载必要的包
library(ggplot2)
library(gridExtra)
library(CADET)
library(ClusterR)
library(MASS)
library(gamlss.dist)
library(doSNOW)


# 参数设置（与t_dist_normtest.R保持一致）
n <- 400   # 样本数
p <- 200   # 总特征数
p1 <- 40   # 信号特征数
K <- 2     # 聚类数
n_repetitions <- 100  # 重复次数

# t分布参数
df_t <- 8  # t分布自由度参数
mu_1 <- 0  # 第一类均值
mu_2 <- 4  # 第二类均值
sigma_1 <- 16  # 第一类标准差
sigma_2 <- 1  # 第二类标准差
rho <- 0.2  # 相关系数

# 维度选择（在循环外固定）
# 除了第p1+1维之外的所有维度
available_dims <- setdiff(1:p, p1+1)
# 从可用维度中随机选择p1个维度作为信号维度
signal_dims <- sample(available_dims, size = p1, replace = FALSE)
# 确保第p1+1维是噪声（不在信号维度中）
noise_dims <- setdiff(1:p, signal_dims)

# 构建协方差矩阵
Sigma <- rho^abs(matrix(rep(1:p, each=p), p, p, byrow = F) - matrix(rep(1:p, each=p), p, p, byrow = T))

# 创建方差向量：信号特征使用sigma_1，噪声特征使用sigma_2
variance_vector <- rep(sigma_2, p)  # 默认所有特征使用sigma_2
variance_vector[signal_dims] <- sigma_1  # 信号特征使用sigma_1

Sigma <- sqrt(diag(variance_vector)) %*% Sigma %*% sqrt(diag(variance_vector))

cat("信号维度:", signal_dims, "\n")
cat("噪声维度:", noise_dims, "\n")
cat("信号特征方差 (sigma_1):", sigma_1, "\n")
cat("噪声特征方差 (sigma_2):", sigma_2, "\n\n")

# 从gaussian_signal.R复制的sym_threshold函数
sym_threshold <- function(psi, alpha, option = c('1', '0')){
  opt = as.numeric(match.arg(option))
  sw = sort(unique(c(0, abs(psi))), decreasing = FALSE)
  efdp = sapply(sw, function(t){
    (opt + sum(psi < -t))/max(1, sum(psi > t))
  }) #计算所有的score对应的值
  L = ifelse(sum(efdp <= alpha), sw[min(which(efdp <= alpha))], +Inf) #如果没有一个小于alpha，则记为正无穷，否则取最小的那个
  disc = which(psi > L) #选出显著的对应的位置
  return(list(disc = disc, L = L))
}

# 从gaussian_signal.R复制的IVS_cluster函数
IVS_cluster <- function(X, K = 2){
  n = nrow(X) #样本的数量
  p = ncol(X) #特征的数量
  
  psi_kj = matrix(NA, nrow = K, ncol = p)
  beta_kj = matrix(NA, nrow = K, ncol = p)
  permus = sample(1:n, n, replace = FALSE) #1到n permute的函数
  folds = lapply(1:K, function(k) permus[seq(k, n, by = K)] ) #样本分割的函数
  for (kk in 1:K) {
    foldIdx = folds[[kk]]
    Xtrain = X[-foldIdx, ]#挑选出一份训练集
    center_cluster = ClusterR::KMeans_arma(Xtrain, clusters = 2, seed = 0)
    # mean(apply(Xtrain - center_cluster[ClusterR::predict_KMeans(Xtrain, center_cluster), ], 1, function(X){ sqrt(sum(X^2)) }))
    for (jj in 1:p){
      idxj = sample(foldIdx, ceiling(length(foldIdx)/2), replace = FALSE)#将测试集分成两份
      Xtest0 = X[idxj, ]#挑选出ori
      Xj0 = Xtest0[, jj]#ori的第jj个变量
      Xtest1 = X[setdiff(foldIdx, idxj), ]#挑选出syn
      Xj1 = Xtest1[, jj]#syn的第jj个变量
      psi_kj[kk, jj] = mean((Xj1 - center_cluster[3-ClusterR::predict_KMeans(Xtest1, center_cluster), jj])^2) -
        mean((Xj0 - center_cluster[ClusterR::predict_KMeans(Xtest0, center_cluster), jj])^2)
      beta_kj[kk, jj] = mean((Xj1 - center_cluster[3-ClusterR::predict_KMeans(Xtest1, center_cluster), jj])^2+
                               (Xj0 - center_cluster[3-ClusterR::predict_KMeans(Xtest0, center_cluster), jj])^2)
    }
  }
  beta = colMeans(beta_kj)
  psi = colMeans(psi_kj)
  bpsi = colMeans(psi_kj)*beta
  return(list(psi = psi, bpsi = bpsi))
}

# t_dist_SI_multi函数（从si_correlation_analysis.R复制）
t_dist_SI_multi <- function(X, j){
  # 计算第j维的p值
  feat <- j
  result <- CADET::kmeans_inference_1f(X, k = 2, 1, 2, feat = feat, iso = F, covMat = cov(X), seed = 2021, iter.max = 30)
  return(result$pval)
}

# 主循环：比较两种方法
cat("开始t_dist_SI_multi与IVS_cluster的Power和FDR比较分析...\n")
cat("参数设置: n =", n, ", p =", p, ", p1 =", p1, ", K =", K, ", 重复次数 =", n_repetitions, "\n")
cat("t分布参数: df =", df_t, ", mu1 =", mu_1, ", mu2 =", mu_2, ", rho =", rho, "\n\n")

# 存储结果
results <- data.frame(
  repetition = 1:n_repetitions,
  SI_power = numeric(n_repetitions),
  SI_fdr = numeric(n_repetitions),
  IVS_psi_power = numeric(n_repetitions),
  IVS_psi_fdr = numeric(n_repetitions),
  IVS_bpsi_power = numeric(n_repetitions),
  IVS_bpsi_fdr = numeric(n_repetitions)
)

for (rep in 1:n_repetitions) {
  if (rep %% 10 == 0) {
    cat("完成重复次数:", rep, "/", n_repetitions, "\n")
  }
  
  # 使用与t_dist_normtest.R相同的数据生成方式
  set.seed(rep + 1000)  # 确保每次重复使用不同的随机种子
  
  # 计算协方差矩阵的1/2次方（Cholesky分解）
  Sigma_sqrt <- chol(Sigma)
  
  # 设置样本分配
  npor <- 0.4
  alpha <- 0.1  # 使用0.05作为显著性水平
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
  
  # === 方法1: t_dist_SI_multi ===
  # 对所有特征计算p值
  all_pvalues <- numeric(p)
  for (j in 1:p) {
    all_pvalues[j] <- t_dist_SI_multi(X, j)
  }
  
  # 使用BH算法进行多重检验校正
  bh_adjusted <- p.adjust(all_pvalues, method = "BH")
  rejected_features_SI <- which(bh_adjusted < alpha)
  
  # 计算SI方法的Power和FDR
  # Power = |拒绝集合 ∩ 信号维度| / |信号维度|
  results$SI_power[rep] <- length(intersect(rejected_features_SI, signal_dims)) / length(signal_dims)
  
  # FDR = |拒绝集合 ∩ 噪声维度| / |拒绝集合|
  if (length(rejected_features_SI) > 0) {
    results$SI_fdr[rep] <- length(intersect(rejected_features_SI, noise_dims)) / length(rejected_features_SI)
  } else {
    results$SI_fdr[rep] <- 0  # 如果没有拒绝任何特征，FDR为0
  }
  
  # === 方法2: IVS_cluster ===
  # 调用IVS_cluster函数
  ivs_result <- IVS_cluster(X, K = K)
  
  # 使用sym_threshold处理psi结果（未增强）
  ivs_psi_result <- sym_threshold(ivs_result$psi, alpha, "1")
  rejected_features_IVS_psi <- ivs_psi_result$disc
  
  # 使用sym_threshold处理bpsi结果（增强）
  ivs_bpsi_result <- sym_threshold(ivs_result$bpsi, alpha, "1")
  rejected_features_IVS_bpsi <- ivs_bpsi_result$disc
  
  # 计算IVS_psi方法的Power和FDR
  results$IVS_psi_power[rep] <- length(intersect(rejected_features_IVS_psi, signal_dims)) / length(signal_dims)
  if (length(rejected_features_IVS_psi) > 0) {
    results$IVS_psi_fdr[rep] <- length(intersect(rejected_features_IVS_psi, noise_dims)) / length(rejected_features_IVS_psi)
  } else {
    results$IVS_psi_fdr[rep] <- 0
  }
  
  # 计算IVS_bpsi方法的Power和FDR
  results$IVS_bpsi_power[rep] <- length(intersect(rejected_features_IVS_bpsi, signal_dims)) / length(signal_dims)
  if (length(rejected_features_IVS_bpsi) > 0) {
    results$IVS_bpsi_fdr[rep] <- length(intersect(rejected_features_IVS_bpsi, noise_dims)) / length(rejected_features_IVS_bpsi)
  } else {
    results$IVS_bpsi_fdr[rep] <- 0
  }
}

cat("完成所有重复次数!\n\n")

# 输出结果摘要
cat("=== Power和FDR比较结果 ===\n")
cat("t_dist_SI_multi方法:\n")
cat("  平均Power:", round(mean(results$SI_power), 4), "±", round(sd(results$SI_power), 4), "\n")
cat("  平均FDR:", round(mean(results$SI_fdr), 4), "±", round(sd(results$SI_fdr), 4), "\n\n")

cat("IVS_cluster (psi, 未增强)方法:\n")
cat("  平均Power:", round(mean(results$IVS_psi_power), 4), "±", round(sd(results$IVS_psi_power), 4), "\n")
cat("  平均FDR:", round(mean(results$IVS_psi_fdr), 4), "±", round(sd(results$IVS_psi_fdr), 4), "\n\n")

cat("IVS_cluster (bpsi, 增强)方法:\n")
cat("  平均Power:", round(mean(results$IVS_bpsi_power), 4), "±", round(sd(results$IVS_bpsi_power), 4), "\n")
cat("  平均FDR:", round(mean(results$IVS_bpsi_fdr), 4), "±", round(sd(results$IVS_bpsi_fdr), 4), "\n\n")

# 创建可视化
library(reshape2)

# 准备数据用于可视化
power_data <- data.frame(
  Method = rep(c("SI", "IVS_psi", "IVS_bpsi"), each = n_repetitions),
  Power = c(results$SI_power, results$IVS_psi_power, results$IVS_bpsi_power)
)

fdr_data <- data.frame(
  Method = rep(c("SI", "IVS_psi", "IVS_bpsi"), each = n_repetitions),
  FDR = c(results$SI_fdr, results$IVS_psi_fdr, results$IVS_bpsi_fdr)
)

# Power比较箱线图
p_power <- ggplot(power_data, aes(x = Method, y = Power, fill = Method)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(title = "Power比较: t_dist_SI_multi vs IVS_cluster",
       x = "方法", y = "Power") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  scale_fill_manual(values = c("SI" = "lightblue", "IVS_psi" = "lightgreen", "IVS_bpsi" = "lightcoral"))

# FDR比较箱线图
p_fdr <- ggplot(fdr_data, aes(x = Method, y = FDR, fill = Method)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(title = "FDR比较: t_dist_SI_multi vs IVS_cluster",
       x = "方法", y = "FDR") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  scale_fill_manual(values = c("SI" = "lightblue", "IVS_psi" = "lightgreen", "IVS_bpsi" = "lightcoral"))

# 合并图形
combined_plot <- grid.arrange(p_power, p_fdr, ncol = 2)

# 保存图形
ggsave("t_dist_SI_vs_IVS_comparison.png", combined_plot, width = 12, height = 6, dpi = 300)
cat("比较图已保存为 t_dist_SI_vs_IVS_comparison.png\n")

# 创建详细的统计摘要表
summary_stats <- data.frame(
  Method = c("SI", "IVS_psi", "IVS_bpsi"),
  Power_Mean = c(mean(results$SI_power), mean(results$IVS_psi_power), mean(results$IVS_bpsi_power)),
  Power_SD = c(sd(results$SI_power), sd(results$IVS_psi_power), sd(results$IVS_bpsi_power)),
  Power_Median = c(median(results$SI_power), median(results$IVS_psi_power), median(results$IVS_bpsi_power)),
  FDR_Mean = c(mean(results$SI_fdr), mean(results$IVS_psi_fdr), mean(results$IVS_bpsi_fdr)),
  FDR_SD = c(sd(results$SI_fdr), sd(results$IVS_psi_fdr), sd(results$IVS_bpsi_fdr)),
  FDR_Median = c(median(results$SI_fdr), median(results$IVS_psi_fdr), median(results$IVS_bpsi_fdr))
)

print(summary_stats)

# 保存结果到CSV文件
write.csv(results, "t_dist_SI_vs_IVS_detailed_results.csv", row.names = FALSE)
write.csv(summary_stats, "t_dist_SI_vs_IVS_summary_stats.csv", row.names = FALSE)

cat("详细结果已保存为 t_dist_SI_vs_IVS_detailed_results.csv\n")
cat("统计摘要已保存为 t_dist_SI_vs_IVS_summary_stats.csv\n")

# 进行统计检验比较
cat("\n=== 统计检验比较 ===\n")

# Power比较
cat("Power比较 (t检验):\n")
power_test_1 <- t.test(results$SI_power, results$IVS_psi_power)
cat("SI vs IVS_psi: p-value =", round(power_test_1$p.value, 4), "\n")

power_test_2 <- t.test(results$SI_power, results$IVS_bpsi_power)
cat("SI vs IVS_bpsi: p-value =", round(power_test_2$p.value, 4), "\n")

power_test_3 <- t.test(results$IVS_psi_power, results$IVS_bpsi_power)
cat("IVS_psi vs IVS_bpsi: p-value =", round(power_test_3$p.value, 4), "\n\n")

# FDR比较
cat("FDR比较 (t检验):\n")
fdr_test_1 <- t.test(results$SI_fdr, results$IVS_psi_fdr)
cat("SI vs IVS_psi: p-value =", round(fdr_test_1$p.value, 4), "\n")

fdr_test_2 <- t.test(results$SI_fdr, results$IVS_bpsi_fdr)
cat("SI vs IVS_bpsi: p-value =", round(fdr_test_2$p.value, 4), "\n")

fdr_test_3 <- t.test(results$IVS_psi_fdr, results$IVS_bpsi_fdr)
cat("IVS_psi vs IVS_bpsi: p-value =", round(fdr_test_3$p.value, 4), "\n")

cat("\n分析完成!\n")
