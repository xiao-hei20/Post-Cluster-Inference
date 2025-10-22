rm(list = ls());
library(class)
library(doSNOW)
library(gamlss.dist)
library(MASS)
library(dbscan)
library(CADET)
library(tidyr)
library(Seurat)
library(scDesign3)
library(SingleCellExperiment)
library(ggplot2)
library(gridExtra)
library(moments)

# DBSCAN版本的norm_CTT函数
norm_CTT_dbscan <- function(X, p1, eps = 15, minPts = 3){
  # 保存当前随机种子状态
  old_seed <- .Random.seed
  
  # 使用DBSCAN进行聚类
  dbscan_result <- dbscan::dbscan(X, eps = eps, minPts = minPts)
  clu <- dbscan_result$cluster
  
  # 恢复随机种子状态
  .Random.seed <<- old_seed
  
  # 检查是否找到了聚类
  unique_clusters <- unique(clu[clu != 0])  # 排除噪声点(标记为0)
  if(length(unique_clusters) < 2) {
    # 如果聚类数不足，返回0
    return(0)
  }
  
  # 计算每个聚类的样本数量
  cluster_sizes <- sapply(unique_clusters, function(c) sum(clu == c))
  
  # 按样本数量排序，选择前两个最大的聚类
  sorted_indices <- order(cluster_sizes, decreasing = TRUE)
  cluster1 <- unique_clusters[sorted_indices[1]]
  cluster2 <- unique_clusters[sorted_indices[2]]
  
  ifel<-function(a,b,c){
    if(a){
      b
    }else{
      c
    }
  }
  # 只计算第p1+1列
  jj = p1 + 1
  
  # 获取两个聚类的数据点
  clu1_indices <- which(clu == cluster1)
  clu2_indices <- which(clu == cluster2)
  
  # 计算前两类样本数量之和（当前未使用，但保留以备将来扩展）
  # n_top2 <- length(clu1_indices) + length(clu2_indices)
  
  # 计算统计量
  numerator<-ifel(length(clu1_indices)==1,X[clu1_indices,jj],mean(X[clu1_indices,jj]))-ifel(length(clu2_indices)==1,X[clu2_indices,jj],mean(X[clu2_indices,jj]))
  denominator<-sqrt(ifel(length(clu1_indices)==1,0,var(X[clu1_indices,jj]))/length(clu1_indices)+ifel(length(clu2_indices)==1,0,var(X[clu2_indices,jj]))/length(clu2_indices))
  t<-numerator/(sqrt(2)*denominator)
  return(t)
}

# DBSCAN版本的norm_DT函数
norm_DT_dbscan <- function(X, p1, epi = 0.5, eps_original = 15, minPts_original = 3, eps_train = 15, minPts_train = 3){
  # 保存当前随机种子状态
  old_seed <- .Random.seed
  
  # 使用DBSCAN进行聚类（原始数据参数）
  dbscan_result <- dbscan::dbscan(X, eps = eps_original, minPts = minPts_original)
  clu <- dbscan_result$cluster
  
  # 恢复随机种子状态
  .Random.seed <<- old_seed
  
  # 检查是否找到了聚类
  unique_clusters <- unique(clu[clu != 0])  # 排除噪声点
  if(length(unique_clusters) < 2) {
    return(0)
  }
  
  # 计算每个聚类的样本数量
  cluster_sizes <- sapply(unique_clusters, function(c) sum(clu == c))
  
  # 按样本数量排序，选择前两个最大的聚类
  sorted_indices <- order(cluster_sizes, decreasing = TRUE)
  cluster1 <- unique_clusters[sorted_indices[1]]
  cluster2 <- unique_clusters[sorted_indices[2]]
  
  # 计算协方差矩阵
  clu1_indices <- which(clu == cluster1)
  clu2_indices <- which(clu == cluster2)
  
  if(length(clu1_indices) == 0 || length(clu2_indices) == 0) {
    return(0)
  }
  
  sigma <- (var(X[clu1_indices,]) + var(X[clu2_indices,])) / 2
  X_test <- mvrnorm(n = nrow(X), mu = rep(0, ncol(X)), Sigma = epi * (1 - epi) * sigma) + epi * X
  X_train <- X - X_test
  
  # 在训练集上重新进行DBSCAN聚类（使用训练集参数）
  dbscan_train <- dbscan::dbscan(X_train, eps = eps_train, minPts = minPts_train)
  clu_test <- dbscan_train$cluster
  
  # 检查测试集聚类结果
  unique_test_clusters <- unique(clu_test[clu_test != 0])
  if(length(unique_test_clusters) < 2) {
    return(0)
  }
  
  # 计算测试集每个聚类的样本数量
  test_cluster_sizes <- sapply(unique_test_clusters, function(c) sum(clu_test == c))
  
  # 按样本数量排序，选择前两个最大的聚类
  test_sorted_indices <- order(test_cluster_sizes, decreasing = TRUE)
  cluster1_test <- unique_test_clusters[test_sorted_indices[1]]
  cluster2_test <- unique_test_clusters[test_sorted_indices[2]]
  
  # 只计算第p1+1维
  jj <- p1 + 1
  
  # 检查聚类结果，处理空聚类的情况
  clu1_indices <- which(clu_test == cluster1_test)
  clu2_indices <- which(clu_test == cluster2_test)
  
  # 计算第一类的均值和方差
  if(length(clu1_indices) == 0) {
    mean1 <- 0
    var1 <- 0
  } else if(length(clu1_indices) == 1) {
    mean1 <- X_test[clu1_indices, jj]
    var1 <- 0
  } else {
    mean1 <- mean(X_test[clu1_indices, jj])
    var1 <- var(X_test[clu1_indices, jj])
  }
  
  # 计算第二类的均值和方差
  if(length(clu2_indices) == 0) {
    mean2 <- 0
    var2 <- 0
  } else if(length(clu2_indices) == 1) {
    mean2 <- X_test[clu2_indices, jj]
    var2 <- 0
  } else {
    mean2 <- mean(X_test[clu2_indices, jj])
    var2 <- var(X_test[clu2_indices, jj])
  }
  
  # 计算t统计量
  numerator <- mean1 - mean2
  denominator <- sqrt(var1 / max(length(clu1_indices), 1) + var2 / max(length(clu2_indices), 1))
  
  t <- numerator / denominator
  return(t)
}

# DBSCAN版本的norm_PCI函数
norm_PCI_dbscan <- function(X, p1, K = 2, eps = 15, minPts = 3){
  n = nrow(X) #样本的数量
  p = ncol(X) #特征的数量
  
  if (p1 + 1 > p) {
    stop("p1 + 1 exceeds the number of columns in X")
  }
  
  # 存储所有fold的平方差，不取均值
  Xj1_squared_diff = c() # 存储Xj1的平方差
  Xj0_squared_diff = c() # 存储Xj0的平方差
  n_top2_total = 0 # 存储所有fold中前两类样本数量的总和
  permus = sample(1:n, n, replace = FALSE) #1到n permute的函数
  folds = lapply(1:K, function(k) permus[seq(k, n, by = K)] ) #样本分割的函数
  for (kk in 1:K) {
    foldIdx = folds[[kk]]
    Xtrain = X[-foldIdx, ]#挑选出一份训练集
    
    # 保存当前随机种子状态
    old_seed <- .Random.seed
    
    # 使用DBSCAN进行聚类
    dbscan_result <- dbscan::dbscan(Xtrain, eps = eps, minPts = minPts)
    clu_train <- dbscan_result$cluster
    
    # 恢复随机种子状态
    .Random.seed <<- old_seed
    
    # 检查聚类结果
    unique_clusters <- unique(clu_train[clu_train != 0])
    if(length(unique_clusters) < 2) {
      # 如果聚类数不足，跳过这个fold
      next
    }
    
    # 计算每个聚类的样本数量
    cluster_sizes <- sapply(unique_clusters, function(c) sum(clu_train == c))
    
    # 按样本数量排序，选择前两个最大的聚类
    sorted_indices <- order(cluster_sizes, decreasing = TRUE)
    cluster1 <- unique_clusters[sorted_indices[1]]
    cluster2 <- unique_clusters[sorted_indices[2]]
    
    # 计算当前fold中前两类的样本数量
    n_top2_current <- cluster_sizes[sorted_indices[1]] + cluster_sizes[sorted_indices[2]]
    n_top2_total <- n_top2_total + n_top2_current
    
    # 只计算第p1+1列
    jj = p1 + 1
    idxj = sample(foldIdx, ceiling(length(foldIdx)/2), replace = FALSE)#将测试集分成两份
    Xtest0 = X[idxj, ]#挑选出ori
    Xj0 = Xtest0[, jj]#ori的第jj个变量
    Xtest1 = X[setdiff(foldIdx, idxj), ]#挑选出syn
    Xj1 = Xtest1[, jj]#syn的第jj个变量
    
    # 预测测试集的聚类标签
    # 对于DBSCAN，我们需要重新计算距离并分配聚类
    # 这里使用简单的最近邻方法
    clu0 <- rep(0, nrow(Xtest0))
    clu1 <- rep(0, nrow(Xtest1))
    
    # 计算到聚类中心的距离
    center1 <- colMeans(Xtrain[clu_train == cluster1, , drop = FALSE])
    center2 <- colMeans(Xtrain[clu_train == cluster2, , drop = FALSE])
    
    # 分配聚类标签
    for(i in seq_len(nrow(Xtest0))) {
      dist1 <- sum((Xtest0[i, ] - center1)^2)
      dist2 <- sum((Xtest0[i, ] - center2)^2)
      clu0[i] <- ifelse(dist1 < dist2, cluster1, cluster2)
    }
    
    for(i in seq_len(nrow(Xtest1))) {
      dist1 <- sum((Xtest1[i, ] - center1)^2)
      dist2 <- sum((Xtest1[i, ] - center2)^2)
      clu1[i] <- ifelse(dist1 < dist2, cluster1, cluster2)
    }
    
    # 计算平方差但不取均值，存储到向量中
    Xj1_diff = (Xj1 - ifelse(clu1 == cluster1, center1[jj], center2[jj]))^2
    Xj0_diff = (Xj0 - ifelse(clu0 == cluster1, center1[jj], center2[jj]))^2
    
    # 将当前fold的结果添加到总向量中
    Xj1_squared_diff = c(Xj1_squared_diff, Xj1_diff)
    Xj0_squared_diff = c(Xj0_squared_diff, Xj0_diff)
  }
  
  # 检查是否有有效数据
  if(length(Xj1_squared_diff) == 0 || length(Xj0_squared_diff) == 0) {
    return(0)
  }
  
  # 计算前两类样本数量之和（用于标准化）
  # 使用所有fold中前两类样本数量的平均值
  n_top2_avg <- n_top2_total / K
  
  # 计算分子：Xj1_squared_diff的均值减去Xj0_squared_diff的均值
  numerator = sqrt(n_top2_avg/4)*(mean(Xj1_squared_diff) - mean(Xj0_squared_diff))
  
  # 计算分母：两个向量方差的均值
  denominator = sqrt((var(Xj1_squared_diff) + var(Xj0_squared_diff)) / 2)
  
  # 返回比值
  return(numerator / denominator)
}

res_list<-list()
seq_tim <- 1
p1_seq     <- seq(from = 20, to = 20, length.out = seq_tim)
p_seq      <- seq(from = 100, to = 200, length.out = seq_tim)
n_seq      <- seq(from = 150, to = 150, length.out = seq_tim)
mu_1_seq   <- seq(from = 0, to = 0, length.out = seq_tim)
mu_2_seq   <- seq(from = 3, to = 3, length.out = seq_tim)
sigma_1_seq <- seq(from = 1, to = 1, length.out = seq_tim)
sigma_2_seq <- seq(from = 1, to = 1, length.out = seq_tim)
rho_seq    <- seq(from = 0, to = 0, length.out = seq_tim)
param_list <- list(
  p1_seq = p1_seq,
  p_seq = p_seq,
  n_seq = n_seq,
  mu_1_seq = mu_1_seq,
  mu_2_seq = mu_2_seq,
  sigma_1_seq = sigma_1_seq,
  sigma_2_seq = sigma_2_seq,
  rho_seq = rho_seq
)
str(param_list)
res_list<-list()

# DBSCAN参数设置
eps_original <- 13      # 原始数据邻域半径
minPts_original <- 5   # 原始数据最小点数
eps_train <- 10        # 训练集邻域半径
minPts_train <- 5      # 训练集最小点数

for (i in 1:seq_tim) {
  ran<-sample(1:100000,1)
  param_names <- names(param_list)
  base_names <- sub("_seq$", "", param_names)
  for (j in seq_along(param_names)) {
    assign(base_names[j], param_list[[param_names[j]]][i], envir = environment())
  }
  Sigma=rho^abs(matrix(rep(1:p,each=p),p,p,byrow = F)-matrix(rep(1:p,each=p),p,p,byrow = T))
  Sigma=sqrt(diag(c(rep(sigma_1,p1),rep(sigma_2,p-p1))))%*%Sigma%*%sqrt(diag(c(rep(sigma_1,p1),rep(sigma_2,p-p1))))
  npor=0.2;alpha=0.2
  n1 = ceiling(n*npor) #number
  n0 = n - n1
  loc = 1:p1
  center0 = rep(mu_1, p);
  center1 = rep(mu_1, p); center1[loc] = mu_2
  n1 = ceiling(n*npor)
  n0 = n - n1
  Z = c(rep(1, n1), rep(0, n0))
  
  # 设置重复次数
  n_repetitions <- 100
  
  # 调用DBSCAN版本的函数n_repetitions次，每次重新生成X，存储结果
  result_PCI <- numeric(n_repetitions)
  result_CTT <- numeric(n_repetitions)
  result_DT <- numeric(n_repetitions)
  
  for(l in 1:n_repetitions) {
    # 重新生成X矩阵
    set.seed(l+ran)
    X1 = t(t(matrix(rnorm(n1*p), n1, p) %*% chol(Sigma))+ center1)
    X0 = t(t(matrix(rnorm(n0*p), n0, p) %*% chol(Sigma))+ center0) 
    X = rbind(X1, X0)
    
    result_PCI[l] <- norm_PCI_dbscan(X, p1, K=5, eps = eps_original, minPts = minPts_original)
    result_CTT[l] <- norm_CTT_dbscan(X, p1, eps = eps_original, minPts = minPts_original)
    result_DT[l] <- norm_DT_dbscan(X, p1, epi = 0.5, eps_original = eps_original, minPts_original = minPts_original, eps_train = eps_train, minPts_train = minPts_train)
  }
  
  # 验证result是否为正态分布
  # 使用ggplot2绘制图形
  
  # 准备数据
  df <- data.frame(
    values = c(result_PCI, result_CTT, result_DT),
    method = rep(c("norm_PCI_dbscan", "norm_CTT_dbscan", "norm_DT_dbscan"), each = n_repetitions)
  )
  
  # 1. 直方图 + 核密度估计 + 正态分布曲线
  p1_PCI <- ggplot(df[df$method == "norm_PCI_dbscan", ], aes(x = values)) +
    geom_histogram(aes(y = ..density..), bins = 20, 
                   fill = "lightblue", color = "black", alpha = 0.7) +
    geom_density(color = "red", size = 1) +
    stat_function(fun = dnorm, 
                  args = list(mean = mean(result_PCI), sd = sd(result_PCI)),
                  color = "blue", linetype = "dashed", size = 1) +
    labs(title = "Distribution of norm_PCI_dbscan results",
         x = "norm_PCI_dbscan values", 
         y = "Density") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  p1_CTT <- ggplot(df[df$method == "norm_CTT_dbscan", ], aes(x = values)) +
    geom_histogram(aes(y = ..density..), bins = 20, 
                   fill = "lightgreen", color = "black", alpha = 0.7) +
    geom_density(color = "red", size = 1) +
    stat_function(fun = dnorm, 
                  args = list(mean = mean(result_CTT), sd = sd(result_CTT)),
                  color = "blue", linetype = "dashed", size = 1) +
    labs(title = "Distribution of norm_CTT_dbscan results",
         x = "norm_CTT_dbscan values", 
         y = "Density") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  p1_DT <- ggplot(df[df$method == "norm_DT_dbscan", ], aes(x = values)) +
    geom_histogram(aes(y = ..density..), bins = 20, 
                   fill = "lightcoral", color = "black", alpha = 0.7) +
    geom_density(color = "red", size = 1) +
    stat_function(fun = dnorm, 
                  args = list(mean = mean(result_DT), sd = sd(result_DT)),
                  color = "blue", linetype = "dashed", size = 1) +
    labs(title = "Distribution of norm_DT_dbscan results",
         x = "norm_DT_dbscan values", 
         y = "Density") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  # 2. Q-Q图
  p2_PCI <- ggplot(df[df$method == "norm_PCI_dbscan", ], aes(sample = values)) +
    stat_qq() +
    stat_qq_line(color = "red", size = 1) +
    labs(title = "Q-Q Plot: norm_PCI_dbscan vs Normal Distribution",
         x = "Theoretical Quantiles",
         y = "Sample Quantiles") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  p2_CTT <- ggplot(df[df$method == "norm_CTT_dbscan", ], aes(sample = values)) +
    stat_qq() +
    stat_qq_line(color = "red", size = 1) +
    labs(title = "Q-Q Plot: norm_CTT_dbscan vs Normal Distribution",
         x = "Theoretical Quantiles",
         y = "Sample Quantiles") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  p2_DT <- ggplot(df[df$method == "norm_DT_dbscan", ], aes(sample = values)) +
    stat_qq() +
    stat_qq_line(color = "red", size = 1) +
    labs(title = "Q-Q Plot: norm_DT_dbscan vs Normal Distribution",
         x = "Theoretical Quantiles",
         y = "Sample Quantiles") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  # 并排显示图形 - 2行3列布局
  combined_plot <- grid.arrange(p1_PCI, p1_CTT, p1_DT,
                                p2_PCI, p2_CTT, p2_DT, 
                                ncol = 3, nrow = 2)
  
  # 保存图形到文件
  ggsave("dbscan_normality_comparison.png", combined_plot, width = 15, height = 10, dpi = 300)
  cat("图形已保存为 dbscan_normality_comparison.png\n")
  
  # 3. 正态性检验
  cat("=== norm_PCI_dbscan 正态性检验 ===\n")
  cat("Shapiro-Wilk normality test:\n")
  shapiro_test_PCI <- shapiro.test(result_PCI)
  print(shapiro_test_PCI)
  
  cat("\nKolmogorov-Smirnov test:\n")
  ks_test_PCI <- ks.test(result_PCI, "pnorm", mean = mean(result_PCI), sd = sd(result_PCI))
  print(ks_test_PCI)
  
  cat("\n=== norm_CTT_dbscan 正态性检验 ===\n")
  cat("Shapiro-Wilk normality test:\n")
  shapiro_test_CTT <- shapiro.test(result_CTT)
  print(shapiro_test_CTT)
  
  cat("\nKolmogorov-Smirnov test:\n")
  ks_test_CTT <- ks.test(result_CTT, "pnorm", mean = mean(result_CTT), sd = sd(result_CTT))
  print(ks_test_CTT)
  
  cat("\n=== norm_DT_dbscan 正态性检验 ===\n")
  cat("Shapiro-Wilk normality test:\n")
  shapiro_test_DT <- shapiro.test(result_DT)
  print(shapiro_test_DT)
  
  cat("\nKolmogorov-Smirnov test:\n")
  ks_test_DT <- ks.test(result_DT, "pnorm", mean = mean(result_DT), sd = sd(result_DT))
  print(ks_test_DT)
  
  # 4. 基本统计信息
  cat("\n=== norm_PCI_dbscan 基本统计信息 ===\n")
  cat("Mean:", mean(result_PCI), "\n")
  cat("SD:", sd(result_PCI), "\n")
  cat("Skewness:", moments::skewness(result_PCI), "\n")
  cat("Kurtosis:", moments::kurtosis(result_PCI), "\n")
  
  cat("\n=== norm_CTT_dbscan 基本统计信息 ===\n")
  cat("Mean:", mean(result_CTT), "\n")
  cat("SD:", sd(result_CTT), "\n")
  cat("Skewness:", moments::skewness(result_CTT), "\n")
  cat("Kurtosis:", moments::kurtosis(result_CTT), "\n")
  
  cat("\n=== norm_DT_dbscan 基本统计信息 ===\n")
  cat("Mean:", mean(result_DT), "\n")
  cat("SD:", sd(result_DT), "\n")
  cat("Skewness:", moments::skewness(result_DT), "\n")
  cat("Kurtosis:", moments::kurtosis(result_DT), "\n")
  
  # 5. DBSCAN参数信息
  cat("\n=== DBSCAN 参数设置 ===\n")
  cat("原始数据参数: eps =", eps_original, ", minPts =", minPts_original, "\n")
  cat("训练集参数: eps =", eps_train, ", minPts =", minPts_train, "\n")
}
