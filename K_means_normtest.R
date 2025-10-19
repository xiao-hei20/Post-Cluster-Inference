rm(list = ls());
library(class)
library(doSNOW)
library(gamlss.dist)
library(MASS)
library(ClusterR)
library(CADET)
library(tidyr)
library(Seurat)
library(scDesign3)
library(SingleCellExperiment)
library(ggplot2)
library(gridExtra)

norm_CTT <- function(X, p1){
  # 保存当前随机种子状态
  old_seed <- .Random.seed
  
  center_cluster = ClusterR::KMeans_arma(X, clusters = 2)
  
  # 恢复随机种子状态
  .Random.seed <<- old_seed
  
  clu<-ClusterR::predict_KMeans(X, center_cluster)
  ifel<-function(a,b,c){
    if(a){
      b
    }else{
      c
    }
  }
  # 只计算第p1+1列
  jj = p1 + 1
  numerator<-ifel(length(which(clu==1))==1,X[which(clu==1),jj],mean(X[which(clu==1),jj]))-ifel(length(which(clu==2))==1,X[which(clu==2),jj],mean(X[which(clu==2),jj]))
  denominator<-sqrt(ifel(length(which(clu==1))==1,0,var(X[which(clu==1),jj]))/length(which(clu==1))+ifel(length(which(clu==2))==1,0,var(X[which(clu==2),jj]))/length(which(clu==2)))
  t<-numerator/(sqrt(2)*denominator)
  return(t)
}


# 优化版本：只计算第p1+1维的p值
norm_SI <- function(X, p1){
  # 只计算第p1+1维的p值
  feat <- p1 + 1
  result <- CADET::kmeans_inference_1f(X, k = 2, 1, 2, feat = feat, iso = F, covMat = cov(X), seed = 2021, iter.max = 30)
  return(result$pval)
}

# 基于norm_split的优化版本：只计算第p1+1维
norm_DT <- function(X, p1, epi = 0.5){
  # 保存当前随机种子状态
  old_seed <- .Random.seed
  
  center_cluster = ClusterR::KMeans_arma(X, clusters = 2)
  
  # 恢复随机种子状态
  .Random.seed <<- old_seed
  
  clu <- ClusterR::predict_KMeans(X, center_cluster)
  sigma <- (var(X[clu==1,]) + var(X[clu==2,])) / 2
  X_test <- mvrnorm(n = nrow(X), mu = rep(0, ncol(X)), Sigma = epi * (1 - epi) * sigma) + epi * X
  X_train <- X - X_test
  
  center_cluster = ClusterR::KMeans_arma(X_train, clusters = 2)
  clu <- ClusterR::predict_KMeans(X_test, center_cluster)
  
  
  # 只计算第p1+1维
  jj <- p1 + 1
  
  # 检查聚类结果，处理空聚类的情况
  clu1_indices <- which(clu==1)
  clu2_indices <- which(clu==2)
  
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

norm_PCI <- function(X, p1, K = 2){
  n = nrow(X) #样本的数量
  p = ncol(X) #特征的数量
  
  if (p1 + 1 > p) {
    stop("p1 + 1 exceeds the number of columns in X")
  }
  
  # 存储所有fold的平方差，不取均值
  Xj1_squared_diff = c() # 存储Xj1的平方差
  Xj0_squared_diff = c() # 存储Xj0的平方差
  permus = sample(1:n, n, replace = FALSE) #1到n permute的函数
  folds = lapply(1:K, function(k) permus[seq(k, n, by = K)] ) #样本分割的函数
  for (kk in 1:K) {
    foldIdx = folds[[kk]]
    Xtrain = X[-foldIdx, ]#挑选出一份训练集
    
    # 保存当前随机种子状态
    old_seed <- .Random.seed
    
    center_cluster = ClusterR::KMeans_arma(Xtrain, clusters = 2)
    
    # 恢复随机种子状态
    .Random.seed <<- old_seed
    
    # 只计算第p1+1列
    jj = p1 + 1
    idxj = sample(foldIdx, ceiling(length(foldIdx)/2), replace = FALSE)#将测试集分成两份
    Xtest0 = X[idxj, ]#挑选出ori
    Xj0 = Xtest0[, jj]#ori的第jj个变量
    Xtest1 = X[setdiff(foldIdx, idxj), ]#挑选出syn
    Xj1 = Xtest1[, jj]#syn的第jj个变量
    
    # 计算平方差但不取均值，存储到向量中
    Xj1_diff = (Xj1 - center_cluster[3-ClusterR::predict_KMeans(Xtest1, center_cluster), jj])^2
    Xj0_diff = (Xj0 - center_cluster[ClusterR::predict_KMeans(Xtest0, center_cluster), jj])^2
    
    # 将当前fold的结果添加到总向量中
    Xj1_squared_diff = c(Xj1_squared_diff, Xj1_diff)
    Xj0_squared_diff = c(Xj0_squared_diff, Xj0_diff)
  }
  
  # 计算分子：Xj1_squared_diff的均值减去Xj0_squared_diff的均值
  numerator = sqrt(n/4)*(mean(Xj1_squared_diff) - mean(Xj0_squared_diff))
  
  # 计算分母：两个向量方差的均值
  denominator = sqrt((var(Xj1_squared_diff) + var(Xj0_squared_diff)) / 2)
  
  # 返回比值
  return(numerator / denominator)
}


res_list<-list()
seq_tim <- 1
p1_seq     <- seq(from = 20, to = 20, length.out = seq_tim)
p_seq      <- seq(from = 100, to = 100, length.out = seq_tim)
n_seq      <- seq(from = 400, to = 400, length.out = seq_tim)
mu_1_seq   <- seq(from = 0, to = 0, length.out = seq_tim)
mu_2_seq   <- seq(from = 0.5, to = 0.5,length.out = seq_tim)
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
  # X_LU<-matrix(rnbinom(n = n0*p1, size = r1, prob = prob1), nrow = n0, ncol = p1)
  # X_LD<-matrix(rnbinom(n = n1*p1, size = r2, prob = prob2), nrow = n1, ncol = p1)
  # X_R<-matrix(rnbinom(n = n*(p-p1), size = r3, prob = prob3), nrow = n, ncol = p-p1)
  # X<-cbind(rbind(X_LU,X_LD),X_R)
  # 设置重复次数
  n_repetitions <- 100
  
  # 调用norm_PCI、norm_CTT、norm_SI和norm_DT函数n_repetitions次，每次重新生成X，存储结果
  result_PCI <- numeric(n_repetitions)
  result_CTT <- numeric(n_repetitions)
  result_SI <- numeric(n_repetitions)
  result_DT <- numeric(n_repetitions)
  for(l in 1:n_repetitions) {
  # 重新生成X矩阵
  set.seed(l+ran)
  X1 = t(t(2*matrix(rnorm(n1*p), n1, p) %*% chol(Sigma))+ center1)
  X0 = t(t(2*matrix(rnorm(n0*p), n0, p) %*% chol(Sigma))+ center0) 
  X = rbind(X1, X0)
  result_PCI[l] <- norm_PCI(X, p1, K=4)
  result_CTT[l] <- norm_CTT(X, p1)
  
  # 计算norm_SI的p值并转换为标准正态值
  p_value_si <- norm_SI(X, p1)
  result_SI[l] <- qnorm(p_value_si)
  
  # 计算norm_DT的t统计量
  result_DT[l] <- norm_DT(X, p1)
  }
  
  # 验证result是否为正态分布
  # 使用ggplot2绘制图形
  
  # 准备数据
  df <- data.frame(
    values = c(result_PCI, result_CTT, result_SI, result_DT),
    method = rep(c("norm_PCI", "norm_CTT", "norm_SI", "norm_DT"), each = n_repetitions)
  )
  
  # 1. 直方图 + 核密度估计 + 正态分布曲线
  p1_PCI <- ggplot(df[df$method == "norm_PCI", ], aes(x = values)) +
    geom_histogram(aes(y = ..density..), bins = 20, 
                   fill = "lightblue", color = "black", alpha = 0.7) +
    geom_density(color = "red", size = 1) +
    stat_function(fun = dnorm, 
                  args = list(mean = mean(result_PCI), sd = sd(result_PCI)),
                  color = "blue", linetype = "dashed", size = 1) +
    labs(title = "Distribution of norm_PCI results",
         x = "norm_PCI values", 
         y = "Density") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  p1_CTT <- ggplot(df[df$method == "norm_CTT", ], aes(x = values)) +
    geom_histogram(aes(y = ..density..), bins = 20, 
                   fill = "lightgreen", color = "black", alpha = 0.7) +
    geom_density(color = "red", size = 1) +
    stat_function(fun = dnorm, 
                  args = list(mean = mean(result_CTT), sd = sd(result_CTT)),
                  color = "blue", linetype = "dashed", size = 1) +
    labs(title = "Distribution of norm_CTT results",
         x = "norm_CTT values", 
         y = "Density") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  p1_SI <- ggplot(df[df$method == "norm_SI", ], aes(x = values)) +
    geom_histogram(aes(y = ..density..), bins = 20, 
                   fill = "lightcoral", color = "black", alpha = 0.7) +
    geom_density(color = "red", size = 1) +
    stat_function(fun = dnorm, 
                  args = list(mean = mean(result_SI), sd = sd(result_SI)),
                  color = "blue", linetype = "dashed", size = 1) +
    labs(title = "Distribution of norm_SI results",
         x = "norm_SI values", 
         y = "Density") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  p1_DT <- ggplot(df[df$method == "norm_DT", ], aes(x = values)) +
    geom_histogram(aes(y = ..density..), bins = 20, 
                   fill = "lightgreen", color = "black", alpha = 0.7) +
    geom_density(color = "red", size = 1) +
    stat_function(fun = dnorm, 
                  args = list(mean = mean(result_DT), sd = sd(result_DT)),
                  color = "blue", linetype = "dashed", size = 1) +
    labs(title = "Distribution of norm_DT results",
         x = "norm_DT values", 
         y = "Density") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  # 2. Q-Q图
  p2_PCI <- ggplot(df[df$method == "norm_PCI", ], aes(sample = values)) +
    stat_qq() +
    stat_qq_line(color = "red", size = 1) +
    labs(title = "Q-Q Plot: norm_PCI vs Normal Distribution",
         x = "Theoretical Quantiles",
         y = "Sample Quantiles") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  p2_CTT <- ggplot(df[df$method == "norm_CTT", ], aes(sample = values)) +
    stat_qq() +
    stat_qq_line(color = "red", size = 1) +
    labs(title = "Q-Q Plot: norm_CTT vs Normal Distribution",
         x = "Theoretical Quantiles",
         y = "Sample Quantiles") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  p2_SI <- ggplot(df[df$method == "norm_SI", ], aes(sample = values)) +
    stat_qq() +
    stat_qq_line(color = "red", size = 1) +
    labs(title = "Q-Q Plot: norm_SI vs Normal Distribution",
         x = "Theoretical Quantiles",
         y = "Sample Quantiles") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  p2_DT <- ggplot(df[df$method == "norm_DT", ], aes(sample = values)) +
    stat_qq() +
    stat_qq_line(color = "red", size = 1) +
    labs(title = "Q-Q Plot: norm_DT vs Normal Distribution",
         x = "Theoretical Quantiles",
         y = "Sample Quantiles") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  # 并排显示图形 - 2行4列布局
  combined_plot <- grid.arrange(p1_PCI, p1_CTT, p1_SI, p1_DT,
                                p2_PCI, p2_CTT, p2_SI, p2_DT, 
                                ncol = 4, nrow = 2)
  
  # 保存图形到文件
  ggsave("normality_comparison.png", combined_plot, width = 12, height = 10, dpi = 300)
  cat("图形已保存为 normality_comparison.png\n")
  
  # 3. 正态性检验
  cat("=== norm_PCI 正态性检验 ===\n")
  cat("Shapiro-Wilk normality test:\n")
  shapiro_test_PCI <- shapiro.test(result_PCI)
  print(shapiro_test_PCI)
  
  cat("\nKolmogorov-Smirnov test:\n")
  ks_test_PCI <- ks.test(result_PCI, "pnorm", mean = mean(result_PCI), sd = sd(result_PCI))
  print(ks_test_PCI)
  
  cat("\n=== norm_CTT 正态性检验 ===\n")
  cat("Shapiro-Wilk normality test:\n")
  shapiro_test_CTT <- shapiro.test(result_CTT)
  print(shapiro_test_CTT)
  
  cat("\nKolmogorov-Smirnov test:\n")
  ks_test_CTT <- ks.test(result_CTT, "pnorm", mean = mean(result_CTT), sd = sd(result_CTT))
  print(ks_test_CTT)
  
  cat("\n=== norm_SI 正态性检验 ===\n")
  cat("Shapiro-Wilk normality test:\n")
  shapiro_test_SI <- shapiro.test(result_SI)
  print(shapiro_test_SI)
  
  cat("\nKolmogorov-Smirnov test:\n")
  ks_test_SI <- ks.test(result_SI, "pnorm", mean = mean(result_SI), sd = sd(result_SI))
  print(ks_test_SI)
  
  cat("\n=== norm_DT 正态性检验 ===\n")
  cat("Shapiro-Wilk normality test:\n")
  shapiro_test_DT <- shapiro.test(result_DT)
  print(shapiro_test_DT)
  
  cat("\nKolmogorov-Smirnov test:\n")
  ks_test_DT <- ks.test(result_DT, "pnorm", mean = mean(result_DT), sd = sd(result_DT))
  print(ks_test_DT)
  
  # 4. 基本统计信息
  cat("\n=== norm_PCI 基本统计信息 ===\n")
  cat("Mean:", mean(result_PCI), "\n")
  cat("SD:", sd(result_PCI), "\n")
  cat("Skewness:", moments::skewness(result_PCI), "\n")
  cat("Kurtosis:", moments::kurtosis(result_PCI), "\n")
  
  cat("\n=== norm_CTT 基本统计信息 ===\n")
  cat("Mean:", mean(result_CTT), "\n")
  cat("SD:", sd(result_CTT), "\n")
  cat("Skewness:", moments::skewness(result_CTT), "\n")
  cat("Kurtosis:", moments::kurtosis(result_CTT), "\n")
  
  cat("\n=== norm_SI 基本统计信息 ===\n")
  cat("Mean:", mean(result_SI), "\n")
  cat("SD:", sd(result_SI), "\n")
  cat("Skewness:", moments::skewness(result_SI), "\n")
  cat("Kurtosis:", moments::kurtosis(result_SI), "\n")
  
  cat("\n=== norm_DT 基本统计信息 ===\n")
  cat("Mean:", mean(result_DT), "\n")
  cat("SD:", sd(result_DT), "\n")
  cat("Skewness:", moments::skewness(result_DT), "\n")
  cat("Kurtosis:", moments::kurtosis(result_DT), "\n")
  
  # 注意：ggplot2不需要恢复图形参数
}







