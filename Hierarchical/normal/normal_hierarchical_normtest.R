rm(list = ls());
library(class)
library(doSNOW)
library(gamlss.dist)
library(MASS)
library(tidyr)
library(cluster)
library(ggplot2)
library(gridExtra)
library(moments)

# 层次聚类版本的CTT函数
norm_CTT_hierarchical <- function(X, p1){
  # 保存当前随机种子状态
  old_seed <- .Random.seed
  
  # 使用DIANA分裂聚类替代K-means
  diana_result <- diana(X)
  clu <- cutree(diana_result, k = 2)
  
  # 恢复随机种子状态
  .Random.seed <<- old_seed
  
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

# 层次聚类版本的DT函数
norm_DT_hierarchical <- function(X, p1, epi = 0.5){
  # 保存当前随机种子状态
  old_seed <- .Random.seed
  
  # 使用DIANA分裂聚类
  diana_result <- diana(X)
  clu <- cutree(diana_result, k = 2)
  
  # 计算协方差矩阵
  sigma <- (var(X[clu==1,]) + var(X[clu==2,])) / 2
  
  # 生成测试数据
  X_test <- mvrnorm(n = nrow(X), mu = rep(0, ncol(X)), Sigma = epi * (1 - epi) * sigma) + epi * X
  X_train <- X - X_test
  
  # 对训练数据进行DIANA分裂聚类
  diana_result_train <- diana(X_train)
  clu_train <- cutree(diana_result_train, k = 2)
  
  # 恢复随机种子状态
  .Random.seed <<- old_seed
  
  # 只计算第p1+1维
  jj <- p1 + 1
  
  # 检查聚类结果，处理空聚类的情况
  clu1_indices <- which(clu_train==1)
  clu2_indices <- which(clu_train==2)
  
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

# 层次聚类版本的PCI函数
norm_PCI_hierarchical <- function(X, p1, K = 2){
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
    
    # 使用DIANA分裂聚类替代K-means
    diana_result <- diana(Xtrain)
    cluster_centers <- matrix(0, nrow = 2, ncol = ncol(Xtrain))
    
    # 计算聚类中心
    clu_train <- cutree(diana_result, k = 2)
    for(i in 1:2) {
      if(sum(clu_train == i) > 0) {
        cluster_centers[i, ] <- colMeans(Xtrain[clu_train == i, , drop = FALSE])
      }
    }
    
    # 恢复随机种子状态
    .Random.seed <<- old_seed
    
    # 只计算第p1+1列
    jj = p1 + 1
    idxj = sample(foldIdx, ceiling(length(foldIdx)/2), replace = FALSE)#将测试集分成两份
    Xtest0 = X[idxj, ]#挑选出ori
    Xj0 = Xtest0[, jj]#ori的第jj个变量
    Xtest1 = X[setdiff(foldIdx, idxj), ]#挑选出syn
    Xj1 = Xtest1[, jj]#syn的第jj个变量
    
    # 预测测试数据的聚类标签
    # 计算到聚类中心的距离
    dist_to_center1 <- sqrt(rowSums((Xtest0 - matrix(cluster_centers[1,], nrow=nrow(Xtest0), ncol=ncol(Xtest0), byrow=TRUE))^2))
    dist_to_center2 <- sqrt(rowSums((Xtest0 - matrix(cluster_centers[2,], nrow=nrow(Xtest0), ncol=ncol(Xtest0), byrow=TRUE))^2))
    clu_test0 <- ifelse(dist_to_center1 < dist_to_center2, 1, 2)
    
    dist_to_center1_syn <- sqrt(rowSums((Xtest1 - matrix(cluster_centers[1,], nrow=nrow(Xtest1), ncol=ncol(Xtest1), byrow=TRUE))^2))
    dist_to_center2_syn <- sqrt(rowSums((Xtest1 - matrix(cluster_centers[2,], nrow=nrow(Xtest1), ncol=ncol(Xtest1), byrow=TRUE))^2))
    clu_test1 <- ifelse(dist_to_center1_syn < dist_to_center2_syn, 1, 2)
    
    # 计算平方差但不取均值，存储到向量中
    Xj1_diff = (Xj1 - cluster_centers[3-clu_test1, jj])^2
    Xj0_diff = (Xj0 - cluster_centers[clu_test0, jj])^2
    
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
  # 设置重复次数
  n_repetitions <- 100
  
  # 调用层次聚类版本的函数n_repetitions次，每次重新生成X，存储结果
  result_PCI <- numeric(n_repetitions)
  result_CTT <- numeric(n_repetitions)
  result_DT <- numeric(n_repetitions)
  for(l in 1:n_repetitions) {
  # 重新生成X矩阵 - 使用正态分布
  set.seed(l+ran)
  X1 = t(t(2*matrix(rnorm(n1*p), n1, p) %*% chol(Sigma))+ center1)
  X0 = t(t(2*matrix(rnorm(n0*p), n0, p) %*% chol(Sigma))+ center0) 
  X = rbind(X1, X0)
  result_PCI[l] <- norm_PCI_hierarchical(X, p1, K=4)
  result_CTT[l] <- norm_CTT_hierarchical(X, p1)
  
  # 计算norm_DT的t统计量
  result_DT[l] <- norm_DT_hierarchical(X, p1)
  }
  
  # 验证result是否为正态分布
  # 使用ggplot2绘制图形
  
  # 准备数据
  df <- data.frame(
    values = c(result_PCI, result_CTT, result_DT),
    method = rep(c("norm_PCI_hierarchical", "norm_CTT_hierarchical", "norm_DT_hierarchical"), each = n_repetitions)
  )
  
  # 1. 直方图 + 核密度估计 + 正态分布曲线
  p1_PCI <- ggplot(df[df$method == "norm_PCI_hierarchical", ], aes(x = values)) +
    geom_histogram(aes(y = ..density..), bins = 20, 
                   fill = "lightblue", color = "black", alpha = 0.7) +
    geom_density(color = "red", size = 1) +
    stat_function(fun = dnorm, 
                  args = list(mean = mean(result_PCI), sd = sd(result_PCI)),
                  color = "blue", linetype = "dashed", size = 1) +
    labs(title = "Distribution of norm_PCI_hierarchical results",
         x = "norm_PCI_hierarchical values", 
         y = "Density") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  p1_CTT <- ggplot(df[df$method == "norm_CTT_hierarchical", ], aes(x = values)) +
    geom_histogram(aes(y = ..density..), bins = 20, 
                   fill = "lightgreen", color = "black", alpha = 0.7) +
    geom_density(color = "red", size = 1) +
    stat_function(fun = dnorm, 
                  args = list(mean = mean(result_CTT), sd = sd(result_CTT)),
                  color = "blue", linetype = "dashed", size = 1) +
    labs(title = "Distribution of norm_CTT_hierarchical results",
         x = "norm_CTT_hierarchical values", 
         y = "Density") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  p1_DT <- ggplot(df[df$method == "norm_DT_hierarchical", ], aes(x = values)) +
    geom_histogram(aes(y = ..density..), bins = 20, 
                   fill = "lightcoral", color = "black", alpha = 0.7) +
    geom_density(color = "red", size = 1) +
    stat_function(fun = dnorm, 
                  args = list(mean = mean(result_DT), sd = sd(result_DT)),
                  color = "blue", linetype = "dashed", size = 1) +
    labs(title = "Distribution of norm_DT_hierarchical results",
         x = "norm_DT_hierarchical values", 
         y = "Density") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  # 2. Q-Q图
  p2_PCI <- ggplot(df[df$method == "norm_PCI_hierarchical", ], aes(sample = values)) +
    stat_qq() +
    stat_qq_line(color = "red", size = 1) +
    labs(title = "Q-Q Plot: norm_PCI_hierarchical vs Normal Distribution",
         x = "Theoretical Quantiles",
         y = "Sample Quantiles") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  p2_CTT <- ggplot(df[df$method == "norm_CTT_hierarchical", ], aes(sample = values)) +
    stat_qq() +
    stat_qq_line(color = "red", size = 1) +
    labs(title = "Q-Q Plot: norm_CTT_hierarchical vs Normal Distribution",
         x = "Theoretical Quantiles",
         y = "Sample Quantiles") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  p2_DT <- ggplot(df[df$method == "norm_DT_hierarchical", ], aes(sample = values)) +
    stat_qq() +
    stat_qq_line(color = "red", size = 1) +
    labs(title = "Q-Q Plot: norm_DT_hierarchical vs Normal Distribution",
         x = "Theoretical Quantiles",
         y = "Sample Quantiles") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  # 并排显示图形 - 2行3列布局
  combined_plot <- grid.arrange(p1_PCI, p1_CTT, p1_DT,
                                p2_PCI, p2_CTT, p2_DT, 
                                ncol = 3, nrow = 2)
  
  # 保存图形到文件
  ggsave("Hierarchical/normal/normality_comparison.png", combined_plot, width = 12, height = 10, dpi = 300)
  cat("图形已保存为 normality_comparison.png\n")

}

