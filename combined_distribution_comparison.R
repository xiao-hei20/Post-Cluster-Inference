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
library(cluster)
library(ggplot2)
library(gridExtra)
library(moments)

# K-means版本的函数（适用于所有分布）
CTT_km <- function(X, p1){
  old_seed <- .Random.seed
  center_cluster = ClusterR::KMeans_arma(X, clusters = 2)
  .Random.seed <<- old_seed
  clu <- ClusterR::predict_KMeans(X, center_cluster)
  
  ifel <- function(a,b,c) if(a) b else c
  jj = p1 + 1
  numerator <- ifel(length(which(clu==1))==1, X[which(clu==1),jj], mean(X[which(clu==1),jj])) - 
    ifel(length(which(clu==2))==1, X[which(clu==2),jj], mean(X[which(clu==2),jj]))
  denominator <- sqrt(ifel(length(which(clu==1))==1, 0, var(X[which(clu==1),jj]))/length(which(clu==1)) +
                      ifel(length(which(clu==2))==1, 0, var(X[which(clu==2),jj]))/length(which(clu==2)))
  t <- numerator/(sqrt(2)*denominator)
  return(t)
}

DT_km <- function(X, p1, epi = 0.5){
  old_seed <- .Random.seed
  center_cluster = ClusterR::KMeans_arma(X, clusters = 2)
  .Random.seed <<- old_seed
  clu <- ClusterR::predict_KMeans(X, center_cluster)
  sigma <- (var(X[clu==1,]) + var(X[clu==2,])) / 2
  X_test <- mvrnorm(n = nrow(X), mu = rep(0, ncol(X)), Sigma = epi * (1 - epi) * sigma) + epi * X
  X_train <- X - X_test
  center_cluster = ClusterR::KMeans_arma(X_train, clusters = 2)
  clu <- ClusterR::predict_KMeans(X_train, center_cluster)
  
  jj <- p1 + 1
  clu1_indices <- which(clu==1)
  clu2_indices <- which(clu==2)
  
  if(length(clu1_indices) == 0) { mean1 <- 0; var1 <- 0
  } else if(length(clu1_indices) == 1) { mean1 <- X_test[clu1_indices, jj]; var1 <- 0
  } else { mean1 <- mean(X_test[clu1_indices, jj]); var1 <- var(X_test[clu1_indices, jj]) }
  
  if(length(clu2_indices) == 0) { mean2 <- 0; var2 <- 0
  } else if(length(clu2_indices) == 1) { mean2 <- X_test[clu2_indices, jj]; var2 <- 0
  } else { mean2 <- mean(X_test[clu2_indices, jj]); var2 <- var(X_test[clu2_indices, jj]) }
  
  numerator <- mean1 - mean2
  denominator <- sqrt(var1 / max(length(clu1_indices), 1) + var2 / max(length(clu2_indices), 1))
  return(numerator / denominator)
}

PCI_km <- function(X, p1, K = 2){
  n = nrow(X); p = ncol(X)
  if (p1 + 1 > p) stop("p1 + 1 exceeds the number of columns in X")
  
  Xj1_squared_diff = c(); Xj0_squared_diff = c()
  permus = sample(1:n, n, replace = FALSE)
  folds = lapply(1:K, function(k) permus[seq(k, n, by = K)])
  
  for (kk in 1:K) {
    foldIdx = folds[[kk]]
    Xtrain = X[-foldIdx, ]
    old_seed <- .Random.seed
    center_cluster = ClusterR::KMeans_arma(Xtrain, clusters = 2)
    .Random.seed <<- old_seed
    
    jj = p1 + 1
    idxj = sample(foldIdx, ceiling(length(foldIdx)/2), replace = FALSE)
    Xtest0 = X[idxj, ]; Xj0 = Xtest0[, jj]
    Xtest1 = X[setdiff(foldIdx, idxj), ]; Xj1 = Xtest1[, jj]
    
    Xj1_diff = (Xj1 - center_cluster[3-ClusterR::predict_KMeans(Xtest1, center_cluster), jj])^2
    Xj0_diff = (Xj0 - center_cluster[ClusterR::predict_KMeans(Xtest0, center_cluster), jj])^2
    
    Xj1_squared_diff = c(Xj1_squared_diff, Xj1_diff)
    Xj0_squared_diff = c(Xj0_squared_diff, Xj0_diff)
  }
  
  numerator = sqrt(n/4)*(mean(Xj1_squared_diff) - mean(Xj0_squared_diff))
  denominator = sqrt((var(Xj1_squared_diff) + var(Xj0_squared_diff)) / 2)
  return(numerator / denominator)
}

SI_km <- function(X, p1){
  feat <- p1 + 1
  result <- CADET::kmeans_inference_1f(X, k = 2, 1, 2, feat = feat, iso = F, covMat = cov(X), seed = 2021, iter.max = 30)
  return(result$pval)
}

# 层次聚类版本的函数（适用于所有分布）
CTT_h <- function(X, p1){
  old_seed <- .Random.seed
  diana_result <- diana(X)
  clu <- cutree(diana_result, k = 2)
  .Random.seed <<- old_seed
  
  ifel <- function(a,b,c) if(a) b else c
  jj = p1 + 1
  numerator <- ifel(length(which(clu==1))==1, X[which(clu==1),jj], mean(X[which(clu==1),jj])) - 
    ifel(length(which(clu==2))==1, X[which(clu==2),jj], mean(X[which(clu==2),jj]))
  denominator <- sqrt(ifel(length(which(clu==1))==1, 0, var(X[which(clu==1),jj]))/length(which(clu==1)) +
                      ifel(length(which(clu==2))==1, 0, var(X[which(clu==2),jj]))/length(which(clu==2)))
  t <- numerator/(sqrt(2)*denominator)
  return(t)
}

DT_h <- function(X, p1, epi = 0.5){
  old_seed <- .Random.seed
  diana_result <- diana(X)
  clu <- cutree(diana_result, k = 2)
  sigma <- (var(X[clu==1,]) + var(X[clu==2,])) / 2
  X_test <- mvrnorm(n = nrow(X), mu = rep(0, ncol(X)), Sigma = epi * (1 - epi) * sigma) + epi * X
  X_train <- X - X_test
  diana_result_train <- diana(X_train)
  clu_train <- cutree(diana_result_train, k = 2)
  .Random.seed <<- old_seed
  
  jj <- p1 + 1
  clu1_indices <- which(clu_train==1)
  clu2_indices <- which(clu_train==2)
  
  if(length(clu1_indices) == 0) { mean1 <- 0; var1 <- 0
  } else if(length(clu1_indices) == 1) { mean1 <- X_test[clu1_indices, jj]; var1 <- 0
  } else { mean1 <- mean(X_test[clu1_indices, jj]); var1 <- var(X_test[clu1_indices, jj]) }
  
  if(length(clu2_indices) == 0) { mean2 <- 0; var2 <- 0
  } else if(length(clu2_indices) == 1) { mean2 <- X_test[clu2_indices, jj]; var2 <- 0
  } else { mean2 <- mean(X_test[clu2_indices, jj]); var2 <- var(X_test[clu2_indices, jj]) }
  
  numerator <- mean1 - mean2
  denominator <- sqrt(var1 / max(length(clu1_indices), 1) + var2 / max(length(clu2_indices), 1))
  return(numerator / denominator)
}

PCI_h <- function(X, p1, K = 2){
  n = nrow(X); p = ncol(X)
  if (p1 + 1 > p) stop("p1 + 1 exceeds the number of columns in X")
  
  Xj1_squared_diff = c(); Xj0_squared_diff = c()
  permus = sample(1:n, n, replace = FALSE)
  folds = lapply(1:K, function(k) permus[seq(k, n, by = K)])
  
  for (kk in 1:K) {
    foldIdx = folds[[kk]]
    Xtrain = X[-foldIdx, ]
    old_seed <- .Random.seed
    diana_result <- diana(Xtrain)
    cluster_centers <- matrix(0, nrow = 2, ncol = ncol(Xtrain))
    clu_train <- cutree(diana_result, k = 2)
    for(i in 1:2) {
      if(sum(clu_train == i) > 0) {
        cluster_centers[i, ] <- colMeans(Xtrain[clu_train == i, , drop = FALSE])
      }
    }
    .Random.seed <<- old_seed
    
    jj = p1 + 1
    idxj = sample(foldIdx, ceiling(length(foldIdx)/2), replace = FALSE)
    Xtest0 = X[idxj, ]; Xj0 = Xtest0[, jj]
    Xtest1 = X[setdiff(foldIdx, idxj), ]; Xj1 = Xtest1[, jj]
    
    dist_to_center1 <- sqrt(rowSums((Xtest0 - matrix(cluster_centers[1,], nrow=nrow(Xtest0), ncol=ncol(Xtest0), byrow=TRUE))^2))
    dist_to_center2 <- sqrt(rowSums((Xtest0 - matrix(cluster_centers[2,], nrow=nrow(Xtest0), ncol=ncol(Xtest0), byrow=TRUE))^2))
    clu_test0 <- ifelse(dist_to_center1 < dist_to_center2, 1, 2)
    
    dist_to_center1_syn <- sqrt(rowSums((Xtest1 - matrix(cluster_centers[1,], nrow=nrow(Xtest1), ncol=ncol(Xtest1), byrow=TRUE))^2))
    dist_to_center2_syn <- sqrt(rowSums((Xtest1 - matrix(cluster_centers[2,], nrow=nrow(Xtest1), ncol=ncol(Xtest1), byrow=TRUE))^2))
    clu_test1 <- ifelse(dist_to_center1_syn < dist_to_center2_syn, 1, 2)
    
    Xj1_diff = (Xj1 - cluster_centers[3-clu_test1, jj])^2
    Xj0_diff = (Xj0 - cluster_centers[clu_test0, jj])^2
    
    Xj1_squared_diff = c(Xj1_squared_diff, Xj1_diff)
    Xj0_squared_diff = c(Xj0_squared_diff, Xj0_diff)
  }
  
  numerator = sqrt(n/4)*(mean(Xj1_squared_diff) - mean(Xj0_squared_diff))
  denominator = sqrt((var(Xj1_squared_diff) + var(Xj0_squared_diff)) / 2)
  return(numerator / denominator)
}

# 参数设置
p1 = 20; p = 100; n = 400; mu_1 = 0; mu_2 = 0.5; rho = 0

# 两种方差配置
sigma_1_config1 = 3; sigma_2_config1 = 1  # 配置1：方差(3,1)
sigma_1_config2 = 1; sigma_2_config2 = 1  # 配置2：方差(1,1)

n_repetitions <- 50

# 存储所有结果 - 使用命名列表
all_results <- list(
  norm_hetero = list(), norm_homo = list(),
  t_hetero = list(), t_homo = list()
)

# 生成数据并计算结果 - 配置1：异方差(hetero, σ(3,1))
cat("计算异方差配置...\n")
sigma_1 = sigma_1_config1; sigma_2 = sigma_2_config1

# 固定种子随机选择信号维度，确保p1+1维是噪声
old_seed <- .Random.seed
set.seed(12345)  # 固定种子，在整个实验中使用相同的信号位置
available_dims <- setdiff(1:p, p1+1)
signal_dims <- sample(available_dims, size = p1, replace = FALSE)
noise_dims <- setdiff(1:p, signal_dims)
.Random.seed <<- old_seed

# 构建协方差矩阵，根据信号和噪声维度分配方差
Sigma = rho^abs(matrix(rep(1:p,each=p),p,p,byrow = F)-matrix(rep(1:p,each=p),p,p,byrow = T))
variance_vector <- rep(sigma_2, p)  # 默认所有特征使用sigma_2
variance_vector[signal_dims] <- sigma_1  # 信号特征使用sigma_1
Sigma <- sqrt(diag(variance_vector)) %*% Sigma %*% sqrt(diag(variance_vector))

npor = 0.2
n1 = ceiling(n*npor)
n0 = n - n1
center0 = rep(mu_1, p)
center1 = rep(mu_1, p)
center1[signal_dims] = mu_2

# 正态分布异方差
cat("计算正态分布异方差...\n")
result_norm_PCI_km_hetero <- numeric(n_repetitions)
result_norm_CTT_km_hetero <- numeric(n_repetitions)
result_norm_SI_km_hetero <- numeric(n_repetitions)
result_norm_DT_km_hetero <- numeric(n_repetitions)
result_norm_PCI_h_hetero <- numeric(n_repetitions)
result_norm_CTT_h_hetero <- numeric(n_repetitions)
result_norm_DT_h_hetero <- numeric(n_repetitions)

for(l in 1:n_repetitions) {
  set.seed(l+31*1000)
  X1 = t(t(2*matrix(rnorm(n1*p), n1, p) %*% chol(Sigma))+ center1)
  X0 = t(t(2*matrix(rnorm(n0*p), n0, p) %*% chol(Sigma))+ center0) 
  X = rbind(X1, X0)
  
  result_norm_PCI_km_hetero[l] <- PCI_km(X, p1, K=4)
  result_norm_CTT_km_hetero[l] <- CTT_km(X, p1)
  p_value_si <- SI_km(X, p1)
  result_norm_SI_km_hetero[l] <- qnorm(p_value_si)
  result_norm_DT_km_hetero[l] <- DT_km(X, p1)
  result_norm_PCI_h_hetero[l] <- PCI_h(X, p1, K=4)
  result_norm_CTT_h_hetero[l] <- CTT_h(X, p1)
  result_norm_DT_h_hetero[l] <- DT_h(X, p1)
}

# t分布异方差
cat("计算t分布异方差...\n")
df_t = 8
result_t_PCI_km_hetero <- numeric(n_repetitions)
result_t_CTT_km_hetero <- numeric(n_repetitions)
result_t_SI_km_hetero <- numeric(n_repetitions)
result_t_DT_km_hetero <- numeric(n_repetitions)
result_t_PCI_h_hetero <- numeric(n_repetitions)
result_t_CTT_h_hetero <- numeric(n_repetitions)
result_t_DT_h_hetero <- numeric(n_repetitions)

for(l in 1:n_repetitions) {
  set.seed(l+31*2000)
  Sigma_sqrt <- chol(Sigma)
  X1 = matrix(rt(n1*p, df = df_t), n1, p)
  X0 = matrix(rt(n0*p, df = df_t), n0, p)
  X1 = X1 %*% t(Sigma_sqrt)
  X0 = X0 %*% t(Sigma_sqrt)
  X1 = t(t(X1) + center1)
  X0 = t(t(X0) + center0) 
  X = rbind(X1, X0)
  
  result_t_PCI_km_hetero[l] <- PCI_km(X, p1, K=4)
  result_t_CTT_km_hetero[l] <- CTT_km(X, p1)
  p_value_si <- SI_km(X, p1)
  result_t_SI_km_hetero[l] <- qnorm(p_value_si)
  result_t_DT_km_hetero[l] <- DT_km(X, p1)
  result_t_PCI_h_hetero[l] <- PCI_h(X, p1, K=4)
  result_t_CTT_h_hetero[l] <- CTT_h(X, p1)
  result_t_DT_h_hetero[l] <- DT_h(X, p1)
}

# 配置2：同方差(homo, σ(1,1))
cat("计算同方差配置...\n")
sigma_1 = sigma_1_config2; sigma_2 = sigma_2_config2

# 使用同样的信号维度，但构建新的协方差矩阵（方差相同）
# signal_dims已在配置1中定义
Sigma = rho^abs(matrix(rep(1:p,each=p),p,p,byrow = F)-matrix(rep(1:p,each=p),p,p,byrow = T))
variance_vector <- rep(sigma_2, p)  # 同方差：所有特征使用sigma_2
variance_vector[signal_dims] <- sigma_1  # 信号特征使用sigma_1，但这里sigma_1 = sigma_2 = 1，所以是同方差
Sigma <- sqrt(diag(variance_vector)) %*% Sigma %*% sqrt(diag(variance_vector))

# 正态分布同方差
cat("计算正态分布同方差...\n")
result_norm_PCI_km_homo <- numeric(n_repetitions)
result_norm_CTT_km_homo <- numeric(n_repetitions)
result_norm_SI_km_homo <- numeric(n_repetitions)
result_norm_DT_km_homo <- numeric(n_repetitions)
result_norm_PCI_h_homo <- numeric(n_repetitions)
result_norm_CTT_h_homo <- numeric(n_repetitions)
result_norm_DT_h_homo <- numeric(n_repetitions)

for(l in 1:n_repetitions) {
  set.seed(l+11*1000)
  X1 = t(t(2*matrix(rnorm(n1*p), n1, p) %*% chol(Sigma))+ center1)
  X0 = t(t(2*matrix(rnorm(n0*p), n0, p) %*% chol(Sigma))+ center0) 
  X = rbind(X1, X0)
  
  result_norm_PCI_km_homo[l] <- PCI_km(X, p1, K=4)
  result_norm_CTT_km_homo[l] <- CTT_km(X, p1)
  p_value_si <- SI_km(X, p1)
  result_norm_SI_km_homo[l] <- qnorm(p_value_si)
  result_norm_DT_km_homo[l] <- DT_km(X, p1)
  result_norm_PCI_h_homo[l] <- PCI_h(X, p1, K=4)
  result_norm_CTT_h_homo[l] <- CTT_h(X, p1)
  result_norm_DT_h_homo[l] <- DT_h(X, p1)
}

# t分布同方差
cat("计算t分布同方差...\n")
result_t_PCI_km_homo <- numeric(n_repetitions)
result_t_CTT_km_homo <- numeric(n_repetitions)
result_t_SI_km_homo <- numeric(n_repetitions)
result_t_DT_km_homo <- numeric(n_repetitions)
result_t_PCI_h_homo <- numeric(n_repetitions)
result_t_CTT_h_homo <- numeric(n_repetitions)
result_t_DT_h_homo <- numeric(n_repetitions)

for(l in 1:n_repetitions) {
  set.seed(l+11*2000)
  Sigma_sqrt <- chol(Sigma)
  X1 = matrix(rt(n1*p, df = df_t), n1, p)
  X0 = matrix(rt(n0*p, df = df_t), n0, p)
  X1 = X1 %*% t(Sigma_sqrt)
  X0 = X0 %*% t(Sigma_sqrt)
  X1 = t(t(X1) + center1)
  X0 = t(t(X0) + center0) 
  X = rbind(X1, X0)
  
  result_t_PCI_km_homo[l] <- PCI_km(X, p1, K=4)
  result_t_CTT_km_homo[l] <- CTT_km(X, p1)
  p_value_si <- SI_km(X, p1)
  result_t_SI_km_homo[l] <- qnorm(p_value_si)
  result_t_DT_km_homo[l] <- DT_km(X, p1)
  result_t_PCI_h_homo[l] <- PCI_h(X, p1, K=4)
  result_t_CTT_h_homo[l] <- CTT_h(X, p1)
  result_t_DT_h_homo[l] <- DT_h(X, p1)
}

# 存储结果
all_results <- list(
  norm_hetero = list(  # 异方差 (heteroscedastic)
    PCI_km = result_norm_PCI_km_hetero,
    CTT_km = result_norm_CTT_km_hetero,
    SI_km = result_norm_SI_km_hetero,
    DT_km = result_norm_DT_km_hetero,
    PCI_h = result_norm_PCI_h_hetero,
    CTT_h = result_norm_CTT_h_hetero,
    DT_h = result_norm_DT_h_hetero
  ),
  norm_homo = list(  # 同方差 (homoscedastic)
    PCI_km = result_norm_PCI_km_homo,
    CTT_km = result_norm_CTT_km_homo,
    SI_km = result_norm_SI_km_homo,
    DT_km = result_norm_DT_km_homo,
    PCI_h = result_norm_PCI_h_homo,
    CTT_h = result_norm_CTT_h_homo,
    DT_h = result_norm_DT_h_homo
  ),
  t_hetero = list(  # 异方差 (heteroscedastic)
    PCI_km = result_t_PCI_km_hetero,
    CTT_km = result_t_CTT_km_hetero,
    SI_km = result_t_SI_km_hetero,
    DT_km = result_t_DT_km_hetero,
    PCI_h = result_t_PCI_h_hetero,
    CTT_h = result_t_CTT_h_hetero,
    DT_h = result_t_DT_h_hetero
  ),
  t_homo = list(  # 同方差 (homoscedastic)
    PCI_km = result_t_PCI_km_homo,
    CTT_km = result_t_CTT_km_homo,
    SI_km = result_t_SI_km_homo,
    DT_km = result_t_DT_km_homo,
    PCI_h = result_t_PCI_h_homo,
    CTT_h = result_t_CTT_h_homo,
    DT_h = result_t_DT_h_homo
  )
)

# 生成所有图的函数
create_histogram <- function(values, title, mean_val, sd_val) {
  df <- data.frame(values = values)
  ggplot(df, aes(x = values)) +
    geom_histogram(aes(y = after_stat(density)), bins = 20, 
                   fill = "lightblue", color = "black", alpha = 0.7) +
    geom_density(color = "red", size = 1) +
    stat_function(fun = dnorm, 
                  args = list(mean = mean_val, sd = sd_val),
                  color = "blue", linetype = "dashed", size = 1) +
    labs(title = title, x = "Value", y = "Density") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 8))
}

# 创建4行8列的图
# 行1: 正态分布, 异方差(hetero, σ(3,1)), 左边4列K-means (PCI, CTT, SI, DT), 右边3列层次聚类 (PCI, CTT, DT), 中间空1列
# 行2: 正态分布, 同方差(homo, σ(1,1)), 左边4列K-means, 右边3列层次聚类, 中间空1列
# 行3: t分布, 异方差(hetero, σ(3,1)), 左边4列K-means, 右边3列层次聚类, 中间空1列
# 行4: t分布, 同方差(homo, σ(1,1)), 左边4列K-means, 右边3列层次聚类, 中间空1列

plot_list <- list()

# 第1行：正态分布异方差(hetero, σ(3,1))
plot_list[[1]] <- create_histogram(all_results$norm_hetero$PCI_km, "Normal Hetero PCI-Km", 
                                    mean(all_results$norm_hetero$PCI_km), sd(all_results$norm_hetero$PCI_km))
plot_list[[2]] <- create_histogram(all_results$norm_hetero$CTT_km, "Normal Hetero CTT-Km", 
                                    mean(all_results$norm_hetero$CTT_km), sd(all_results$norm_hetero$CTT_km))
plot_list[[3]] <- create_histogram(all_results$norm_hetero$SI_km, "Normal Hetero SI-Km", 
                                    mean(all_results$norm_hetero$SI_km), sd(all_results$norm_hetero$SI_km))
plot_list[[4]] <- create_histogram(all_results$norm_hetero$DT_km, "Normal Hetero DT-Km", 
                                    mean(all_results$norm_hetero$DT_km), sd(all_results$norm_hetero$DT_km))
plot_list[[5]] <- ggplot() + theme_void()
plot_list[[6]] <- create_histogram(all_results$norm_hetero$PCI_h, "Normal Hetero PCI-H", 
                                    mean(all_results$norm_hetero$PCI_h), sd(all_results$norm_hetero$PCI_h))
plot_list[[7]] <- create_histogram(all_results$norm_hetero$CTT_h, "Normal Hetero CTT-H", 
                                    mean(all_results$norm_hetero$CTT_h), sd(all_results$norm_hetero$CTT_h))
plot_list[[8]] <- create_histogram(all_results$norm_hetero$DT_h, "Normal Hetero DT-H", 
                                    mean(all_results$norm_hetero$DT_h), sd(all_results$norm_hetero$DT_h))

# 第2行：正态分布同方差(homo, σ(1,1))
plot_list[[9]] <- create_histogram(all_results$norm_homo$PCI_km, "Normal Homo PCI-Km", 
                                    mean(all_results$norm_homo$PCI_km), sd(all_results$norm_homo$PCI_km))
plot_list[[10]] <- create_histogram(all_results$norm_homo$CTT_km, "Normal Homo CTT-Km", 
                                    mean(all_results$norm_homo$CTT_km), sd(all_results$norm_homo$CTT_km))
plot_list[[11]] <- create_histogram(all_results$norm_homo$SI_km, "Normal Homo SI-Km", 
                                    mean(all_results$norm_homo$SI_km), sd(all_results$norm_homo$SI_km))
plot_list[[12]] <- create_histogram(all_results$norm_homo$DT_km, "Normal Homo DT-Km", 
                                    mean(all_results$norm_homo$DT_km), sd(all_results$norm_homo$DT_km))
plot_list[[13]] <- ggplot() + theme_void()
plot_list[[14]] <- create_histogram(all_results$norm_homo$PCI_h, "Normal Homo PCI-H", 
                                    mean(all_results$norm_homo$PCI_h), sd(all_results$norm_homo$PCI_h))
plot_list[[15]] <- create_histogram(all_results$norm_homo$CTT_h, "Normal Homo CTT-H", 
                                    mean(all_results$norm_homo$CTT_h), sd(all_results$norm_homo$CTT_h))
plot_list[[16]] <- create_histogram(all_results$norm_homo$DT_h, "Normal Homo DT-H", 
                                    mean(all_results$norm_homo$DT_h), sd(all_results$norm_homo$DT_h))

# 第3行：t分布异方差(hetero, σ(3,1))
plot_list[[17]] <- create_histogram(all_results$t_hetero$PCI_km, "t Hetero PCI-Km", 
                                     mean(all_results$t_hetero$PCI_km), sd(all_results$t_hetero$PCI_km))
plot_list[[18]] <- create_histogram(all_results$t_hetero$CTT_km, "t Hetero CTT-Km", 
                                     mean(all_results$t_hetero$CTT_km), sd(all_results$t_hetero$CTT_km))
plot_list[[19]] <- create_histogram(all_results$t_hetero$SI_km, "t Hetero SI-Km", 
                                     mean(all_results$t_hetero$SI_km), sd(all_results$t_hetero$SI_km))
plot_list[[20]] <- create_histogram(all_results$t_hetero$DT_km, "t Hetero DT-Km", 
                                     mean(all_results$t_hetero$DT_km), sd(all_results$t_hetero$DT_km))
plot_list[[21]] <- ggplot() + theme_void()
plot_list[[22]] <- create_histogram(all_results$t_hetero$PCI_h, "t Hetero PCI-H", 
                                     mean(all_results$t_hetero$PCI_h), sd(all_results$t_hetero$PCI_h))
plot_list[[23]] <- create_histogram(all_results$t_hetero$CTT_h, "t Hetero CTT-H", 
                                     mean(all_results$t_hetero$CTT_h), sd(all_results$t_hetero$CTT_h))
plot_list[[24]] <- create_histogram(all_results$t_hetero$DT_h, "t Hetero DT-H", 
                                     mean(all_results$t_hetero$DT_h), sd(all_results$t_hetero$DT_h))

# 第4行：t分布同方差(homo, σ(1,1))
plot_list[[25]] <- create_histogram(all_results$t_homo$PCI_km, "t Homo PCI-Km", 
                                     mean(all_results$t_homo$PCI_km), sd(all_results$t_homo$PCI_km))
plot_list[[26]] <- create_histogram(all_results$t_homo$CTT_km, "t Homo CTT-Km", 
                                     mean(all_results$t_homo$CTT_km), sd(all_results$t_homo$CTT_km))
plot_list[[27]] <- create_histogram(all_results$t_homo$SI_km, "t Homo SI-Km", 
                                     mean(all_results$t_homo$SI_km), sd(all_results$t_homo$SI_km))
plot_list[[28]] <- create_histogram(all_results$t_homo$DT_km, "t Homo DT-Km", 
                                     mean(all_results$t_homo$DT_km), sd(all_results$t_homo$DT_km))
plot_list[[29]] <- ggplot() + theme_void()
plot_list[[30]] <- create_histogram(all_results$t_homo$PCI_h, "t Homo PCI-H", 
                                     mean(all_results$t_homo$PCI_h), sd(all_results$t_homo$PCI_h))
plot_list[[31]] <- create_histogram(all_results$t_homo$CTT_h, "t Homo CTT-H", 
                                     mean(all_results$t_homo$CTT_h), sd(all_results$t_homo$CTT_h))
plot_list[[32]] <- create_histogram(all_results$t_homo$DT_h, "t Homo DT-H", 
                                     mean(all_results$t_homo$DT_h), sd(all_results$t_homo$DT_h))

# 组合所有图
combined_plot <- grid.arrange(
  plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]], plot_list[[5]], plot_list[[6]], plot_list[[7]], plot_list[[8]],
  plot_list[[9]], plot_list[[10]], plot_list[[11]], plot_list[[12]], plot_list[[13]], plot_list[[14]], plot_list[[15]], plot_list[[16]],
  plot_list[[17]], plot_list[[18]], plot_list[[19]], plot_list[[20]], plot_list[[21]], plot_list[[22]], plot_list[[23]], plot_list[[24]],
  plot_list[[25]], plot_list[[26]], plot_list[[27]], plot_list[[28]], plot_list[[29]], plot_list[[30]], plot_list[[31]], plot_list[[32]],
  ncol = 8, nrow = 4
)

ggsave("combined_distribution_comparison.png", combined_plot, width = 20, height = 12, dpi = 300)
cat("图形已保存为 combined_distribution_comparison.png\n")

