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
p1 = 80; p = 400; n = 800; mu_1 = 0; mu_2 = 3; rho = 0.2

# 两种方差配置
sigma_1_config1 = 2; sigma_2_config1 = 1  # 配置1：方差(3,1)
sigma_1_config2 = 1; sigma_2_config2 = 1  # 配置2：方差(1,1)

n_repetitions <- 100

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

npor = 0.5
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
  set.seed(l+35*1000)
  X1 = t(t(matrix(rnorm(n1*p), n1, p) %*% chol(Sigma))+ center1)
  X0 = t(t(matrix(rnorm(n0*p), n0, p) %*% chol(Sigma))+ center0) 
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
  set.seed(l+35*2000)
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
  X1 = t(t(matrix(rnorm(n1*p), n1, p) %*% chol(Sigma))+ center1)
  X0 = t(t(matrix(rnorm(n0*p), n0, p) %*% chol(Sigma))+ center0) 
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

# 创建一个函数来生成图
create_plot <- function(data_homo, data_hetero, top_title, right_title, filename, has_si = TRUE, show_top_title = TRUE, show_right_title = TRUE, 
                         show_facet_row_title = TRUE, show_facet_col_title_top = TRUE, show_facet_col_title_bottom = FALSE) {
  # 确定方法列表
  if (has_si) {
    methods <- c("PCI", "CTT", "SI", "DT")
    method_values_homo <- c(data_homo$PCI_km, data_homo$CTT_km, data_homo$SI_km, data_homo$DT_km)
    method_values_hetero <- c(data_hetero$PCI_km, data_hetero$CTT_km, data_hetero$SI_km, data_hetero$DT_km)
  } else {
    methods <- c("PCI", "CTT", "DT")
    method_values_homo <- c(data_homo$PCI_km, data_homo$CTT_km, data_homo$DT_km)
    method_values_hetero <- c(data_hetero$PCI_km, data_hetero$CTT_km, data_hetero$DT_km)
  }
  
  # 创建Homo数据框
  homo_data <- data.frame(
    value = method_values_homo,
    method = rep(methods, each = length(data_homo$PCI_km)),
    condition = "Homo"
  )
  
  # 创建Hetero数据框
  hetero_data <- data.frame(
    value = method_values_hetero,
    method = rep(methods, each = length(data_hetero$PCI_km)),
    condition = "Hetero"
  )
  
  # 生成正态曲线数据
  generate_curves <- function(data) {
    curves <- data.frame()
    for (i in seq_along(unique(data$method))) {
      method_name <- unique(data$method)[i]
      data_subset <- data[data$method == method_name, ]
      
    mean_val <- 0
    sd_val <- 1
    
    x_seq <- seq(min(data_subset$value), max(data_subset$value), length.out = 100)
    y_seq <- dnorm(x_seq, mean = mean_val, sd = sd_val)
    
      curves <- rbind(curves, data.frame(
      value = x_seq,
      density = y_seq,
      method = method_name,
        condition = data$condition[1]
      ))
    }
    return(curves)
  }
  
  homo_normal_curves <- generate_curves(homo_data)
  hetero_normal_curves <- generate_curves(hetero_data)
  
  # 创建Homo图 - 根据参数控制标题显示
  homo_theme <- theme_minimal() +
  theme(
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    axis.ticks = element_line(size = 0.5),
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_rect(fill = "gray80", color = "black", linewidth = 0.8),
    strip.placement = "outside",
    panel.spacing = unit(0.5, "cm"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.0),
    panel.grid.major = element_line(color = "gray70", linewidth = 0.5),
      panel.grid.minor = element_line(color = "gray85", linewidth = 0.3)
    )
  
  # 根据参数设置标题
  if (!show_facet_col_title_top) {
    homo_theme <- homo_theme + theme(strip.text.x = element_blank())
  }
  if (!show_facet_row_title) {
    homo_theme <- homo_theme + theme(strip.text.y = element_blank())
  }
  
  homo_plot <- ggplot(homo_data, aes(x = value)) +
  geom_histogram(aes(y = after_stat(density)), bins = 20, 
                 fill = "lightblue", color = "black", alpha = 0.7) +
  geom_density(color = "red", size = 1) +
    geom_line(data = homo_normal_curves, aes(y = density), color = "blue", 
            linetype = "dashed", size = 1) +
  labs(x = "", y = "") +
    homo_theme +
    facet_grid(condition ~ method, scales = "free_x")
  
  # 创建Hetero图 - 根据参数控制标题显示
  hetero_theme <- theme_minimal() +
  theme(
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    axis.ticks = element_line(size = 0.5),
      strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_rect(fill = "gray80", color = "black", linewidth = 0.8),
    strip.placement = "outside",
    panel.spacing = unit(0.5, "cm"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.0),
    panel.grid.major = element_line(color = "gray70", linewidth = 0.5),
    panel.grid.minor = element_line(color = "gray85", linewidth = 0.3)
    )
  
  # 根据参数设置标题
  if (!show_facet_col_title_bottom) {
    hetero_theme <- hetero_theme + theme(strip.text.x = element_blank())
  }
  if (!show_facet_row_title) {
    hetero_theme <- hetero_theme + theme(strip.text.y = element_blank())
  }
  
  hetero_plot <- ggplot(hetero_data, aes(x = value)) +
    geom_histogram(aes(y = after_stat(density)), bins = 20, 
                   fill = "lightblue", color = "black", alpha = 0.7) +
    geom_density(color = "red", size = 1) +
    geom_line(data = hetero_normal_curves, aes(y = density), color = "blue", 
              linetype = "dashed", size = 1) +
    labs(x = "", y = "") +
    hetero_theme +
    facet_grid(condition ~ method, scales = "free_x")
  
  # 拼接两个图
combined_plot <- grid.arrange(
    homo_plot,
    hetero_plot,
  nrow = 2,
  heights = c(1, 1)
)

  # 根据参数决定是否添加标题
  if (show_top_title && show_right_title) {
    # 两个标题都要
    title_top <- grid::textGrob(top_title, gp = grid::gpar(fontsize = 16, fontface = "bold"))
    title_right <- grid::textGrob(right_title, rot = 270, gp = grid::gpar(fontsize = 16, fontface = "bold"))
final_plot <- grid.arrange(
  title_top,
  grid.arrange(combined_plot, title_right, ncol = 2, widths = c(20, 1)),
  nrow = 2,
  heights = c(0.5, 10)
    )
  } else if (show_top_title && !show_right_title) {
    # 只要顶部标题
    title_top <- grid::textGrob(top_title, gp = grid::gpar(fontsize = 16, fontface = "bold"))
    final_plot <- grid.arrange(
      title_top,
      combined_plot,
      nrow = 2,
      heights = c(0.5, 10)
    )
  } else if (!show_top_title && show_right_title) {
    # 只要右侧标题
    title_right <- grid::textGrob(right_title, rot = 270, gp = grid::gpar(fontsize = 16, fontface = "bold"))
    final_plot <- grid.arrange(
      combined_plot, title_right, ncol = 2, widths = c(20, 1)
    )
  } else {
    # 都不要
    final_plot <- combined_plot
  }
  
  # 返回图形对象
  return(final_plot)
}

# 生成4张图并保存到变量
# 1. K-means + Normal - 只要顶部标题，不要右侧标题；不显示facet行标题
p1 <- create_plot(all_results$norm_homo, all_results$norm_hetero, "K-means", "Normal", "km_normal.png", 
            show_top_title = TRUE, show_right_title = FALSE, 
            show_facet_row_title = FALSE, show_facet_col_title_top = TRUE, show_facet_col_title_bottom = FALSE)

# 2. K-means + t分布 - 两个标题都不要；facet行标题和列标题都不显示
p2 <- create_plot(all_results$t_homo, all_results$t_hetero, "K-means", "t-distribution", "km_t.png", 
            show_top_title = FALSE, show_right_title = FALSE,
            show_facet_row_title = FALSE, show_facet_col_title_top = FALSE, show_facet_col_title_bottom = FALSE)

# 3. Hierarchical + Normal (只使用PCI, CTT, DT，没有SI) - 两个标题都要；facet行列标题都显示
hierarchy_norm_homo <- list(
  PCI_km = all_results$norm_homo$PCI_h, 
  CTT_km = all_results$norm_homo$CTT_h, 
  DT_km = all_results$norm_homo$DT_h
)
hierarchy_norm_hetero <- list(
  PCI_km = all_results$norm_hetero$PCI_h, 
  CTT_km = all_results$norm_hetero$CTT_h, 
  DT_km = all_results$norm_hetero$DT_h
)
p3 <- create_plot(hierarchy_norm_homo, hierarchy_norm_hetero, "Hierarchical", "Normal", "h_normal.png", 
            has_si = FALSE, show_top_title = TRUE, show_right_title = TRUE,
            show_facet_row_title = TRUE, show_facet_col_title_top = TRUE, show_facet_col_title_bottom = FALSE)

# 4. Hierarchical + t分布 - 只要右侧标题；facet列标题不显示
hierarchy_t_homo <- list(
  PCI_km = all_results$t_homo$PCI_h, 
  CTT_km = all_results$t_homo$CTT_h, 
  DT_km = all_results$t_homo$DT_h
)
hierarchy_t_hetero <- list(
  PCI_km = all_results$t_hetero$PCI_h, 
  CTT_km = all_results$t_hetero$CTT_h, 
  DT_km = all_results$t_hetero$DT_h
)
p4 <- create_plot(hierarchy_t_homo, hierarchy_t_hetero, "Hierarchical", "t-distribution", "h_t.png", 
            has_si = FALSE, show_top_title = FALSE, show_right_title = TRUE,
            show_facet_row_title = TRUE, show_facet_col_title_top = FALSE, show_facet_col_title_bottom = FALSE)

# 创建空白间隔
blank_v <- grid::rectGrob(gp = grid::gpar(col = "white", fill = "white"))  # 垂直间隔
blank_h <- grid::rectGrob(gp = grid::gpar(col = "white", fill = "white"))  # 水平间隔

# 拼接四张图，添加间隔
final_arrangement <- grid.arrange(
  p1, blank_v, p3,  # 第一行：Normal (K-means, 空白, Hierarchical)
  blank_h, blank_h, blank_h,  # 行间隔
  p2, blank_v, p4,  # 第二行：t-distribution (K-means, 空白, Hierarchical)
  ncol = 3,
  nrow = 3,
  heights = c(1, 0, 1),  # 中间行高度较小作为间隔
  widths = c(4, 0.3, 3)   # 中间列宽度较小作为间隔
)

# 保存最终图
ggsave("combined_4_panels.png", final_arrangement, width = 20, height = 14, dpi = 300)
cat("最终拼接图已保存为 combined_4_panels.png\n")

