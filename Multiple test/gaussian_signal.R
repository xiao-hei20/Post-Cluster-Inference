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

progress <- function(nfin){
  cat(sprintf('%s: tasks completed: %d.\n', Sys.time(), nfin))
}

opts <- list(progress = progress)

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

sta<-function(X, alpha = NULL, df = NULL){
  center_cluster = ClusterR::KMeans_arma(X, clusters = 2)
  clu<-ClusterR::predict_KMeans(X, center_cluster)
  ifel<-function(a,b,c){
    if(a){
      b
    }else{
      c
    }
  }
  numerator<-ifel(length(which(clu==1))==1,X[which(clu==1),],colMeans(X[which(clu==1),]))-ifel(length(which(clu==2))==1,X[which(clu==2),],colMeans(X[which(clu==2),]))
  denominator<-sqrt(ifel(length(which(clu==1))==1,0,apply(X[which(clu==1),], 2, var))/length(which(clu==1))+ifel(length(which(clu==2))==1,0,apply(X[which(clu==2),], 2, var))/length(which(clu==2)))
  t<-numerator/denominator
  
  if(!is.null(alpha) && !is.null(df)){
    p_values <- 2 * (1 - pt(abs(t), df = df))
    return(which(p.adjust(p_values, method = "BH") < alpha))
  }
  return(t)
}

syn_MVNM<-function(X){
  mu<-colMeans(X)
  center_cluster = ClusterR::KMeans_arma(X, clusters = 2)
  clu<-ClusterR::predict_KMeans(X, center_cluster)
  sigma<-(var(X[clu==1,])+var(X[clu==2,]))/2
  syn <- mvrnorm(n,mu ,sigma)
  return(syn)
}

naive_split<-function(X, alpha = NULL, df = NULL){
  n = nrow(X) #样本的数量
  p = ncol(X) #特征的数量
  permus = sample(1:n, n, replace = FALSE) #1到n permute的函数
  folds = lapply(1:2, function(k) permus[seq(k, n, by = 2)] )
  X_train<-X[folds[[1]],]
  X_test<-X[folds[[2]],]
  center_cluster = ClusterR::KMeans_arma(X_train, clusters = 2)
  clu<-ClusterR::predict_KMeans(X_test, center_cluster)
  ifel<-function(a,b,c){
    if(a){
      b
    }else{
      c
    }
  }
  numerator<-ifel(length(which(clu==1))==1,X_test[which(clu==1),],colMeans(X_test[which(clu==1),]))-ifel(length(which(clu==2))==1,X_test[which(clu==2),],colMeans(X_test[which(clu==2),]))
  denominator<-sqrt(ifel(length(which(clu==1))==1,0,apply(X_test[which(clu==1),], 2, var))/length(which(clu==1))+ifel(length(which(clu==2))==1,0,apply(X_test[which(clu==2),], 2, var))/length(which(clu==2)))
  t<-numerator/denominator
  
  if(!is.null(alpha) && !is.null(df)){
    p_values <- 2 * (1 - pt(abs(t), df = df))
    return(which(p.adjust(p_values, method = "BH") < alpha))
  }
  return(t)
}

DS<-function(X,alpha=0.1){
  index<-sample(1:nrow(X),floor(0.5*nrow(X)),replace = F)
  X_1<-X[index,];X_2<-X[-index,]
  sta1<-sta(X_1);sta2<-sta(X_2)
  psi<-sign(sum(sta1*sta2))*sign(sta1*sta2)*(abs(sta1)+abs(sta2))
  return(sym_threshold(psi,alpha,"1")$disc)
}

MDS<-function(X,alpha=0.1,K=5){
  I<-rep(0,ncol(X))
  for (k in 1:K) {
    index<-DS(X,alpha)
    I_k<-rep(0,ncol(X))
    if(length(index)!=0){I_k[index]<-1/length(index)}
    I<-I+I_k
  }
  I<-I/K
  re<-which(I>max(I[sapply(I, function(x){sum(I[I<=x])})<1-alpha]))
  return(re)
}

IVS_cluster <- function(X, K = 2, alpha = NULL, option = NULL, return_both = FALSE){
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
  
  if(!is.null(alpha)){
    if(return_both){
      return(list(
        psi = sym_threshold(psi, alpha, "1")$disc,
        bpsi = sym_threshold(bpsi, alpha, "1")$disc
      ))
    }
    if(option == "psi"){
      return(sym_threshold(psi, alpha, "1")$disc)
    } else if(option == "bpsi"){
      return(sym_threshold(bpsi, alpha, "1")$disc)
    }
  }
  return(list(psi = psi, bpsi = bpsi))
}

IVS_naive <- function(X, K = 2){
  n = nrow(X) #样本的数量
  p = ncol(X) #特征的数量
  Xtrain = X#挑选出一份训练集
  center_cluster = ClusterR::KMeans_arma(Xtrain, clusters = 2, seed = 0)
  psi<-rep(NA,p)
  beta<-rep(NA,p)
  # mean(apply(Xtrain - center_cluster[ClusterR::predict_KMeans(Xtrain, center_cluster), ], 1, function(X){ sqrt(sum(X^2)) }))
    for (jj in 1:p){
      idxj = sample(1:n,floor(n/2), replace = FALSE)#将测试集分成两份
      Xtest0 = X[idxj, ]#挑选出ori
      Xj0 = Xtest0[, jj]#ori的第jj个变量
      Xtest1 = X[setdiff(1:n, idxj), ]#挑选出syn
      Xj1 = Xtest1[, jj]#syn的第jj个变量
      psi[jj] = mean((Xj1 - center_cluster[3-ClusterR::predict_KMeans(Xtest1, center_cluster), jj])^2) -
        mean((Xj0 - center_cluster[ClusterR::predict_KMeans(Xtest0, center_cluster), jj])^2)
      beta[jj] = mean((Xj1 - center_cluster[3-ClusterR::predict_KMeans(Xtest1, center_cluster), jj])^2+
                               (Xj0 - center_cluster[3-ClusterR::predict_KMeans(Xtest0, center_cluster), jj])^2)
    }
  bpsi = psi*beta
  return(list(psi = psi, bpsi = bpsi))
}

IVS_split_cluster<-function(X,K=5){
  n = nrow(X) #样本的数量
  p = ncol(X) #特征的数量
  permus = sample(1:n, n, replace = FALSE) #1到n permute的函数
  folds = lapply(1:K, function(k) permus[seq(k, n, by = K)] ) #样本分割的函数
  for (k_1 in 1:(K-1)) {
    for (k_2 in (k_1+1):K){
      Xtrain = X[-c(folds[[k_1]],folds[[k_2]]), ]#挑选出一份训练集
      center_cluster = ClusterR::KMeans_arma(Xtrain, clusters = 2, seed = 0)
      psi_k<-rep(NA,p)
      for (jj in 1:p) {
        idk1_1 = sample(folds[[k_1]], ceiling(length(folds[[k_1]])/2), replace = FALSE)
        idk1_2 = setdiff(folds[[k_1]],idk1_1)
        Xk1_1=X[idk1_1,]
        Xjk1_1= Xk1_1[, jj]
        Xk1_2=X[idk1_2,]
        Xjk1_2= Xk1_2[, jj]
        psi_k1_jj = mean((Xjk1_1 - center_cluster[3-ClusterR::predict_KMeans(Xk1_1, center_cluster), jj])^2) -
          mean((Xjk1_2 - center_cluster[ClusterR::predict_KMeans(Xk1_2, center_cluster), jj])^2)
        idk2_1 = sample(folds[[k_2]], ceiling(length(folds[[k_2]])/2), replace = FALSE)
        idk2_2 = setdiff(folds[[k_2]],idk2_1)
        Xk2_1=X[idk2_1,]
        Xjk2_1= Xk2_1[, jj]
        Xk2_2=X[idk2_2,]
        Xjk2_2= Xk2_2[, jj]
        psi_k2_jj = mean((Xjk2_1 - center_cluster[3-ClusterR::predict_KMeans(Xk2_1, center_cluster), jj])^2) -
          mean((Xjk2_2 - center_cluster[ClusterR::predict_KMeans(Xk2_2, center_cluster), jj])^2)
        psi_k[jj]<- psi_k1_jj*psi_k2_jj
      }
      if(k_1==1 & k_2==2){
        psi<-psi_k
      }else{
        psi<-cbind(psi,psi_k)
      }
    }
  }
  psi<-rowMeans(psi)
  return(psi)
}

CADET<-function(X, alpha = NULL){
  p<-ncol(X)
  p_value<-rep(NA,p)
  naive_p<-rep(NA,p)
  for (i in 1:p) {
    list<-CADET::kmeans_inference_1f(X,k = 2, 1, 2,feat = i, iso = F,covMat = cov(X), seed = 2021,iter.max = 30)
    p_value[i]<-list$pval
    naive_p[i]<-list$p_naive
  }
  
  if(!is.null(alpha)){
    return(which(p.adjust(p_value, method = "BH") < alpha))
  }
  return(list(p_value=p_value,naive_p=naive_p))
}

norm_split<-function(X,epi=0.5, alpha = NULL, df = NULL){
  # sigma<-diag(rep(1,ncol(X)))
  center_cluster = ClusterR::KMeans_arma(X, clusters = 2)
  clu<-ClusterR::predict_KMeans(X, center_cluster)
  sigma<-(var(X[clu==1,])+var(X[clu==2,]))/2
  # sigma<-var(X)
  X_test<-mvrnorm(n=nrow(X),mu=rep(0,ncol(X)),Sigma = epi*(1-epi)*sigma)+epi*X
  X_train<-X-X_test
  center_cluster = ClusterR::KMeans_arma(X_train, clusters = 2)
  clu<-ClusterR::predict_KMeans(X_train, center_cluster)
  ifel<-function(a,b,c){
    if(a){
      b
    }else{
      c
    }
  }
  numerator<-ifel(length(which(clu==1))==1,X_test[which(clu==1),],colMeans(X_test[which(clu==1),]))-ifel(length(which(clu==2))==1,X_test[which(clu==2),],colMeans(X_test[which(clu==2),]))
  denominator<-sqrt(ifel(length(which(clu==1))==1,0,apply(X_test[which(clu==1),], 2, var))/length(which(clu==1))+ifel(length(which(clu==2))==1,0,apply(X_test[which(clu==2),], 2, var))/length(which(clu==2)))
  t<-numerator/denominator
  
  if(!is.null(alpha) && !is.null(df)){
    p_values <- 2 * (1 - pt(abs(t), df = df))
    return(which(p.adjust(p_values, method = "BH") < alpha))
  }
  return(t)
}

cluster_DE<-function(X, alpha){
  n = nrow(X)
  p = ncol(X)
  syn_X<-syn_MVNM(X)
  C<-log10(2*(1-pt(abs(sta(syn_X)),df=n-2))+runif(p,0,0.00001))-log10(2*(1-pt(abs(sta(X)),df=n-2))+runif(p,0,0.00001))
  return(sym_threshold(C,alpha,"1")$disc)
}
# npor=0.6;num_sim=30;alpha=0.2
# n=2000;p=400;p1=80
# mu1=0.4;mu2=0;Sigma<-diag(rep(1,p))

# npor=0.6;alpha=0.2 
# n=2000;p=400;p1=10 
# prob1=0.1;prob2=0.2;prob3=0.3
# r1=1;r2=3;r3=3 
# n1 = ceiling(n*npor) #1号细胞群的数量 n0 = n - n1
# 
# times=20
# power_IVS<-rep(NA,times);power_en<-rep(NA,times);power_split<-rep(NA,times)
# fdr_IVS<-rep(NA,times);fdr_en<-rep(NA,times);fdr_split<-rep(NA,times)
# for (i in 1:times) {
#   set.seed(20251018 + i*10000)
#   # loc = 1:p1
#   # center0 = rep(mu1, p);
#   # center1 = rep(mu1, p); center1[loc] = mu2
#   # n1 = ceiling(n*npor)
#   # n0 = n - n1
#   # Z = c(rep(1, n1), rep(0, n0))
#   # X1 = t(t(matrix(rnorm(n1*p), n1, p) %*% chol(Sigma)) + center1)
#   # X0 = t(t(matrix(rnorm(n0*p), n0, p) %*% chol(Sigma)) + center0)
#   # X = rbind(X1, X0)
#   X_LU<-matrix(rnbinom(n = n0*p1, size = r1, prob = prob1), nrow = n0, ncol = p1)
#   X_LD<-matrix(rnbinom(n = n1*p1, size = r2, prob = prob2), nrow = n1, ncol = p1)
#   X_R<-matrix(rnbinom(n = n*(p-p1), size = r3, prob = prob3), nrow = n, ncol = p-p1)
#   X<-cbind(rbind(X_LU,X_LD),X_R)
#   lis<-IVS_cluster(X,K=10)
#   psi_split<-IVS_split_cluster(X,K=10)
#   IVS<-sym_threshold(lis$psi,alpha,c("1"))$disc
#   power_IVS[i]<-sum(IVS<=p1)/p1;fdr_IVS[i]<-sum(IVS>p1)/max(length(IVS),1)
#   IVS_en<-sym_threshold(lis$bpsi,alpha,c("1"))$disc
#   power_en[i]<-sum(IVS_en<=p1)/p1;fdr_en[i]<-sum(IVS_en>p1)/max(length(IVS_en),1)
#   IVS_split<-sym_threshold(psi_split,alpha,c("1"))$disc
#   power_split[i]<-sum(IVS_split<=p1)/p1;fdr_split[i]<-sum(IVS_split>p1)/max(length(IVS_split),1)
# }
# table<-matrix(rep(NA,6),2,3) 
# colnames(table)<-c("IVS","IVS_en","IVS_split")
# rownames(table)<-c("FDR","Power")
# table[1,]<-c(mean(fdr_IVS),mean(fdr_en),mean(fdr_split))
# table[2,]<-c(mean(power_IVS),mean(power_en),mean(power_split)) 
# table
# cat(sprintf('%s: tasks begin...\n', Sys.time()))
a=sample(1:100000,1)
b=sample(1:100000,1)


res_list<-list()
seq_tim <- 4
p1_seq     <- seq(from = 75, to = 75, length.out = seq_tim)
p_seq      <- seq(from = 400, to = 400, length.out = seq_tim)
n_seq      <- seq(from = 800, to = 800, length.out = seq_tim)
mu_1_seq   <- seq(from = 1, to = 1, length.out = seq_tim)
mu_2_seq   <- seq(from = 2, to = 3.2, length.out = seq_tim)
sigma_1_seq <- seq(from = 2, to = 2, length.out = seq_tim)
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

# 固定信号位置 loc（在所有参数组合和所有模拟中保持一致）
set.seed(2025)
loc <- sort(sample(1:p_seq[1], p1_seq[1], replace = FALSE))
cat("Fixed signal locations (loc):", loc[seq_len(min(20, length(loc)))], ifelse(length(loc)>20, "...", ""), "\n")

# 自动检测变化参数
var_flags <- sapply(param_list, function(v) length(unique(v)) > 1)
x_param <- names(param_list)[which(var_flags)[1]]
x_values <- param_list[[x_param]]
cat("Varying parameter:", x_param, "with values:", paste(x_values, collapse=", "), "\n")

for (i in 1:seq_tim) {
  param_names <- names(param_list)
  base_names <- sub("_seq$", "", param_names)
  for (j in seq_along(param_names)) {
    assign(base_names[j], param_list[[param_names[j]]][i], envir = environment())
  }
  cl <- makeSOCKcluster(10)
  registerDoSNOW(cl)
  num_sim=50
  res = foreach (num=1:num_sim, .combine = 'cbind', .multicombine = TRUE, .options.snow = opts,.export=c("mvrnorm","sym_threshold","kmeans_inference_1f","IVS_cluster","sta","naive_split","norm_split","syn_MVNM","cluster_DE","CADET","loc")) %dopar% {
    set.seed(a*i+b*num)
    Sigma=rho^abs(matrix(rep(1:p,each=p),p,p,byrow = F)-matrix(rep(1:p,each=p),p,p,byrow = T))
    Sigma=sqrt(diag(c(rep(sigma_1,p1),rep(sigma_2,p-p1))))%*%Sigma%*%sqrt(diag(c(rep(sigma_1,p1),rep(sigma_2,p-p1))))
    npor=0.6;alpha=0.2
    n1 = ceiling(n*npor) #number
    n0 = n - n1
    K=8
    # loc 已在主循环外固定，使用固定的 loc
    center0 = rep(mu_1, p);
    center1 = rep(mu_1, p); center1[loc] = mu_2
    n1 = ceiling(n*npor)
    n0 = n - n1
    Z = c(rep(1, n1), rep(0, n0))
    X1 = t(t(matrix(rnorm(n1*p), n1, p) %*% chol(Sigma))+ center1)
    X0 = t(t(matrix(rnorm(n0*p), n0, p) %*% chol(Sigma))+ center0) 
    X = rbind(X1, X0)
    # X_LU<-matrix(rnbinom(n = n0*p1, size = r1, prob = prob1), nrow = n0, ncol = p1)
    # X_LD<-matrix(rnbinom(n = n1*p1, size = r2, prob = prob2), nrow = n1, ncol = p1)
    # X_R<-matrix(rnbinom(n = n*(p-p1), size = r3, prob = prob3), nrow = n, ncol = p-p1)
    # X<-cbind(rbind(X_LU,X_LD),X_R)
    res_PCI<-IVS_cluster(X, K=K, alpha=alpha, option="bpsi")
    # psi_split<-IVS_split_cluster(X,K=K)
    # res_split<-sym_threshold(psi_split,alpha,"1")$disc
    res_MDS<-MDS(X,alpha,K=K)
    res_DS<-DS(X,alpha)
    res_CAD<-CADET(X, alpha)
    res_nai_1<-sta(X, alpha=alpha, df=n-2)
    #######cluster DE#######
    res_CD<-cluster_DE(X, alpha)
    #######norm_split########
    res_countsplit<-norm_split(X, alpha=alpha, df=n-2)
    return(list(PCI=res_PCI,count_split=res_countsplit,DS=res_DS,MDS=res_MDS,CAD=res_CAD,Cluster_DE=res_CD,naive_1=res_nai_1))
  } 
  stopCluster(cl)
  # save(res, file = file_name)
  cat(sprintf('%s: tasks end and results saved.\n\n', Sys.time()))
  
  # nmthd = 6
  # out = matrix(NA, 2, nmthd)
  # loc=1:p1
  # for (mm in 1:nmthd){
  #   out[, mm] = rowMeans(sapply(1:num_sim, function(num){
  #     symm = res[[mm, num]]
  #     c(length(setdiff(symm, loc))/max(1, length(symm)), 
  #       length(intersect(symm, loc))/max(1, p1))
  #   }))
  # }
  # colnames(out)<-c("IVS","IVS_en","MDS","CADET","Naive")
  res_list[[i]]<-res
  
  save(res_list,file = "res_list.rds")
}

name<-row.names(res_list[[1]])
nmthd=length(name)

# 重写 FDR/Power 计算，使用固定的 loc
cat("Calculating FDR and Power with fixed loc...\n")
stat_list <- vector('list', length=seq_tim)

for (i in 1:seq_tim) {
  res<-res_list[[i]]
  p1<-length(loc)  # 使用固定的 loc 长度
  
  # 计算每个方法的 FDR 和 Power
  fdr_mean <- sapply(1:nmthd, function(m){
    mean(sapply(1:num_sim, function(num){
      symm <- res[[m, num]]
      length(setdiff(symm, loc))/max(1, length(symm))
    }))
  })
  
  power_mean <- sapply(1:nmthd, function(m){
    mean(sapply(1:num_sim, function(num){
      symm <- res[[m, num]]
      length(intersect(symm, loc))/max(1, p1)
    }))
  })
  
  mat <- rbind(FDR=fdr_mean, Power=power_mean)
  colnames(mat) <- name
  stat_list[[i]] <- mat
  
  cat(sprintf("  Parameter set %d/%d done\n", i, seq_tim))
}

cat("FDR and Power calculation completed.\n")

# 生成长表（长格式，包含 Metric 列）
cat("Generating long-format dataframe...\n")
library(tidyr)
library(dplyr)

long_df <- do.call(rbind, lapply(1:seq_tim, function(i){
  df <- as.data.frame(t(stat_list[[i]]))  # 转置：方法为行，FDR/Power 为列
  df$Method <- rownames(df)
  rownames(df) <- NULL
  
  # 使用 pivot_longer 转换为长格式
  df <- df %>% tidyr::pivot_longer(cols=c('FDR','Power'), names_to='Metric', values_to='Value')
  df$ParamName <- x_param
  df$ParamValue <- x_values[i]
  
  # 调整列顺序
  df <- df[, c('ParamName', 'ParamValue', 'Method', 'Metric', 'Value')]
  df
}))

cat("Long-format dataframe created with", nrow(long_df), "rows.\n")
cat("Preview of long_df:\n")
print(head(long_df))

# 保存长表到 CSV
write.csv(long_df, "power_fdr_long.csv", row.names = FALSE)
cat("Long-format results saved to power_fdr_long.csv\n")

# 绘图：左右并列的 FDR/Power 折线图
cat("Creating visualization plots...\n")
library(ggplot2)
library(patchwork)

# 提取 alpha 值（用于 FDR 基准线）
alpha_val <- 0.2

# Power 图
p_power <- ggplot(long_df %>% filter(Metric=='Power'), 
                  aes(x=ParamValue, y=Value, color=Method)) +
  geom_line(size=1) + 
  geom_point(size=2) + 
  labs(x=x_param, y='Power', title='Power Comparison') + 
  theme_bw() +
  theme(legend.position='bottom', plot.title = element_text(hjust = 0.5))

# FDR 图
p_fdr <- ggplot(long_df %>% filter(Metric=='FDR'), 
                aes(x=ParamValue, y=Value, color=Method)) +
  geom_line(size=1) + 
  geom_point(size=2) + 
  geom_hline(yintercept=alpha_val, linetype='dashed', color='red', alpha=0.7) +
  labs(x=x_param, y='FDR', title='FDR Comparison') + 
  theme_bw() +
  theme(legend.position='bottom', plot.title = element_text(hjust = 0.5))

# 组合两个图
combined_plot <- p_power | p_fdr
combined_plot <- combined_plot + plot_annotation(
  title = paste('Method Comparison: Power and FDR vs', x_param),
  theme = theme(plot.title = element_text(hjust = 0.5, size = 14, face = 'bold'))
)

# 保存图形
plot_filename <- paste0('power_fdr_plot_', x_param, '.png')
ggsave(plot_filename, combined_plot, width = 16, height = 8, dpi = 300)
cat("Plot saved to", plot_filename, "\n")

# 显示图形
print(combined_plot)

cat("\n===== Summary =====\n")
cat("stat_list: A list with", length(stat_list), "elements (one per parameter set)\n")
cat("Each element is a 2 x", nmthd, "matrix (rows: FDR, Power; columns: methods)\n")
cat("long_df: A long-format dataframe with", nrow(long_df), "rows\n")
cat("Column structure:", paste(colnames(long_df), collapse=", "), "\n")
cat("===================\n")

stat_list

# # 保持原有代码（用于保存 wide 格式结果）
# param_names <- names(param_list)
# param_df <- as.data.frame(param_list)
# 
# wide_result <- do.call(rbind, lapply(1:seq_tim, function(i) {
#   res_i <- as.data.frame(stat_list[[i]])
#   names(res_i) <- colnames(stat_list[[i]])
#   param_row <- param_df[i, , drop = FALSE]
#   sim_row <- cbind(a=a,b=b,times=num_sim)
#   # 把两行分别拉成一列（metric 为列名）
#   df<-as.data.frame(t(res_i))
#   df$Method <-colnames(res_i)
#   rownames(df)<-NULL
#   colnames(df)<-c("FDR","Power","Method")
#   
#   # 加上参数列
#   df <- cbind(param_row[rep(1, nrow(df)),], df)
#   df <- cbind(df,sim_row[rep(1, nrow(df)),])
#   return(df)
# }))
# 
# if(file.exists("gausigtestresultwide.csv")){
#   res_pr<-read.csv("gausigtestresultwide.csv")[,-1]
#   wide_result<-rbind(wide_result,res_pr)
#   write.csv(wide_result,"gausigtestresultwide.csv")
# }else{
#   write.csv(wide_result,"gausigtestresultwide.csv")
# }
