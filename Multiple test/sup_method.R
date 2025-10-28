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

sta<-function(X){
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

naive_split<-function(X){
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

CADET<-function(X){
  p<-ncol(X)
  p_value<-rep(NA,p)
  naive_p<-rep(NA,p)
  for (i in 1:p) {
    list<-CADET::kmeans_inference_1f(X,k = 2, 1, 2,feat = i, iso = F,covMat = cov(X), seed = 2021,iter.max = 30)
    p_value[i]<-list$pval
    naive_p[i]<-list$p_naive
  }
  return(list(p_value=p_value,naive_p=naive_p))
}

norm_split<-function(X,epi=0.5){
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
  return(t)
}

a=sample(1:100000,1)
b=sample(1:100000,1)


res_list<-list()
seq_tim <- 4
p1_seq     <- seq(from = 75, to = 75, length.out = seq_tim)
p_seq      <- seq(from = 400, to = 400, length.out = seq_tim)
n_seq      <- seq(from = 200, to = 200, length.out = seq_tim)
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

for (i in 1:seq_tim) {
  param_names <- names(param_list)
  base_names <- sub("_seq$", "", param_names)
  for (j in seq_along(param_names)) {
    assign(base_names[j], param_list[[param_names[j]]][i], envir = environment())
  }
  cl <- makeSOCKcluster(10)
  registerDoSNOW(cl)
  num_sim=100
  res = foreach (num=1:num_sim, .combine = 'cbind', .multicombine = TRUE, .options.snow = opts,.export=c("mvrnorm","sym_threshold","kmeans_inference_1f")) %dopar% {
    set.seed(a*i+b*num)
    Sigma=rho^abs(matrix(rep(1:p,each=p),p,p,byrow = F)-matrix(rep(1:p,each=p),p,p,byrow = T))
    Sigma=sqrt(diag(c(rep(sigma_1,p1),rep(sigma_2,p-p1))))%*%Sigma%*%sqrt(diag(c(rep(sigma_1,p1),rep(sigma_2,p-p1))))
    npor=0.6;alpha=0.2
    n1 = ceiling(n*npor) #number
    n0 = n - n1
    K=8
    loc = 1:p1
    center0 = rep(mu_1, p);
    center1 = rep(mu_1, p); center1[loc] = mu_2
    n1 = ceiling(n*npor)
    n0 = n - n1
    Z = c(rep(1, n1), rep(0, n0))
    X1 = t(t(2*matrix(rnorm(n1*p), n1, p) %*% chol(Sigma))+ center1)
    X0 = t(t(2*matrix(rnorm(n0*p), n0, p) %*% chol(Sigma))+ center0) 
    X = rbind(X1, X0)
    p_CAD<-CADET(X)
    res_CAD<-which(p.adjust(p_CAD$p_value,method = "BH")<alpha)
    #######norm_split########
    p_countsplit<-2*(1-pt(abs(norm_split(X,epi = 0.5)),df=n-2))
    res_countsplit<-which(p.adjust(p_countsplit,method = "BH")<alpha)
    return(list(count_split=res_countsplit,CAD=res_CAD))
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
  
  save(res_list,file = "res_sup_list.rds")
}

name<-row.names(res_list[[1]])

FDR<-list()
Power<-list()
nmthd=length(name)
for (i in 1:seq_tim) {
  fdr_i<-matrix(NA,nrow = num_sim,ncol = nmthd)
  power_i<-matrix(NA,nrow = num_sim,ncol = nmthd)
  for (method in 1:nmthd) {
    p1<-p1_seq[i]
    res<-res_list[[i]]
    loc=1:p1
    res_method<-sapply(1:num_sim, function(num){
      symm = res[[method, num]]
      c(length(setdiff(symm, loc))/max(1, length(symm)),
        length(intersect(symm, loc))/max(1, p1))
    })
    fdr_i[,method]<-res_method[1,]
    power_i[,method]<-res_method[2,]
  }
  
  FDR[[i]]<-fdr_i
  Power[[i]]<-power_i
}
result<-list()
for (i in 1:seq_tim) {
  resulti<-rbind(colMeans(FDR[[i]]),colMeans(Power[[i]]))
  # colnames(resulti)<-c("IVS","IVS_en","count_split","DS","MDS","CADET","Cluster_DE","Naive_1","Naive_2")
  colnames(resulti)<-name
  result[[i]]<-resulti
}
result


FDR<-matrix(NA,nrow=seq_tim,ncol=nmthd)
# colnames(FDR)<-c("IVS","IVS_en","count_split","DS","MDS","CADET","Cluster_DE","Naive_1","Naive_2")
colnames(FDR)<-name
Power<-matrix(NA,nrow=seq_tim,ncol=nmthd)
# colnames(Power)<-c("IVS","IVS_en","count_split","DS","MDS","CADET","Cluster_DE","Naive_1","Naive_2")
colnames(Power)<-name
for (i in 1:seq_tim) {
  result[[i]][1,]->FDR[i,]
  result[[i]][2,]->Power[i,]
}


# 参数准备
param_names <- names(param_list)
param_df <- as.data.frame(param_list)

# 绑定结果和参数（每组参数 -> 一行 FDR、一行 Power → 我们 reshape 为 wide）
library(dplyr)

wide_result <- do.call(rbind, lapply(1:seq_tim, function(i) {
  res_i <- as.data.frame(result[[i]])
  names(res_i) <- colnames(result[[i]])
  param_row <- param_df[i, , drop = FALSE]
  sim_row <- cbind(a=a,b=b,times=num_sim)
  # 把两行分别拉成一列（metric 为列名）
  df<-as.data.frame(t(res_i))
  df$Method <-colnames(res_i)
  rownames(df)<-NULL
  colnames(df)<-c("FDR","Power","Method")
  
  # 加上参数列
  df <- cbind(param_row[rep(1, nrow(df)),], df)
  df <- cbind(df,sim_row[rep(1, nrow(df)),])
  return(df)
}))

filename<-"gausig_sup_resultwide.csv"
if(file.exists(filename)){
  res_pr<-read.csv(filename)[,-1]
  wide_result<-rbind(wide_result,res_pr)
  write.csv(wide_result,filename)
}else{
  write.csv(wide_result,filename)
}
