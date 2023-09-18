###Cat X, cat Y, cat Z with pearsons chi sq and mh coefficient.
###Source the last para of this code for a line diagram on the number of rejections.
rm(list = ls())
set.seed(109)
require(sgof)
library(GUniFrac)
require(dplyr)
require(foreach)
library(doParallel)
library(doRNG)
load("/Users/asmitaroy/Documents/OneDrive - Texas A&M University/Documents/2dFDR/Real_data/Throatdata.rdata")
library(compositions)

X = throat.meta$SmokingStatus
Z = throat.meta$Sex

##Filter out the OTUs with a large number of zeros
index_mat = rowSums(throat.otu.tab == 0) < 54
Y_new = throat.otu.tab[index_mat,]



Y = Rarefy(t(Y_new))$otu.tab.rff
Y = as.matrix((Y>0) + 0)
Y = t(Y)
Y = Y[rowSums(Y) > 0 && rowSums(Y) < 60,]
dim(Y)


  m = nrow(Y)
  n = ncol(Y)
  
 
  ##Generating additional Xb
  B = 500
  if(is.factor(X))
  {
    propensity = glm(as.factor(X)~as.factor(Z), family = binomial(link = "logit"))$fitted.values
    Xb = matrix(rbinom(n*B,1,propensity), nrow = n, ncol = B) #maybe as.factor Xb too? check.
    
  }
  
  ##Computing the statistic
  
  T1b <- T2b <- matrix(NA, nrow = m, ncol = B)
  T1 <- T2 <- vector()
  
  for(j in 1:m)
  {
    print(j)
    y = Y[j,]
    T2[j] = chisq.test(as.factor(y), X)$statistic
    T1[j] = mantelhaen.test(as.factor(y), X, Z)$statistic
    
    for(b in 1:B)
    {
      T2b[j, b] = chisq.test(as.factor(y), as.factor(Xb[,b]))$statistic
      T1b[j, b] = mantelhaen.test(as.factor(y), as.factor(Xb[,b]), Z)$statistic
    }
    
  }
  T1b <- cbind(T1b, as.vector(unlist(T1)))
  T2b <- cbind(T2b, as.vector(unlist(T2)))
  #####Calculating the null proportion
  lambda = mean(T1)
  
  pi0= min(1, mean((T1 <= lambda)/ if_else(rep(sum(rowMeans(T1b <= lambda)), m) > rep(0, m), rowMeans(T1b <= lambda), rep(1,m))))
  
  
  ########################################################################################
  ##Grid search with the two statistics.
  
  ngrid = 200
  
  rejection_fun_smoking = function(X, Y, Z, T1, T2, T1b, T2b, level)
  {
    
    p = vector()
    for(i in 1: nrow(Y))
    {
      p[i] = mantelhaen.test(x = as.factor(Y[i,]), y = X, z = Z, conf.level = 1-level)$p.value
      
    }
    benj = BH(p, alpha = level)
    ret0 = benj$Rejections
    ##1DFDR
    #####1DFDR
    
    rej.1d <- fdp.1d <- fdp.mar <- rej.mar <-  rep(NA,m)
    #T2 = as.vector(unlist(T2))
    for(k in 1:m)
    {
      Num.1d<- sum(T1b > T1[k])/(B+1)
      Dem.1d <- sum(T1 > T1[k])
      fdp.1d[k] <- Num.1d/max(1,Dem.1d)
      rej.1d[k] <- Dem.1d
      
      Num.mar <- sum(T2b > T2[k])/(B+1)
      Dem.mar <- sum(T2 > T2[k])
      fdp.mar[k] <- Num.mar/max(1, Dem.mar)
      rej.mar[k] <- Dem.mar
        
    }
    
    
    if(sum(fdp.1d < level)>0) 
    {
      ind <- which(rej.1d == max(rej.1d[fdp.1d < level]))
      ind.s <- ind[which.min(fdp.1d[ind])]
      ret.1d <- (T1 > T1[ind.s])
      
    } else {ret.1d = rep(0, m); ind.s = which.max(T1)[1]} 
    
    if(sum(fdp.mar < level) >0)
    {
      ind <- which(rej.mar == max(rej.mar[fdp.mar < level]))
      ind.h.m <- ind[which.min(fdp.mar[ind])]
      
      
    }else {ind.h.m = which(T2 == max(T2))[1]}
    ##Define the grids
    T1 = as.vector(unlist(T1))
    t1t <- quantile(T1[T1<= T1[ind.s]], probs = seq(0, 1, length.out = ngrid))
    
    t2t <- quantile(T2[T2<= T2[ind.h.m]], probs = seq(0, 1, length.out = ngrid))
    
  rej <- fdp <- matrix(nrow = ngrid, ncol = ngrid)
  
  G <- matrix(0, nrow = ngrid, ncol = ngrid)
  
  for(i in 1:ngrid)
  {
    for(j in 1:ngrid)
    {
      num.cor = mean((T1b > t1t[i])*(T2b > t2t[j]))
      den.cor = mean((T1 > t1t[i])*(T2 > t2t[j]))
      fdp[i,j] = pi0 * num.cor/max(1/m, den.cor)
      rej[i,j] = den.cor
      if(fdp[i,j] <= level){
        G[i,j] = 1
      }
      
    }
  
  }
  
  
  if(sum(G) > 0){
    
    temp1 = rej * G
    
    if(sum(temp1) > 0){
      t_star_ind = which(temp1 == max(temp1), arr.ind = T)[1,] ##run until this
      t1.star = t1t[t_star_ind[1]]  ##Check is t1 is rows and t2 is column
      t2.star = t2t[t_star_ind[2]]
      
      ret = (T1> t1.star)*(T2 > t2.star)
      
    }else{ret = rep(0, m)}
  }else{ret = rep(0, m)}
  
  
  
  
  return(c(ret0, sum(ret), sum(ret.1d)))
  
  
  }
  
  rejection_fun_smoking(0.2)
rejectionmat3 = matrix(nrow = 11, ncol = 3)
level_vec = seq(0.1, 0.3, by = 0.02)
cl = makeCluster(5)
registerDoParallel(cl)
rejectionmat3 = foreach( i=1:length(level_vec), .combine = rbind, .packages = "sgof") %dopar% {
  rejection_fun_smoking(X, Y, Z, T1, T2, T1b, T2b, level = level_vec[i])
 
}
stopCluster(cl)
write.csv(rejectionmat3, "smoking_mhpearson1.csv")



###############################################################
##Source from here for the graph
###############################################################
rejectionmat3 = read.csv("smoking_mhpearson1.csv")[,-1]
rej = cbind(level_vec, rejectionmat3)
colnames(rej) = c("alpha", "BH", "2dFDR+", "1dFDR")

rej = as.data.frame(rej)
library(reshape2)
library(ggplot2)
plot_data = melt(rej, id.var = "alpha")
colnames(plot_data) = c("alpha", "method", "value")
ggplot(plot_data, aes(x = alpha, y = value, group = method, colour = method)) +
  geom_point(aes( shape = method, color = method)) + 
  geom_line() + theme_bw() + 
  scale_color_manual("method", labels = c("BH", "2dFDR+", "1dFDR"), values = c("green4", "dodgerblue4", "skyblue1"))+ 
  scale_shape_manual("method", labels = c("BH", "2dFDR+", "1dFDR"), values = c(15, 1,2))+
  labs(x = "FDR", y = "Rejections")

