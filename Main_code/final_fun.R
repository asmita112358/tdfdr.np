library(dplyr)
library(splines)
library(doParallel)
library(foreach)
library(tidyr)
library(MASS)
library(nnet)
library(car)
library(fastDummies)
###sample functions for computing the two statistics

K_matrix = function(X, h)  ###function for Kernel matrix
{
  
  return(exp(-as.matrix(dist(X))^2/(2*h^2)))
  
}

##function for computing multivariate RV
stat = function(ex, ey)
{
  cx = cov(as.matrix(ex))
  cy = cov(as.matrix(ey))
  cxy = cov(as.matrix(ex), as.matrix(ey))
  return(sum(cxy^2)/(sqrt(sum(cx^2)*sum(cy^2))))
}

###function for computing multivariate RV
mult.rv = function(X, Y, Z)
{
  ey = resid(lm(Y~Z))
  ex = resid(lm(X~Z))
  t1 = stat(ey, ex)
  t2 = stat(X,Y)
  return(c(t1,t2))
}

##Takes X, Y and Z and returns the hsic, chsic
##X, Y and Z are vectors of length n (sample size) if they are univariate, 
#X, Y and Z are matrix with n rows and p columns if multi(p-variate)variate

stat.hsic.fun = function(X, Y, Z) 
{
  n = ifelse(is.null(dim(X)), length(X), nrow(X))
  H = diag(rep(1,n)) - 1/n * matrix(rep(1, n^2), nrow = n, ncol = n)             
  hx = median(as.matrix(dist(X)))
  KX = H %*% K_matrix(X, hx) %*% H
  hz = median(as.matrix(dist(Z)))                    
  hxz = median(as.matrix(dist(cbind(X,Z))))
  KZ = H%*% K_matrix(Z, hz) %*% H
  R = 0.001 * solve(KZ + 0.001*diag(rep(1,n)))
  XZ = cbind(X,Z)
  KXZ = R %*% H%*% K_matrix(XZ, hxz)%*%H %*% R
  hy = median(as.matrix(dist(Y)))
  KY = H %*% K_matrix(Y, hy) %*% H
  KYZ = R %*% KY %*% R
  
  return(c(cond = (1/n)*sum(KXZ*KYZ), mar = (1/n)*sum(KX* KY)))
  
  
  
}

stat.rv.fun = function(X, Y, Z,dof =5, response = "continuous", covariate = "binary")  
{
  if(is.factor(Z) == FALSE && is.null(dim(Z))) BZ = bs(Z, df=dof, intercept = FALSE)
  if(is.factor(Z)) BZ =Z
  if(is.null(dim(Z)) == FALSE) BZ = Z
  
  if(response == "poisson") {
   
   propensity = glm(Y~BZ, family = poisson)$fitted.values
   resY = scale(Y - propensity)
   
 } else if(response == "nb") {
   
   propensity = glm.nb(Y ~ BZ)$fitted.values
   resY = scale(Y - propensity)
   
   
 } else if(response == "continuous"){
   
   resY = scale(resid(lm(Y ~ BZ))) 
   
 } else if(response == "binary"){
   
   propensity = glm(Y~ BZ, family = binomial(link = "logit"))$fitted.values
   resY = scale(as.numeric(Y) - propensity)
   
 }else stop(print("Error: out of scope of this function, please use the stat.fun argument to define the statistics"))
  
  
  
  
  ###Finding resX
  if(covariate == "poisson") {
    
    propensity = glm(X ~ BZ, family = poisson)$fitted.values
    resX = scale(X - propensity)
    
  } else if(covariate == "nb") {
    
    if(is.null(dim(X)) == FALSE){
      stop(print("Error: Currently not supporting multivariate Negative Binomial X"))
    }
    propensity = glm.nb(X ~ Z)$fitted.values
    resX = scale(X - propensity)
    
    
  } else if(covariate == "continuous"){
    
      resX = scale(resid(lm(X ~ BZ))) 
    
  } else if(covariate == "binary"){
    
    propensity = glm(X ~ BZ, family = binomial(link = "logit"))$fitted.values
    resX = scale(as.numeric(X) - propensity) 
    
  }else stop(print("Error: out of scope of this function, please use the stat.fun argument to define the statistics"))
  
  if( ncol(resX) == 1 ){
    t1 = abs(mean(resY * resX))
    t2 = abs(cor(as.numeric(X),Y))
  }else{
    t1 = stat(resX, resY)
      t2 = stat(X, Y)
    }
  
 
  return(c(t1, t2))
  
}

##model can be "linear", "glm"
##family can be "binomial", "poisson" or "nb"
##link would be whatever link you'd like to use in the model. eg "log", "logit" etc.
stat.ms.fun = function(X, Y, Z,dof = 5, response = "continuous", covariate = "multinomial",link = "identity") ##takes X, Y and Z and returns the model based statistics
{
  if(is.factor(Z) == FALSE && is.null(dim(Z))) BZ = bs(Z, df=dof, intercept = FALSE)
  if(is.factor(Z)) BZ =Z
  if(is.null(dim(Z)) == FALSE) BZ = Z
  if(is.null(dim(X)) == FALSE){
    q = ncol(X)
  }else if(is.factor(X))
  {
    q = nlevels(X)
  }else{
    q = 1
  }
 if(response == "continuous"){
      fit2 = lm(Y~X)
    fit1 = lm(Y~X+BZ)
      t2 = anova(fit2)$F[1]
      t1 = anova(fit1)$F[1]
  
      
 } else if(response == "binary"){
  # if(is.factor(X)){
    # fit1 = glm(Y ~   X + BZ, family = binomial(link = link))
    # fit2 = glm(Y ~ X, family = binomial(link = link))
    # t1 = qnorm(1-anova(fit1,test="LRT")[[5]][3])
    # t2 = qnorm(1-anova(fit2,test="LRT")[[5]][2])

  # }else{
     t1 = abs(sum(summary(glm(Y ~   X + BZ, family = binomial(link = link)))[["coefficients"]][,"z value"][2:(q +1)]))
     t2 = abs(sum(summary(glm(Y ~ X, family = binomial(link = link)))[["coefficients"]][,"z value"][2:(q+1)]))
     
#   }
        
  
     
   }else if(response == "poisson"){
    
      t1 = abs(summary(glm(Y ~  X + BZ, family = poisson(link = link), control = glm.control(maxit = 100)))[["coefficients"]][,"z value"][2:(q+1)] )
      t2 = abs(summary(glm(Y ~ X, family = poisson(link = link), control = glm.control(maxit = 100)))[["coefficients"]][,"z value"][2:(q+1)])
       
    
    
       
     
     
   }else if(response == "nb"){
     if(q > 1){
     stop(print("Error: Currently negative binomial response is only compatible with univariate, continuous X or binary X"))
     }else{
       t1 = abs(summary(glm.nb(Y ~   X + BZ, control = glm.control(maxit = 500)))[["coefficients"]][,"z value"][2])
       t2 =  abs(summary(glm.nb(Y ~ X, control = glm.control(maxit = 100)))[["coefficients"]][,"z value"][2])
       
       }
     
       
   
      
     }else stop(print("Error: out of scope of this function, please use the stat.fun argument to define the statistics"))

   return(c(t1,t2))
}

###Sample functions for generating data from the conditional distribution

generate.cont = function(X, Z, B) ###Use this if X and Z are continuous
{
  
  if(is.factor(X))return("Error: use the function for generating discrete data")
  ###what happens if X is multivariate? think. Same for Z
  if(is.null(dim(X)))
  {
    n = length(X)
    Xb = matrix(nrow = n, ncol = B)
    Lm = lm(X ~ Z)
    res = resid(Lm)
    index = sapply(rep(n, B), function(x) sample(x))
    Xb  = matrix(res[index], n, B) + Lm$fitted.values
    return(Xb)
  }else{
    p = ncol(Z)
    q = ncol(X)
    n = nrow(X)
    mlm = lm(X ~ Z)
    res = resid(mlm)
    index = sapply(rep(n, B), function(x) sample(x))
    Xb = array(dim = c(B, n, q))
    for(b in 1:B)
    {
      Xb[b,,]  = matrix(res[index,], nrow = n, ncol = q) + mlm$fitted.values 
    }
    return(Xb)
  }
  
}



generate.discrete = function(X, Z, B) #Use this if X is discrete
{
  
  if(is.null(dim(X))){
    if(nlevels(X) > 2) {
      q = nlevels(X)
      model = multinom(X~Z)
      prob = fitted(model)
      Xb = array(dim = c(B,n,q))
      Xb_mn = matrix(nrow = n, ncol = B)
      for(b in 1:B)
      {
        temp = vector()
        for(i in 1:n)
        {
          Xb[b,i,] = rmultinom(1,1,prob[i,])
        }
        index = which(Xb[b,,] ==1, arr.ind = T)
        temp[index[,1]] = as.factor(index[,2]) 
        temp = as.factor(temp)
        Xb_mn[,b] = temp
      }
      return(Xb_mn)
    }else{
      n = length(X)
      propensity = glm(X~ Z, family = binomial(link = "logit"))$fitted.values
      Xb = matrix(rbinom(n*B,1,propensity), nrow = n, ncol = B) 
      return(Xb)
    } 
  }else print("Please write your own data generation function")
  
  
}


###function for 2dfdr-np


tdfdr.np = function(X, Y, Z, ngrid = 50, level = 0.05, B = 50, parallel = TRUE, 
                    ncores = detectCores() - 1, method = "RV", response = "continuous", covariate = "continuous",
                    link = "log", etype = "FDR", stat.fun = stat.rv.fun)
{
  vec = dim(Y)
  m = vec[1]
  n = vec[2]
  l = vec[3]
  
  
  ###Checking the dimensions
  if(is.null(dim(X))) {
    check = (length(X) == n)
    q = 1
  }else {
    q = ncol(X)
    check = (nrow(X) == n)
  }
  if(check == FALSE) stop("Error: X and Y dimensions don't match")
  if(is.null(dim(Z))) {
    check = (length(Z) == n)
    p = 1
  } else {
    p = ncol(Z)
    check = (nrow(Z) == n)
  }
  if(check == FALSE) stop("Error: Z and Y dimensions don't match")
  if(covariate == "multinomial" && method == "RV") stop("Error: Please use method = MS")
  if(covariate == "multinomial" && response != "continuous") stop("Error: Currently only supporting continuous response with multinomial covariates")
  ####
  if(is.factor(X) == FALSE && (covariate == "binary"  || covariate == "multinomial")) stop("Please put X as as.factor()")
  if((covariate != "continuous" || response != "continuous" ||is.factor(Z)) && method == "HSIC") stop(print("HSIC can only be used for continuous X and Y"))
  ###Generating Xb
  if(is.factor(X)){
    Xb = generate.discrete(X, Z, B)
  }else{
    Xb = generate.cont(X, Z, B)
  }
  
  
  
  ###Computing the statistics
  T1 = T2 = vector()
  T1b = T2b = matrix(nrow = m, ncol = B)
  cat('Computing Statistics...\n')
  output = sapply(mclapply(1:m, function(j) {
   # for(j in 1:m)
     # {
    if(is.na(l)){
      y = Y[j,]
    }else{
      stop("Currently only supporting univariate Y")
    }
    if(method == "RV"){
      obj = stat.rv.fun(X, y, Z, dof = 5, response = response, covariate = covariate)
    }else if(method == "HSIC"){
      obj = stat.hsic.fun(X, y, Z)
    }else if(method == "MS"){
      obj = stat.ms.fun(X, y, Z, response = response, covariate = covariate, link = link)
    }else {
      obj = stat.fun(X, y, Z)
    }
    
    T1[j] = obj[1]
    T2[j] = obj[2]
    for(b in 1:B)
    {
      
      if(q == 1){
        
        xb_gen = Xb[,b]
      }else{
        
        xb_gen = Xb[b,,]
      }
      
      if(method == "RV"){
        obj = stat.rv.fun(xb_gen, y, Z, dof = 5, response = response, covariate = covariate)
      }else if(method == "HSIC"){
        obj = stat.hsic.fun(xb_gen, y, Z)
      }else if(method == "MS"){
        obj = stat.ms.fun(xb_gen, y, Z, response = response, covariate = covariate, link = link)
      }else {
        obj = stat.fun(xb_gen, y, Z)
      }
      T1b[j,b] = obj[1]
      T2b[j,b] = obj[2]
    }
    u = c(T1b[j,], T1[j])
    v = c(T2b[j,], T2[j])
    #print(j)
    #}
    return(c(u, v))
  }, mc.cores = ncores), 'c')
  
  T1b = t(output[1:B, ])
  T2b = t(output[(B+2):(2*B + 1), ])
  
  T1 = output[B+1, ]
  T2 = output[2*B + 2, ]
  cat('Searching for thresholds in each dimension...\n')
  ##1dFDR
  fdp.1d = rej.1d = fdp.mar = rej.mar = rep(NA,m)
  fdp = rej = rep(NA,m)
  
  for(k in 1:m)
  {
    t1.c = T1[k]
    t2.c = T2[k]
    
    Num.1d = sum(T1b > t1.c)/(B+1)
    Dem.1d = sum(T1 > t1.c)
    fdp.1d[k] = Num.1d/max(1,Dem.1d)
    rej.1d[k] = Dem.1d
    
    Num.mar = sum(T2b > t2.c)/(B+1)
    Dem.mar = sum(T2 > t2.c)
    fdp.mar[k] = Num.mar/max(1, Dem.mar)
    rej.mar[k] = Dem.mar
    
  }
  if(sum(fdp.mar < level) >0)
  {
    ind = which(rej.mar == max(rej.mar[fdp.mar < level]))
    ind.m = ind[which.min(fdp.mar[ind])]
    
    
  }else {ind.m = which(T2 == max(T2))[1]}  
  
  if(sum(fdp.1d < level)>0) 
  {
    ind = which(rej.1d == max(rej.1d[fdp.1d < level]))
    ind.c = ind[which.min(fdp.1d[ind])]
    
  } else {ind.c = which.max(T1)[1]} 
  
  ##MF-2dFDR  ##Gridsearch
  cat('2D Gridsearch...\n')
  t1c = quantile(T1[T1 <= T1[ind.c]], probs = seq(0, 1, length.out = ngrid))
  t2c = quantile(T2[T2 <= T2[ind.m]], probs = seq(0, 1, length.out = ngrid))
  
  
  rej.c = matrix(0, nrow = ngrid, ncol = ngrid)
  fdp.c = matrix(0, nrow = ngrid, ncol = ngrid)
  G.c = matrix(0, nrow = ngrid, ncol = ngrid)
  
  
  ##Estimating null proportion
  ##=============================
  
  
  lambda = mean(T1)
  pi0 = min(1, mean(T1 <= lambda )/ mean(T1b <= lambda))
  
  ##Gridsearch
  
  for(i in 1:ngrid)
  {
    for(j in 1:ngrid)
    {
      t1.c = t1c[i]
      t2.c = t2c[j]
      
      
      
      num = sum((T1b > t1.c)*(T2b > t2.c))/(B+1)
      den = sum((T1 > t1.c)*(T2 > t2.c))
      
      if(etype == "FDR"){
        fdp.c[i,j] = pi0 * num/max(1, den)
        rej.c[i,j] = den
      }else if(etype == "FWER"){
        fdp.c[i,j] = num
        rej.c[i,j] = den
      }else{
        stop(print("etype must be either of FDR or FWER"))
      }
      
      if(fdp.c[i,j] <= level) G.c[i,j] = 1
      
     
    }
  }
  if(sum(G.c) > 0){
    
    
    temp2 = rej.c * G.c
    if(sum(temp2) > 0) {
      t_star_ind_c = which(temp2 == max(temp2), arr.ind = T)[1,] 
      t1.c.star = t1c[t_star_ind_c[1]]  
      t2.c.star = t2c[t_star_ind_c[2]]
      ret2 = (T1 > t1.c.star)*(T2> t2.c.star)
      
    } else {stop(print("Cutoffs not found1"))}
    
    
  } else {stop(print("Cutoffs not found2")) }
  cat("Done!")
  return(list(cutoffs = c(t1.c.star, t2.c.star), rejections = ret2))
  
}