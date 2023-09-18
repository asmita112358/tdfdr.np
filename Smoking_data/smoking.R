
rm(list = ls())
library(tidyr)
library(ggplot2)
library(limma)
library(nspmix)
library(dplyr)
library(splines)
library(doParallel)
library(foreach)
library(tidyr)
library(doMC)
library(pbivnorm)
library(GUniFrac)

###################################################################
###Continuous Y
###################################################################

###Function for the previous method
fastLM <- function(Y, M) {
  Y <- as.matrix(Y)
  XXI <- solve(t(M) %*% M)
  dof <- ncol(Y) - ncol(M)
  est <- XXI %*% t(M) %*% t(Y)
  resid <- t(Y) - M %*% est
  sigma <- sqrt(colSums(resid^2)/dof)
  pval <- 2 * pt(-abs(t(est/(sqrt(diag(XXI)))) / sigma), dof)
  return(list(pval = pval, sigma = sigma, dof = dof, est = est))
}
tdfdr <- function (y,  x,  z,  
                   est.pi0 = TRUE, lambda = 0.5, alpha = level, etype = c('FDR', 'FWER'), 
                   t1a = NULL, t2a = NULL, ngrid = 50, npeb.grid = 500,
                   parallel = FALSE, cores = NULL, verbose = TRUE) {	
  
  etype <- match.arg(etype)
  
  
  if (is.factor(x)) {
    if (nlevels(x) > 2) {
      stop('x should be a vector of two levels!\n')
    } else {
      x <- as.numeric(x)
    }
  }
  
  if (is.vector(z)) {
    z <- as.matrix(z)
  }
  
  if (is.factor(z)) {
    z <- model.matrix(~ z)[, -1, drop = FALSE]
  }
  
  if (is.data.frame(z)) {
    z <- model.matrix(~., data = z)[, -1, drop = FALSE]
  }
  
  if (sum(z[, 1] != 1) == 0) {
    z <- z[, -1, drop = FALSE]
  }
  
  # scaling to improve the stability
  y <- scale(y)
  z <- scale(z)
  x <- as.vector(scale(x))
  
  qrZ <- qr(z, tol = 1e-07)
  z <- qr.Q(qrZ)
  z <- z[, 1:qrZ$rank, drop = FALSE]
  
  n <- length(x)
  p <- ncol(y)
  
  Pz <- diag(n) - z %*% solve(t(z) %*% z) %*% t(z) 
  Px <- diag(n) - x %*% t(x) / sum(x^2) 
  Oz <- as.numeric(t(x) %*% Pz %*% x / n)
  
  O <- sum(x^2) / n
  G <- t(x) %*% z / n
  
  rho <- sqrt(Oz / O)
  
  cat(date(), '\n')
  # We still include intercept for better performance for extremely small sample size
  if (verbose) cat('Perform linear regression ...\n')
  
  M1 <- cbind(1, x)
  M2 <- cbind(1, x, z)
  
  y <- t(y)
  lm.obj1 <- fastLM(y, M1)
  lm.obj2 <- fastLM(y, M2)
  
  coef.x.u <- lm.obj1$est[2, ]
  coef.x.a <- lm.obj2$est[2, ]
  coef.z.a <- lm.obj2$est[-(1:2), , drop = FALSE]
  sigma.est <- lm.obj2$sigma
  
  # need to adjust for multivariate confounder, and collinearity between Z and X
  out <- squeezeVar(sigma.est^2, lm.obj2$dof)
  sigma.est <- sqrt(out$var.post)
  #	
  Zu <- (sqrt(n * O) * coef.x.u / sigma.est)
  Za <- (sqrt(n * Oz) * coef.x.a / sigma.est)
  
  A <- as.numeric(sqrt(1 / O) * G %*% solve(t(z) %*% Px %*% z / n) %*% t(G) * sqrt(1 / O))
  eta <- as.numeric((1 / sqrt(A)) * sqrt(n / O) * G %*% coef.z.a / sigma.est)
  
  if (verbose) cat('Perform NPEB estimation of eta...\n')
  
  #	mm <- REBayes::GLmix(x = eta)
  #	normalized.prob <- mm$y / sum(mm$y)
  
  mm <- cnm(npnorm(eta), grid = npeb.grid)   
  normalized.prob <- mm$mix$pr
  mmx <- mm$mix$pt
  
  # Here pi0 is a vector, and could be a global estimate pi0 <- mean(pi0)
  if (est.pi0 == TRUE) {
    pi0 <- (abs(Za) <= lambda) / (1 - 2 * pnorm(-lambda))
    
  } else {
    pi0 <- rep(1, p)
  }
  
  # Search algorithmn to be optimized in the future
  fdr.est <- function (t1, t2, etype) {
    
    # Number of positives
    NP <- sum(abs(Zu) >= t1 & abs(Za) >= t2)
    
    if (NP != 0) {
      x1 <- -t1 - mmx * sqrt(A)
      x2 <- t1 - mmx * sqrt(A)
      #			x1 <- -t1 - mm$x * sqrt(A)
      #    		x2 <- t1 - mm$x * sqrt(A)
      y1 <- -t2
      y2 <- t2
      
      A1 <- pbivnorm(x = x1, y = y1, rho = rho)
      A2 <- pbivnorm(x = x2, y = y1, rho = rho)
      A3 <- pbivnorm(x = x2, y = y2, rho = rho)
      A4 <- pbivnorm(x = x1, y = y2, rho = rho)
      
      B1 <- pnorm(x1)
      B2 <- pnorm(x2)
      C1 <- pnorm(y1)
      C2 <- pnorm(y2)
      
      if (etype == 'FDR') {
        FD <- sum(pi0) * sum(normalized.prob * (1 + A1 + B1 + C1 + A3 - A2 - B2 - C2 - A4))
        FDP <- FD / NP
      }
      
      if (etype == 'FWER') {
        FDP <- 1 - exp(sum(log(-(A1 + B1 + C1 + A3 - A2 - B2 - C2 - A4)) + log(normalized.prob * length(normalized.prob))
        ))
      }
    } else {
      FDP <- 0
    }
    
    return(list(FDP = FDP, NP = NP))	
  }
  
  if (verbose) cat('Search the best combination of thresholds ...\n')
  
  if (is.null(t1a)) {
    t1a <- seq(min(abs(Zu)), max(abs(Zu)), len = ngrid)
  }
  
  if (is.null(t2a)) {
    t2a <- seq(min(abs(Za)), max(abs(Za)), len = ngrid)
  }
  
  
  if(isTRUE(parallel)){
    
    registerDoMC(cores = cores)
    rr <- foreach(t1 = rep(t1a, each = length(t2a)), t2 = rep(t2a, length(t1a)), .combine = cbind) %dopar% fdr.est(t1, t2, etype)   
    FDP <- unname(unlist(rr[1, ]))
    NP <- unname(unlist(rr[2, ]))
    t1t2 <- paste(rep(t1a, each = length(t2a)), t2a)
    
  } else{
    
    FDP <- NULL
    NP <- NULL
    t1t2 <- NULL
    
    for (t1 in t1a) {
      for (t2 in t2a) {
        obj <- fdr.est(t1, t2, etype)
        NP <- c(NP, obj$NP)
        FDP <- c(FDP, obj$FDP)
        t1t2 <- c(t1t2, paste(t1, t2))
      }
    }
    
  }
  
  
  ind <- which(FDP <= alpha)
  
  if (length(ind) == 0) {
    pos <- rep(FALSE, p)
  } else {
    NP2 <- NP[ind]
    FDP2 <- FDP[ind]
    names(NP2) <- names(FDP2) <- paste(ind)
    
    ind2 <- as.numeric(names(which.min(FDP2[NP2 == max(NP2)])))
    #		ind2 <- as.numeric(names(which.max(NP2)))
    temp <- unlist(strsplit(t1t2[ind2], ' '))
    t1 <- as.numeric(temp[1])
    t2 <- as.numeric(temp[2])
    
    pos <- abs(Zu) >= t1 & abs(Za) >= t2
  }
  
  
  if (verbose) cat('Completed!\n')
  
  return(list(pos = pos, t1 = t1, t2 = t2, t1t2 = t1t2,
              Zu = Zu, Za = Za, NP = NP, FDP = FDP,
              p.value = lm.obj2$pval[, 2], mm = mm, pi0 = mean(pi0)))
}



load("/Users/asmitaroy/Documents/OneDrive - Texas A&M University/Documents/2dFDR/Real_data/Throatdata.rdata")
library(compositions)

X = throat.meta$SmokingStatus
Z = throat.meta$Sex

##Filter out the OTUs with a large number of zeros
index_mat = rowSums(throat.otu.tab == 0) < 54
Y_new = throat.otu.tab[index_mat,]
Y = as.matrix(clr(Y_new + 0.5))
require(sgof)

m = nrow(Y)
n = ncol(Y)
B = 500
if(is.factor(X))
{
  propensity = glm(as.factor(X)~as.factor(Z), family = binomial(link = "logit"))$fitted.values
  Xb = matrix(rbinom(n*B,1,propensity), nrow = n, ncol = B) #maybe as.factor Xb too? check.
  
}
T1.t <- T2.t <- T1.cor <- T2.cor <- vector()
T1b.t <- T2b.t <- T1b.cor <- T2b.cor <- matrix(NA, nrow = m , ncol = B) 
for(j in 1:m)
  
{
  y = Y[j,]
  resY <- scale(resid(lm(y ~ as.factor(Z)))) 
  resX <- scale(as.numeric(X) - 1 - propensity)  ##(X - pi)
  resXb <- scale(Xb - propensity)
  T1.cor[j] <- abs(mean(resY * resX))
  T2.cor[j] <- abs(cor(as.numeric(X) - 1, y))
  T1.t[j] = abs(summary(lm(y ~ as.factor(X)))[["coefficients"]][,"t value"][2])
  T2.t[j] = abs(summary(lm(y ~ as.factor(X) + as.factor(Z)))[["coefficients"]][,"t value"][2])
  for(b in 1:B)
  {
    T1b.t[j,b] = abs(summary(lm(y ~ as.factor(Xb[,b])))[["coefficients"]][,"t value"][2])
    T2b.t[j,b] = abs(summary(lm(y ~ as.factor(Xb[,b]) + as.factor(Z)))[["coefficients"]][,"t value"][2])
  }
  
  T1b.cor[j,] = abs(colMeans(as.vector(resY) * resXb))
  T2b.cor[j,] = abs(cor(Xb,y))
  
}

T1b.cor <- cbind(T1b.cor, T1.cor)
T2b.cor <- cbind(T2b.cor, T2.cor)

T1b.t <- cbind(T1b.t, T1.t)
T2b.t <- cbind(T2b.t, T2.t)


##BH procedure
rejection_fun = function(Y, X, Z, level = level)
{
  p = vector()
  for(i in 1:nrow(Y))
  {
    p[i] = summary(lm(Y[i,] ~ as.factor(X)+as.factor(Z)))[["coefficients"]][2,4]
  }
  benj = BH(p, alpha = level)
  ret0 = benj$Adjusted.pvalues < level
  
  #1DFDR
  
  fdp.hsic <- rej.hsic <- fdp.cor <- rej.cor <- rep(NA, m)
  fdp.hsic.1d <- rej.hsic.1d <- fdp.cor.1d <- rej.cor.1d <- fdp.hsic.mar <- fdp.cor.mar <- rej.hsic.mar <- rej.cor.mar <- rep(NA, m)
  
  for(k in 1:m)
  { 
    t1.h <- T1.t[k]
    t2.h <- T2.t[k]
    t1.c <- T1.cor[k]
    t2.c <- T2.cor[k]
    
    
    Num.1d.h <- sum(T1b.t > t1.h)/(B+1)
    Dem.1d.h <- sum(T1.t > t1.h)
    fdp.hsic.1d[k] <- Num.1d.h/max(1,Dem.1d.h)
    rej.hsic.1d[k] <- Dem.1d.h
    
    Num.1d.c <- sum(T1b.cor > t1.c)/(B+1)
    Dem.1d.c <- sum(T1.cor > t1.c)
    fdp.cor.1d[k] <- Num.1d.c/max(1,Dem.1d.c)
    rej.cor.1d[k] <- Dem.1d.c
    
    
    Num.mar.h = sum(T2b.t > t2.h)/ (B+1)
    Dem.mar.h = sum(T2.t > t2.h)
    fdp.hsic.mar[k] = Num.mar.h/max(1, Dem.mar.h)
    rej.hsic.mar[k] = Dem.mar.h
    
    Num.mar.c = sum(T2b.cor > t2.c)/(B+1)
    Dem.mar.c = sum(T2.cor > t2.c)
    fdp.cor.mar[k] = Num.mar.c/max(1, Dem.mar.c)
    rej.cor.mar[k] = Dem.mar.c
    
    
  }
  
  if(sum(fdp.hsic.1d < level)>0) 
  {
    ind <- which(rej.hsic.1d == max(rej.hsic.1d[fdp.hsic.1d < level]))
    ind.h <- ind[which.min(fdp.hsic.1d[ind])]
    ret.hsic.1d <- (T1.t > T1.t[ind.h])
     
  } else {ret.hsic.1d <- rep(0,m); ind.h = which.max(T1.t)[1]} 
  
  if(sum(fdp.cor.1d < level)>0) 
  {
    ind <- which(rej.cor.1d == max(rej.cor.1d[fdp.cor.1d < level]))
    ind.c <- ind[which.min(fdp.cor.1d[ind])]
    ret.cor.1d <- (T1.cor > T1.cor[ind.c])
    
  } else {ret.cor.1d <- rep(0,m); ind.c = which.max(T1.cor)[1]} 
  
  
  
  if(sum(fdp.hsic.mar < level) >0)
  {
    ind <- which(rej.hsic.mar == max(rej.hsic.mar[fdp.hsic.mar < level]))
    ind.h.m <- ind[which.min(fdp.hsic.mar[ind])]
    
    
  }else {ind.h.m = which(T2.t == max(T2.t))[1]}
  if(sum(fdp.cor.mar < level) >0)
  {
    ind <- which(rej.cor.mar == max(rej.cor.mar[fdp.cor.mar < level]))
    ind.c.m <- ind[which.min(fdp.cor.mar[ind])]
    
    
  }else {ind.c.m = which(T2.cor == max(T2.cor))[1]}
  
  
  ##Power is coming very bad. Investigate. Tomorrow. Ar parchi na. Booh.
  
  ##########################################################################################
  ##Grid search with the two statistics.
  
  
  
  t1t <- quantile(T1.t[T1.t <= T1.t[ind.h]], probs = seq(0, 1, length.out = ngrid))
  
  t2t <- quantile(T2.t[T2.t <= T2.t[ind.h.m]], probs = seq(0, 1, length.out = ngrid))
  
  t1c <- quantile(T1.cor[T1.cor <= T1.cor[ind.c]], probs = seq(0, 1, length.out = ngrid))
  
  t2c <- quantile(T2.cor[T2.cor <= T2.cor[ind.c.m]], probs = seq(0, 1, length.out = ngrid))
  
  
  rej.cor <- fdp.cor <- rej.reg <-fdp.reg <- matrix(nrow = ngrid, ncol = ngrid)
  
  G.cor <- G.reg <- matrix(0, nrow = ngrid, ncol = ngrid)
  
  ##Estimating null proportion
  ##=============================
  
  
  lambda = mean(T1.cor)
  pi0 = min(1, mean(T1.cor <= lambda )/ mean(T1b.cor <= lambda))
  
  
  ##grid search
  
  for(i in 1:ngrid)
  {
    for(j in 1:ngrid)
    {
      t1.cor = t1c[i]
      t2.cor = t2c[j]
      t1.reg = t1t[i]
      t2.reg = t2t[j]
      
      num.cor = mean((T1b.cor > t1.cor)*(T2b.cor > t2.cor))
      den.cor = mean((T1.cor > t1.cor)*(T2.cor > t2.cor))
      fdp.cor[i,j] = pi0 * num.cor/max(1/m, den.cor)
      rej.cor[i,j] = den.cor
      if(fdp.cor[i,j] <= level){
        G.cor[i,j] = 1
      }
      
      num.reg = mean((T1b.t > t1.reg)*(T2b.t > t2.reg))
      den.reg = mean((T1.t > t1.reg)*(T2.t > t2.reg))
      fdp.reg[i,j] = pi0 * num.reg/max(1/m, den.reg)
      rej.reg[i,j] = den.reg
      if(fdp.reg[i,j] <= level){
        G.reg[i,j] = 1
      }
    }
  }
  
  if(sum(G.reg) > 0){
    
    temp1 = rej.reg * G.reg
    
    if(sum(temp1) > 0){
      t_star_ind_reg = which(temp1 == max(temp1), arr.ind = T)[1,] ##run until this
      t1.reg.star = t1t[t_star_ind_reg[1]]  ##Check is t1 is rows and t2 is column
      t2.reg.star = t2t[t_star_ind_reg[2]]
      
      ret.ms.2d = (T1.t > t1.reg.star)*(T2.t > t2.reg.star)
     
    }else{ret.ms.2d <- rep(0,m)}
    
    
    
    
    
  }else{ret.ms.2d <- rep(0,m)}
  
  if(sum(G.cor) >  0){
    temp2 = rej.cor * G.cor
    if(sum(temp2) > 0){
      t_star_ind_cor = which(temp2 == max(temp2), arr.ind = T)[1,] ##run until this
      t1.cor.star = t1c[t_star_ind_cor[1]]  ##Check is t1 is rows and t2 is column
      t2.cor.star = t2c[t_star_ind_cor[2]]
      
      ret.cor.2d = (T1.cor > t1.cor.star)*(T2.cor > t2.cor.star)
      
    }else{ret.cor.2d <- rep(0,m)}
    
    
  }else{ret.cor.2d <- rep(0,m)}
  
  
  obj = tdfdr(t(Y), X, Z, alpha = level)
  ret3 =  obj$pos 
  #return(rbind(ret0, ret.cor.1d, ret.hsic.1d, ret.ms.2d, ret.cor.2d, ret3))
  return(c(sum(ret0), sum(ret.cor.1d), sum(ret.hsic.1d), sum(ret.ms.2d), sum(ret.cor.2d), sum(ret3)))
  #return(list(benjamini = ret0, mfreg = sum(ret), mfcor = sum(ret2), mfreg1d = sum(ret.1d), mfcor1d = sum(ret2.1d), td = sum(ret3)))
  
   
}
u = rejection_fun(Y,as.factor(X),as.factor(Z), 0.1)
u = ifelse(u ==0, FALSE, TRUE)

data = as.data.frame(t(u[c(1,3,4,6),]))


value = 1:174
library(ggvenn)
data = data.frame(value, data)
colnames(data) = c("value", "BH", "1dFDR-MS", "MF-2dFDR-MS", "2dFDR")
ggvenn(data, c("BH", "1dFDR-MS", "MF-2dFDR-MS", "2dFDR"), fill_color = c("green4", "goldenrod", "goldenrod4", "gray42"))


rejectionmat2 = matrix(nrow = 11, ncol = 6)
level_vec = seq(0.01, 0.21, by = 0.02)
for( i in 1:11)
{
  rejectionmat2[i,] = rejection_fun(Y,as.factor(X),as.factor(Z), level = level_vec[i])
  print(i)
}

#alpha = (1:10)/20

rej = cbind(level_vec, rejectionmat2)
colnames(rej) = c("alpha", "BH",  "RV-1dFDR", "MS-1dFDR","MS-2dFDR+", "RV-2dFDR+", "2DFDR")
write.csv(rej, "smoking_cont.csv")
rej = read.csv("smoking_cont.csv")
rej = as.data.frame(rej)
rej = rej[,-1]
colnames(rej) = c("alpha", "BH",  "RV-1dFDR", "MS-1dFDR","MS-2dFDR+", "RV-2dFDR+", "2DFDR")

library(reshape2)
library(ggplot2)
plot_data = melt(rej, id.var = "alpha")
colnames(plot_data) = c("alpha", "method", "value")


ggplot(plot_data, aes(x = alpha, y = value, group = method, color = method)) +
  geom_point(aes( shape = method, color = method)) + 
  # scale_y_continuous(limits=c(0,30))+
 geom_line() + 
  theme_bw() + 
  # scale_shape_manual("method", values = c("RV-1dFDR" = 21, "RV-2dFDR+" = 22, "MS-1dFDR" = 23, "MS-2dFDR+" = 24, "2dFDR" = 25, "BH" = 10))+
  scale_color_manual("method", label = c("BH",  "RV-1dFDR", "MS-1dFDR", "RV-2dFDR+", "MS-2dFDR+","2dFDR"), values = c("green4","skyblue1","goldenrod", "dodgerblue4", "goldenrod4", "gray42 " ))+ 
   scale_shape_manual("method", label = c("BH",  "RV-1dFDR", "MS-1dFDR", "RV-2dFDR+", "MS-2dFDR+","2dFDR"), values = c(15:18,3,7))+
  
  labs(x = "FDR", y = "Rejections")


