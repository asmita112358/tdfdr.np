###tdfdr
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
library(sgof)
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
                   est.pi0 = TRUE, lambda = 0.5, alpha = 0.05, etype = c('FDR', 'FWER'), 
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



###Metabolome data

filepath = "/Users/asmitaroy/Desktop/TDFDR/RealData/Metabolomics_Figure2bc/Data/"
library(dplyr)
library(textshape)
library(radiant.data)
# loading data
pheno <- read.table(paste0(filepath, "phenotypes.tab"), sep = "\t", head = TRUE, row.names = 1)
head(pheno)


dim(pheno)
metabolic <- read.table(paste0(filepath, 'metabolomic.tab'), sep = '\t', head = TRUE, row.names = 1) 
lipidomic <- read.table(paste0(filepath, 'lipidomic.tab'), sep = '\t', head = TRUE, row.names = 1)
metabolome <- inner_join(metabolic%>% rownames_to_column('id'), lipidomic%>% rownames_to_column('id')) %>% column_to_rownames('id')
dim(metabolome);dim(lipidomic);dim(metabolic)


# filter all NAs, most occurred in Homa.IR
idx <- rownames(pheno[!is.na(pheno$Homa.IR),]);length(idx)
metabolic1 <- metabolic[(rownames(metabolic)%in% idx),] %>% as.matrix();dim(metabolic1);dim(metabolic)
lipidomic1 <- lipidomic[(rownames(lipidomic)%in% idx),] %>% as.matrix();dim(lipidomic1);dim(lipidomic)
metabolome1 <- metabolome[(rownames(metabolome)%in% idx),] %>% as.matrix() 


# distribution plots
par(mfrow=c(1,3))
hist(t(metabolome1));hist(metabolic1);hist(lipidomic1)


# replace 0 with (0.5 *minimal non-zero value)
replace_halfmin <- function(df){
  for(i in 1:nrow(df)){
    idx = which(df[i,] == 0)
    df[i, idx] = min(df[i, df[i,] > 0] * 0.5) 
  }
  return(df)
}

metabolic1 <- replace_halfmin(metabolic1)
lipidomic1 <- replace_halfmin(lipidomic1)
metabolome1 <- replace_halfmin(metabolome1)


# log2 transformation
metabolic1 <- log2(metabolic1);dim(metabolome1)
lipidomic1 <- log2(lipidomic1);dim(metabolome1)
metabolome1 <- log2(metabolome1);dim(metabolome1)
pheno1 <- pheno[!is.na(pheno$Homa.IR),];dim(pheno1)
table(pheno$Diabetes);table(pheno1$Diabetes)
cor.test(pheno1$Homa.IR, pheno1$BMI.kg.m2, method = 'spearman')

par(mfrow=c(1,3))
hist(t(metabolome1));hist(metabolic1);hist(lipidomic1)


#  tdfdr analysis for Insulin resistance metabolomics dataset pooling polar metabolites and molecular lipids


###MF-2DFDR analysis

K_matrix = function(X, h)
{
  
  return(exp(-as.matrix(dist(X))^2/(2*h^2)))
  
}



  Y = t(metabolome1)
  X = pheno1$Homa.IR
  Z = pheno1$BMI.kg.m2
  
  
  m <- nrow(Y)
  n <- ncol(Y)
  if(length(X) != n) return("Error: X and Y dimensions don't match")
  
  # Generate Xb
  dof = 5
  B = 500
  ngrid = 400
  
  ###The bad method
  ###+++++++++++++++++++++++++++++++++++++++++++
  Xb_bad = vector()
  BZ = bs(Z, df = dof, intercept = FALSE)
  obj = lm(X~BZ)
  res = resid(obj)
  index = sample(n)
  Xb_bad = obj$fitted.values + res[index]
  
  ###Alternative generation of X from X|Z: Zhang approved method
  ##++++++++++++++++++++++++++++++++++++++++++++
  Xb = list()
  m1 = mean(X[Z <= 26])
  v1 = var(X[Z <= 26])
  n1 = length(X[Z <= 26])
  BZ = bs(Z[Z<= 26], df=dof, intercept = FALSE)
  lm1 = lm(X[Z <= 26] ~ BZ)
  res1 = resid(lm1)
  index1 = sapply(rep(n1, B), function(x) sample(x))
  Xb[[1]] = lm1$fitted.values + matrix(res1[index1], n1, B)
  
  m2 = mean(X[Z > 26])
  v2 = var(X[ Z > 26])
  n2 = length(X[Z >26])
  BZ2 = bs(Z[Z > 26], df=dof, intercept = FALSE)
  lm2 = lm(X[Z > 26] ~ BZ2)
  res2 = resid(lm2)
  index2 = sapply(rep(n2, B), function(x) sample(x))
  Xb[[2]] = lm2$fitted.values + matrix(res2[index2], n2, B)
  
  Xb_new = do.call(rbind, Xb)
  Z_new = c(Z[Z <= 26], Z[Z> 26])
  X = c(X[Z <= 26], X[Z > 26])
 
  par(mfrow = c(1,3))
  plot(Z_new, X, xlab = "BMI", ylab = "IR", main = "Original Data")
  plot(Z, Xb_bad, xlab = "BMI", ylab = "IR(Resampled)", main = "Residual permutation")
  plot(Z_new, Xb_new[,2], xlab = "BMI", ylab = "IR (Resampled)", main = "Binned Residual permutation")
  cor(X, Z_new)
            ##Does not work
  ##++++++++++++++++++++++++++++++++++++++++++++
  
  
  ##Brute force way, not mathematically sound. 
  ##+++++++++++++++++++++++++++++++++++++++++++
  #Xb = list()
  #m1 = mean(X[Z <= 26])
  #v1 = var(X[Z <= 26])
  #n1 = length(X[Z <= 26])
  #Xb[[1]] = sapply(rep(n1, B), function(x) rgamma(x, 3, 3.5))
  
  #m2 = mean(X[Z > 26])
  #v2 = var(X[ Z > 26])
  #n2 = length(X[Z >26])
  #Xb[[2]] = sapply(rep(n2, B), function(x) rgamma(x, 2, 1))
  
  #Xb_new = do.call(rbind, Xb)
  ###Did not work either
  ##+++++++++++++++++++++++++++++++++++++++++++
  
  # Compute T1.cor, T2.cor, T1.hsic, T2.hsic
  ###Changing Y and X accordingly
  for(i in 1:m)
  {
    Y[i,] = c(Y[i,][Z <= 26], Y[i,][Z > 26])
  }
  X = c(X[Z <= 26], X[Z > 26])
  Z = Z_new
  T1.cor <- T2.cor <- T1.hsic <- T2.hsic <- vector()
  T1b.cor <- T2b.cor <- T1b.hsic <- T2b.hsic <- matrix(NA, nrow = m , ncol = B) 
  BZ = bs(Z_new, df = 5)
  resX <- scale(resid(lm(X ~ BZ)))
  resXb <- scale(resid(lm(Xb_new ~ BZ))) 
  H <- diag(rep(1,n)) - 1/n * matrix(rep(1, n^2), nrow = n, ncol = n)             ##Check the bandwidths
  hx <- median(as.matrix(dist(X)))
  KX <- H %*% K_matrix(X, hx) %*% H
  hz <- median(as.matrix(dist(Z)))                    
  hxz <- median(as.matrix(dist(cbind(X,Z))))
  KZ <- H%*% K_matrix(Z, hz) %*% H
  R <- 0.001 * solve(KZ + 0.001*diag(rep(1,n)))
  XZ <- cbind(X,Z)
  KXZ <- R %*% H%*% K_matrix(XZ, hxz)%*%H %*% R
  
  output <- sapply(mclapply(1:m, function(j){
    
    y <- Y[j,]
    resY <- scale(resid(lm(y ~ BZ))) 
    T1.cor[j] <- abs(mean(resY*resX))
    T2.cor[j] <- abs(cor(X,y))
    T1b.cor[j,] <- abs(colMeans(as.vector(resY) * resXb))
    T2b.cor[j,] <- abs(cor(Xb_new,y))
    
    hy <- median(as.matrix(dist(Y[j,])))
    KY <- H %*% K_matrix(Y[j,], hy) %*% H
    KYZ <- R %*% KY %*% R
    
    T1.hsic[j] <- (1/n)*sum(KXZ*KYZ)
    T2.hsic[j] <- (1/n)*sum(KX* KY)
    
    # mid = Sys.time()
    
    for(b in 1:B)
    {
      hxzb <- median(as.matrix(dist(cbind(Xb_new[,b],Z)))) 
      XZb <- cbind(Xb_new[,b],Z)
      KXZb <-  R %*% H%*% K_matrix(XZb, hxzb)%*%H %*% R
      KXb <- H %*%K_matrix(Xb_new[,b], hx) %*% H
      T1b.hsic[j,b] <- (1/n)*sum(KXZb* KYZ)
      T2b.hsic[j,b] <- (1/n)*sum(KXb*KY)
    }
    
    T1b.hsic <- c(T1b.hsic[j,], T1.hsic[j])
    T2b.hsic <- c(T2b.hsic[j,], T2.hsic[j])
    T1b.cor <- c(T1b.cor[j,], T1.cor[j])
    T2b.cor <- c(T2b.cor[j,], T2.cor[j])
    
    return(c(T1b.hsic, T2b.hsic, T1b.cor, T2b.cor))
  }, mc.cores = 6), 'c') 
  
  T1b.hsic <- t(output[1:B, ])
  T2b.hsic <- t(output[(B+2):(2*B + 1), ])
  T1b.cor <- t(output[(2*B+3):(3*B+2), ])
  T2b.cor <- t(output[(3*B + 4):(4*B + 3), ])
  
  T1.hsic <- output[B+1, ]
  T2.hsic <- output[2*B + 2, ]
  T1.cor <- output[3*B + 3, ]
  T2.cor <- output[4*B + 4, ]
  
  T1b.hsic <-cbind(T1b.hsic, T1.hsic)
  T2b.hsic <-cbind(T2b.hsic, T2.hsic)
  
  T1b.cor <- cbind(T1b.cor, T1.cor)
  T2b.cor <- cbind(T2b.cor, T2.cor)
  
  ###saving the statistics:
  
  write.csv(data.frame(T1b.cor, T2b.cor, T1b.hsic, T2b.hsic), "metabolomics_stat")
  
  
  ####Just running the ngrid code again
  u = read.csv("metabolomics_stat")
  dim(u)
  u = u[,-1]
  T1b.cor = u[,1:501]
  T1.cor = u[,501]
  T2b.cor = u[,502:1002]
  T2.cor = u[,1002]
  T1b.hsic = u[,1003:1503]
  T1.hsic = u[,1503]
  T2b.hsic = u[,1504:2004]
  T2.hsic = u[,2004]
  ##Estimating null proportion
  ##=============================
  
  
  lambda = mean(T1.cor)
  pi0 = min(1, mean(T1.cor <= lambda )/ mean(T1b.cor <= lambda))
  #pi0 = 0.573
  ##=============================
  
  rejection_fun = function( level)
  {
  #old method
  tdfdr.obj <- tdfdr(t(Y), X, Z, alpha = level)
  ret0 = tdfdr.obj$pos
  
  ###Benjamini Hocheberg
  p = vector()
  for(i in 1:m)
  {
    p[i] = summary(lm(Y[i,] ~ X + Z))[["coefficients"]][2,4]
  }
  benj = BH(p, alpha = level)
  ret5 = benj$Adjusted.pvalues < level
  
  
  # 1DFDR
  
  fdp.hsic <- rej.hsic <- fdp.cor <- rej.cor <- rep(NA, m)
  fdp.hsic.1d <- rej.hsic.1d <- fdp.cor.1d <- rej.cor.1d <- fdp.hsic.mar <- fdp.cor.mar <- rej.hsic.mar <- rej.cor.mar <- rep(NA, m)
  
  for(k in 1:m)
  { 
    t1.h <- T1.hsic[k]
    t2.h <- T2.hsic[k]
    t1.c <- T1.cor[k]
    t2.c <- T2.cor[k]
    
    Num.1d.h <- sum(T1b.hsic > t1.h)/(B+1)
    Dem.1d.h <- sum(T1.hsic > t1.h)
    fdp.hsic.1d[k] <- Num.1d.h/max(1,Dem.1d.h)
    rej.hsic.1d[k] <- Dem.1d.h
    
    Num.1d.c <- sum(T1b.cor > t1.c)/(B+1)
    Dem.1d.c <- sum(T1.cor > t1.c)
    fdp.cor.1d[k] <- Num.1d.c/max(1,Dem.1d.c)
    rej.cor.1d[k] <- Dem.1d.c
    
    
    Num.mar.h = sum(T2b.hsic > t2.h)/ (B+1)
    Dem.mar.h = sum(T2.hsic > t2.h)
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
    ret1.hsic <- (T1.hsic > T1.hsic[ind.h])
     
  } else {ret1.hsic = rep(0,m); ind.h = which.max(T1.hsic)[1]} 
  
  if(sum(fdp.cor.1d < level)>0) 
  {
    ind <- which(rej.cor.1d == max(rej.cor.1d[fdp.cor.1d < level]))
    ind.c <- ind[which.min(fdp.cor.1d[ind])]
    ret1.cor <- (T1.cor > T1.cor[ind.c])
    
  } else {ret1.cor = rep(0,m); ind.c = which.max(T1.cor)[1]} 
  
  if(sum(fdp.hsic.mar < level) >0)
  {
    ind <- which(rej.hsic.mar == max(rej.hsic.mar[fdp.hsic.mar < level]))
    ind.h.m <- ind[which.min(fdp.hsic.mar[ind])]
    
    
  }else {ind.h.m = which(T2.hsic == max(T2.hsic))[1]}
  if(sum(fdp.cor.mar < level) >0)
  {
    ind <- which(rej.cor.mar == max(rej.cor.mar[fdp.cor.mar < level]))
    ind.c.m <- ind[which.min(fdp.cor.mar[ind])]
    
    
  }else {ind.c.m = which(T2.cor == max(T2.cor))[1]}
  
  ##MF-2DFDR
  
  # Grid search
  
  t1h <- quantile(T1.hsic[T1.hsic <= T1.hsic[ind.h]], probs = seq(0, 1, length.out = ngrid))
  
  t2h <- quantile(T2.hsic[T2.hsic <= T2.hsic[ind.h.m]], probs = seq(0, 1, length.out = ngrid))
  
  t1c <- quantile(T1.cor[T1.cor <= T1.cor[ind.c]], probs = seq(0, 1, length.out = ngrid))
  
  t2c <- quantile(T2.cor[T2.cor <= T2.cor[ind.c.m]], probs = seq(0, 1, length.out = ngrid))
  
  
  
  
  rej.h <- rej.c <- matrix(0, nrow = ngrid, ncol = ngrid)
  fdp.h <- fdp.c <- matrix(0, nrow = ngrid, ncol = ngrid)
  G.h <- G.c <- matrix(0, nrow = ngrid, ncol = ngrid)
  
  
  
  ##=============================
  
  for(i in 1:ngrid)
  {
    for(j in 1:ngrid)
    {
      t1.c <- t1c[i]
      t2.c <- t2c[j]
      
      t1.h <- t1h[i]
      t2.h <- t2h[j]
      
      num.cor <- sum((T1b.cor > t1.c)*(T2b.cor > t2.c))/(B+1)
      den.cor <- sum((T1.cor > t1.c)*(T2.cor > t2.c))
      
      fdp.c[i,j] <- pi0 * num.cor/max(1, den.cor)
      rej.c[i,j] <- den.cor
      
      if(fdp.c[i,j] <= level) G.c[i,j] = 1
      
      num.hsic <- sum((T1b.hsic > t1.h)*(T2b.hsic > t2.h))/(B+1)
      den.hsic <- sum((T1.hsic > t1.h)*(T2.hsic > t2.h))
      
      fdp.h[i,j] <- pi0 * num.hsic/max(1, den.hsic)
      rej.h[i,j] <- den.hsic
      
      if(fdp.h[i,j] <= level) G.h[i,j] = 1
    }
  }
  
  if(sum(G.h) > 0){
    
    temp1 = rej.h * G.h
    if(sum(temp1) > 0) {
      t_star_ind_h = which(temp1 == max(temp1), arr.ind = T)[1,] 
      t1.h.star = t1h[t_star_ind_h[1]]  
      t2.h.star = t2h[t_star_ind_h[2]]
      ret2.hsic = (T1.hsic > t1.h.star)*(T2.hsic > t2.h.star)
      
    } else {ret2.hsic = rep(0,m) }
    
    
  } else {ret2.hsic = rep(0,m) }
  
  if(sum(G.c) > 0){
    
    
    temp2 = rej.c * G.c
    if(sum(temp2) > 0) {
      t_star_ind_c = which(temp2 == max(temp2), arr.ind = T)[1,] 
      t1.c.star = t1c[t_star_ind_c[1]]  
      t2.c.star = t2c[t_star_ind_c[2]]
      ret2.cor = (T1.cor > t1.c.star)*(T2.cor > t2.c.star)
      
    } else {ret2.cor = rep(0, m)}
    
    
  } else {ret2.cor = rep(0, m) }
  
  
  
  
  
  #return(cbind(ret0, ret1.hsic, ret1.cor, ret2.hsic, ret2.cor, ret5))
  return(c(sum(ret0), sum(ret1.hsic), sum(ret1.cor), sum(ret2.hsic), sum(ret2.cor), sum(ret5)))
  
  }
  
 ###Venn diagram at alpha = 0.05 
  new = rejection_fun(0.05)
  new_fac = ifelse(new == 0, FALSE, TRUE)
  data = as.data.frame(new_fac[,c(1,3,5,6)])
  value = 1:1201
  
  data = data.frame(value, data)
  colnames(data) = c("value", "2dFDR", "1dFDR-RV", "MF-2dFDR-RV","BH")
  ggvenn::ggvenn(data, c("BH", "1dFDR-RV", "MF-2dFDR-RV", "2dFDR"),fill_color = c("green4", "skyblue1","dodgerblue4", "gray42" ))
    #d = data.frame(ret0, ret1.hsic, ret1.cor, ret2.hsic, ret2.cor, ret5)
  #write.csv(d, "rejections_0.05_metabolomics.csv")
rejectionmat2 = matrix(nrow = 10, ncol = 6)
level_vec = seq(0.01, 0.21, by = 0.02)
for( i in 1:10)
{
  rejectionmat2[i,] = rejection_fun(level = level_vec[i])
  print(i)
}
write.csv(rejectionmat2, "metabolomics_output4.csv")
#rejection_fun(Y,X,Z, level = 0.05)



#####################################################################
###Source from here for the graph
#####################################################################
rejectionmat2 = read.csv("metabolomics_output4.csv")[,-1]
level_vec = seq(0.01, 0.21, by = 0.02)
rej = cbind(level_vec[1:10], rejectionmat2)
colnames(rej) = c("alpha", "2dFDR", "HSIC-1dFDR", "RV-1dFDR","HSIC-2dFDR+", "RV-2dFDR+",  "BH")

rej = as.data.frame(rej)
library(reshape2)
plot_data = melt(rej, id.var = "alpha")
colnames(plot_data) = c("alpha", "method", "value")
ggplot(plot_data, aes(x = alpha, y = value, group = method, colour = method)) +
  geom_point(aes( shape = method, color = method)) + 
  geom_line() +  theme_bw() + 
  scale_color_manual("method", labels = c("2dFDR", "HSIC-1dFDR", "RV-1dFDR","HSIC-2dFDR+", "RV-2dFDR+",  "BH"), 
                     values = c("gray42", "indianred1", "skyblue1", "indianred4", "dodgerblue4", "green4")) +
  scale_shape_manual("method", labels = c("2dFDR", "HSIC-1dFDR", "RV-1dFDR","HSIC-2dFDR+", "RV-2dFDR+",  "BH"),
                     values = c(7,8,16,9,18,15))+
  labs(x = "FDR", y = "Rejections")




