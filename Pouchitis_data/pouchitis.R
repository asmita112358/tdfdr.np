####Analysis of pouchitis data
###loading libraries
library(nnet)
library(sgof)
##Reading data
tab <- read.table('gene_expression_meta_file.txt', sep = '\t', head = TRUE, row.names = 1)
geneexpression = read.table('gene_expression.tab', sep = '\t', head = TRUE, row.names = 1)

##Location: PPI
index_ppi = tab$Location == 'PPI'
tab1 <- tab[index_ppi, ]
Y1 = geneexpression[,index_ppi]
dim(Y1)
m = nrow(Y1)
n = ncol(Y1)

X = as.factor(tab1$Outcome)
Z1 = as.factor(tab1$Antibiotics)
Z2 = as.factor(tab1$Gender)
Z3 = as.factor(tab1$Batch)

##Generating Additional Xb from Z1, Z2, Z3.
B = 500
propensity = multinom(X ~ Z1 + Z2 + Z3)$fitted.values
Xb = matrix(nrow = B, ncol = n)
xvec = as.character(levels(X))
for(b in 1:B)
{
  for(i in 1:n)
  {
    Xb[b, i] = xvec[as.logical(rmultinom(1, 1, prob = propensity[i,]))]
  }
 
}

###Computing the statistics.
T1b <- T2b <- matrix(NA, nrow = m, ncol = B)
T1 <- T2 <- vector()
for(j in 1:m)
{
  y = as.vector(unlist(Y1[j,]))
   
  T1[j] = summary(aov(y ~ X + Z1 + Z2 + Z3))[[1]][1,4]
  T2[j] = summary(aov(y~X))[[1]][1,4]
  for(b in 1:B)
  {
    
    T1b[j, b] = summary(aov(y ~ as.factor(Xb[b,]) + Z1 + Z2 + Z3))[[1]][1,4]
    T2b[j, b] = summary(aov(y ~ as.factor(Xb[b,])))[[1]][1,4]
  }
}

####Calculating the null proportion

lambda = mean(T1)

#pi0= min(1, mean((T1 <= lambda)/ if_else(rep(sum(rowMeans(T1b <= lambda)), m) > rep(0, m), rowMeans(T1b <= lambda), rep(1,m))))

pi0 = min(1, mean(T1 <= lambda )/ mean(T1b <= lambda))

########################################################################################
##Grid search with the two statistics.


ngrid = 400


rejection_fun_pouchitis = function(level)
{
  p = vector()
  for(i in 1:m)
  {
    
    y = as.vector(unlist(Y1[i,]))
    
    p[i] = summary(aov(y ~ X + Z1 + Z2 + Z3))[[1]][1,5]
  }
  benj = BH(p, alpha = level)
  ret0 = benj$Rejections
  print("BH done")
  ##2DFDR
  t1t = quantile(T1, probs = seq(0, 1, length.out = ngrid))
  t2t = quantile(as.vector(unlist(T2)), probs = seq(0, 1, length.out = ngrid))
  
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
   print(i) 
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
  
  write.csv(ret, "output_pouchitis_400.csv")
  print("MF-2DFDR done")
  ##1DFDR
  rej.1d <- fdp.1d <- rep(NA,m)
  for(k in 1:m)
  {
    Num.1d<- sum(T1b > T1[k])/(B+1)
    Dem.1d <- sum(T1 > T1[k])
    fdp.1d[k] <- Num.1d/max(1,Dem.1d)
    rej.1d[k] <- Dem.1d
  }
  
  
  if(sum(fdp.1d < level)>0) 
  {
    ind <- which(rej.1d == max(rej.1d[fdp.1d < level]))
    ind.s <- ind[which.min(fdp.1d[ind])]
    ret.1d <- (T1 > T1[ind.s])
    
  } else {ret.1d = rep(0, m)} 
  print("1DFDR done")
  
  return(c(ret0, sum(ret), sum(ret.1d)))
  
  
  
}

u = rejection_fun_pouchitis(0.05)
write.csv(u, "output_pouchitis.csv")
##Location: Pouch
#ndex_pouch = tab$Location == 'Pouch'
#tab2 = tab[index_pouch,]
#Y2 = geneexpression[,index_pouch]
#dim(Y2)


##Analysis


###Data from PPI did not work. However data from Pouch worked very well
###Here is the code for the graph from the pouch data.

data = read.csv("pouchitis_pouch_all.csv")
data = data[,-1]
colnames(data) <- c("alpha", "BH", "MS-2dFDR+", "MS-1dFDR")
data = as.data.frame(data)
library(reshape2)
plot_data = melt(data, id.var = "alpha")
colnames(plot_data) = c("alpha", "method", "value")
library(ggplot2)
ggplot(plot_data, aes(x = alpha, y = value, group = method, colour = method)) +
  geom_point(aes( shape = method, color = method)) + 
  geom_line() + theme_bw() +
  scale_color_manual("method",label = c("BH",  "MS-2dFDR+","MS-1dFDR"), values = c("green4","goldenrod4", "goldenrod") )+ 
  scale_shape_manual("method", label = c("BH",  "MS-2dFDR+", "MS-1dFDR"), values = c(15,3, 17))+
  
  labs(x = "FDR", y = "Rejections")
