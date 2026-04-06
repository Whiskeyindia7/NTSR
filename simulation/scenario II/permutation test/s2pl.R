rm(list = ls())
setwd('/home/whiskeyindia7/theta/s2')
set.seed(123)
library(spatstat);library(fields);library(mvtnorm);library(ggplot2);library(geoR);library(mgcv);library(np);library(expm)
library(doParallel);library(foreach);library(convoSPAT);library(sp);library(nlme)
## Constants
num<-1; n.obs<-100; n.iter<-1000; N.perm<-499; 

source('iterative_test_gls.R')
rej<-matrix(NA, nrow=n.iter, ncol=4)

for(a in 1:n.iter){
  # 1. Data generation
  window <- owin(c(0,1),c(0,1))
  test.points <- runifpoint(n.obs, win=window)   ## Observed locations
  auxA <- attr(rLGCP("exp", mu=0, var=1, scale=0.2, saveLambda=TRUE, win=window),"Lambda") # Covariate 
  auxB <- attr(rLGCP("stable", mu=0, var=1, scale=0.1, alpha=0.5, saveLambda=TRUE, win=window),"Lambda") # Covariate 
  auxC <- attr(rLGCP("matern", mu=0, var=1, scale=0.1, saveLambda=TRUE, win=window, nu=2),"Lambda") # Random effect 
  auxD <- attr(rLGCP("gencauchy", mu=0, var=1, scale=0.1, saveLambda=TRUE, win=window, alpha=2, beta=5),"Lambda") # Random effect 
  auxE <- attr(rLGCP("gau", mu=0, var=1, scale=0.2, saveLambda=TRUE, win=window),"Lambda") # Random effect 
  
  Z1 <- eval.im(log(auxA))
  Z2 <- eval.im(log(auxB))
  Z3 <- eval.im(log(auxC))
  L <- eval.im(log(auxE))
  
  auxD <- attr(rLGCP("exp", mu=0.8*Z1, var=1, scale=0.2, saveLambda=TRUE, win=window),"Lambda") # Random effect 
  Z4 <- eval.im(log(auxD))
  rm(auxA);rm(auxB); rm(auxC); rm(auxD); rm(auxE)
  
  ## Observed covariates
  C1 <- Z1[test.points]; C2 <- Z2[test.points]; C3 <- Z3[test.points]
  C4 <- Z4[test.points]; 
  
  ## Now, construct observed data y
  loc<-cbind(test.points$x,test.points$y)
  W<-L[test.points]
  nugget<-rnorm(n.obs)
  beta0<- -0.5; beta1<- 1; beta2 <- 1; 
  
  X<-cbind(C1,C2,C3,C4);
  X.fields<-list(C1=Z1, C2=Z2, C3=Z3, C4=Z4)
  y <- beta0 + C1 + C2 + C4 + W + nugget
  X<-as.matrix(X)
  test<-iterative.test(X=X, y=y, test.points=test.points)
  rej[a,]<-test[[length(test)]]
  
  if(a%%100==0){
    print(paste0("Iteration ", a, " out of ", n.iter))
    rej.a<-colMeans(rej[1:a,])
    names(rej.a)<-colnames(X)
    filename<-paste0("s2pl2", a, ".RData")
    save(rej.a, file=filename)
  }
}

rej<-colMeans(rej)
names(rej)<-colnames(X)
save(rej, file="s2pl2.RData")
