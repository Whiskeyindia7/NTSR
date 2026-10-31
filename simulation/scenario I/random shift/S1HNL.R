rm(list=ls())
setwd("/home/whiskeyindia7/revised_rs/uncentered")

## Packages
library(spatstat);library(mvtnorm);library(ggplot2);library(geoR);library(mgcv);library(np)
library(fields);library(classInt);library(nlme);library(convoSPAT);library(dcov)

source("iterative_rs.R")

## Constants
n.obs <- 100; n.iter <- 1000; alpha <- 0.05; radius <- 0.50
N.perm <- c(99, 299, 499, 799, 999)

n.K <- length(N.perm)
covariate.names <- c("C1", "C2", "C3", "C4")

array.dim <- c(n.iter, n.K, length(covariate.names))
array.names <- list(iteration=NULL, K=as.character(N.perm),covariate=covariate.names)

rej.test.lin <- array(NA, dim=array.dim, dimnames=array.names)
rej.test.nw  <- array(NA, dim=array.dim, dimnames=array.names)
rej.test.nl  <- array(NA, dim=array.dim, dimnames=array.names)

set.seed(1234)


for(a in 1:n.iter){
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
  Z4 <- eval.im(log(auxD))
  L <- eval.im(log(auxE))
  rm(auxA);rm(auxB); rm(auxC); rm(auxD)
  
  ## Observed covariates
  C1 <- Z1[test.points]; C2 <- Z2[test.points]; C3 <- Z3[test.points]
  C4 <- Z4[test.points]; W <- L[test.points]
  ## Now, construct observed data y
  loc<-cbind(test.points$x,test.points$y)
  
  nugget<-rnorm(n.obs)
  beta0<- -0.5; beta1<- 1; beta2 <- 1; 
  
  X<-cbind(C1,C2,C3,C4);
  X.fields<-list(C1=Z1, C2=Z2, C3=Z3, C4=Z4)
  y <- beta0 + 2*sin(2*C1) + 1.5*(C2^2-1) + 2*tanh(2*C4) + W + nugget
  X<-as.matrix(X)

  # Sample covariance
  reg.model<-"lin.GAM";test.method<-"var";test.statistic<-"cov.uncentered"
  test.lin <- iterative.test(X=X, X.fields=X.fields, y=y, reg.model=reg.model, test.method=test.method, radius=radius, N.perm=N.perm,
                             test.statistic=test.statistic, alpha=alpha, test.points=test.points, seed=a)
  reg.model<-"NW";test.method<-"var";test.statistic<-"cov.uncentered"
  test.nw <- iterative.test(X=X, X.fields=X.fields, y=y, reg.model=reg.model, test.method=test.method, radius=radius, N.perm=N.perm,
                            test.statistic=test.statistic, alpha=alpha, test.points=test.points, seed=a)
  reg.model<-"nl.GAM";test.method<-"var";test.statistic<-"cov.uncentered"
  test.nl <- iterative.test(X=X, X.fields=X.fields, y=y, reg.model=reg.model, test.method=test.method, radius=radius, N.perm=N.perm,
                            test.statistic=test.statistic, alpha=alpha, test.points=test.points, seed=a)
  

  ################################################
  ## Final selection results
  ################################################

  rej.test.lin[a,,] <- test.lin$selection
  rej.test.nw[a,,]  <- test.nw$selection
  rej.test.nl[a,,]  <- test.nl$selection


  if(a %% 50 == 0){
    cat(sprintf("%d of %d trials completed\n", a, n.iter))
  }
}

R <- list(
  rej.test.lin =
    apply(rej.test.lin, c(2,3), mean, na.rm=TRUE),

  rej.test.nw =
    apply(rej.test.nw, c(2,3), mean, na.rm=TRUE),

  rej.test.nl =
    apply(rej.test.nl, c(2,3), mean, na.rm=TRUE)
)

save(
  R,
  rej.test.lin,
  rej.test.nw,
  rej.test.nl,
  file="s1HNL_variousK.RData"
)