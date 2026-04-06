rm(list=ls())
setwd('/home/whiskeyindia7/theta/s2')
#setwd('/home/whiskeyindia7/commuting/y_data')
set.seed(1234)
## Packages
library(spatstat);library(mvtnorm);library(ggplot2);library(geoR);library(mgcv);library(np)
library(fields); library(classInt); library(nlme); library(convoSPAT); library(dcov)

source("iterative_test_rs.R")

# Constant 
n.obs<-100;n.iter<-1000
N.perm<-499;alpha<-0.05;radius<-0.50;theta<-0.75
rej.test.lin<-rej.test.nw<-rej.test.nl<-matrix(NA, nrow=n.iter, ncol=4)
rej.test.lin.tau<-rej.test.nw.tau<-rej.test.nl.tau<-matrix(NA, nrow=n.iter, ncol=4)
rej.test.lin.dcov<-rej.test.nw.dcov<-rej.test.nl.dcov<-matrix(NA, nrow=n.iter, ncol=4)
cov.lin<-cov.nw<-cov.nl<-tau.lin<-tau.nw<-tau.nl<-dcov.lin<-dcov.nw<-dcov.nl<-list()

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
  y <- beta0 + exp(C1)+ C2^2 + C4^3 + W + nugget
  X<-as.matrix(X)

  ## Examples with covariance + variance correction
  # Examples: linear GAM
  reg.model<-"lin.GAM";test.method<-"var";test.statistic<-"cov"
  test.lin <- iterative.test(X=X, X.fields=X.fields, y=y, reg.model=reg.model, test.method=test.method, radius=radius,N.perm=N.perm,
                             test.statistic=test.statistic, alpha=alpha, test.points=test.points,theta=theta)
  cov.lin[[a]]<-test.lin
  # Examples: Nadaraya-Watson
  reg.model<-"NW";test.method<-"var";test.statistic<-"cov"
  test.nw <- iterative.test(X=X, X.fields=X.fields, y=y, reg.model=reg.model, test.method=test.method, radius=radius,N.perm=N.perm,
                            test.statistic=test.statistic, alpha=alpha, test.points=test.points,theta=theta)
  cov.nw[[a]]<-test.nw
  # Examples: nonlinear GAM
  reg.model<-"nl.GAM";test.method<-"var";test.statistic<-"cov"
  test.nl <- iterative.test(X=X, X.fields=X.fields, y=y, reg.model=reg.model, test.method=test.method, radius=radius,N.perm=N.perm,
                            test.statistic=test.statistic, alpha=alpha, test.points=test.points,theta=theta)
  cov.nl[[a]]<-test.nl
  
  ## Examples with tau + variance correction
  # Examples: linear GAM
  reg.model<-"lin.GAM";test.method<-"var";test.statistic<-"tau"
  test.lin.tau <- iterative.test(X=X, X.fields=X.fields, y=y, reg.model=reg.model, test.method=test.method, radius=radius,N.perm=N.perm,
                                 test.statistic=test.statistic, alpha=alpha, test.points=test.points,theta=theta)
  tau.lin[[a]]<-test.lin.tau
  # Examples: Nadaraya-Watson
  reg.model<-"NW";test.method<-"var";test.statistic<-"tau"
  test.nw.tau <- iterative.test(X=X, X.fields=X.fields, y=y, reg.model=reg.model, test.method=test.method, radius=radius,N.perm=N.perm,
                                test.statistic=test.statistic, alpha=alpha, test.points=test.points,theta=theta)
  tau.nw[[a]]<-test.nw.tau
  # Examples: nonlinear GAM
  reg.model<-"nl.GAM";test.method<-"var";test.statistic<-"tau"
  test.nl.tau <- iterative.test(X=X, X.fields=X.fields, y=y, reg.model=reg.model, test.method=test.method, radius=radius,N.perm=N.perm,
                                test.statistic=test.statistic, alpha=alpha, test.points=test.points,theta=theta)
  tau.nl[[a]]<-test.nl.tau
  
  ## Examples with dcov + variance correction
  # Examples: linear GAM
  reg.model<-"lin.GAM";test.method<-"var";test.statistic<-"dcov"
  test.lin.dcov <- iterative.test(X=X, X.fields=X.fields, y=y, reg.model=reg.model, test.method=test.method, radius=radius,N.perm=N.perm,
                                  test.statistic=test.statistic, alpha=alpha, test.points=test.points,theta=theta)
  dcov.lin[[a]]<-test.lin.dcov
  
  # Examples: Nadaraya-Watson
  reg.model<-"NW";test.method<-"var";test.statistic<-"dcov"
  test.nw.dcov <- iterative.test(X=X, X.fields=X.fields, y=y, reg.model=reg.model, test.method=test.method, radius=radius,N.perm=N.perm,
                                 test.statistic=test.statistic, alpha=alpha, test.points=test.points,theta=theta)
  dcov.nw[[a]]<-test.nw.dcov
  
  # Examples: nonlinear GAM
  reg.model<-"nl.GAM";test.method<-"var";test.statistic<-"dcov"
  test.nl.dcov <- iterative.test(X=X, X.fields=X.fields, y=y, reg.model=reg.model, test.method=test.method, radius=radius,N.perm=N.perm,
                                 test.statistic=test.statistic, alpha=alpha, test.points=test.points,theta=theta)
  dcov.nl[[a]]<-test.nl.dcov
  
  
  rej.test.lin[a,]<-test.lin[[length(test.lin)]]
  rej.test.nw[a,]<-test.nw[[length(test.nw)]]
  rej.test.nl[a,]<-test.nl[[length(test.nl)]]
  
  rej.test.lin.tau[a,]<-test.lin.tau[[length(test.lin.tau)]]
  rej.test.nw.tau[a,]<-test.nw.tau[[length(test.nw.tau)]]
  rej.test.nl.tau[a,]<-test.nl.tau[[length(test.nl.tau)]]
  
  rej.test.lin.dcov[a,]<-test.lin.dcov[[length(test.lin.dcov)]]
  rej.test.nw.dcov[a,]<-test.nw.dcov[[length(test.nw.dcov)]]
  rej.test.nl.dcov[a,]<-test.nl.dcov[[length(test.nl.dcov)]]
}

rej.test.lin<-colMeans(rej.test.lin)
rej.test.nw<-colMeans(rej.test.nw)
rej.test.nl<-colMeans(rej.test.nl)
rej.test.lin.tau<-colMeans(rej.test.lin.tau)
rej.test.nw.tau<-colMeans(rej.test.nw.tau)
rej.test.nl.tau<-colMeans(rej.test.nl.tau)
rej.test.lin.dcov<-colMeans(rej.test.lin.dcov)
rej.test.nw.dcov<-colMeans(rej.test.nw.dcov)
rej.test.nl.dcov<-colMeans(rej.test.nl.dcov)

names(rej.test.lin)<-names(rej.test.nw)<-names(rej.test.nl)<-colnames(X)
names(rej.test.lin.tau)<-names(rej.test.nw.tau)<-names(rej.test.nl.tau)<-colnames(X)
names(rej.test.lin.dcov)<-names(rej.test.nw.dcov)<-names(rej.test.nl.dcov)<-colnames(X)

R<-list(rej.test.lin=rej.test.lin,
        rej.test.nw=rej.test.nw,
        rej.test.nl=rej.test.nl,
        rej.test.lin.tau=rej.test.lin.tau,
        rej.test.nw.tau=rej.test.nw.tau,
        rej.test.nl.tau=rej.test.nl.tau,
        rej.test.lin.dcov=rej.test.lin.dcov,
        rej.test.nw.dcov=rej.test.nw.dcov,
        rej.test.nl.dcov=rej.test.nl.dcov)


save(R, cov.lin, cov.nw, cov.nl,
     tau.lin, tau.nw, tau.nl,
     dcov.lin, dcov.nw, dcov.nl,
     file="s2_nl_theta075.RData")
