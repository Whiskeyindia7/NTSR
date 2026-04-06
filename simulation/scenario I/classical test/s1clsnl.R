rm(list=ls())
setwd('/Users/kimhyoeun/Desktop/code revised/data/classical')
#setwd('/Users/whiskeyindia/Library/Mobile Documents/com~apple~CloudDocs/spatialregression/realdata/georgia_census')
#setwd('/home/whiskeyindia7/commuting/logy_data')
set.seed(123)
## Packages
library(spatstat);library(mvtnorm);library(ggplot2);library(geoR);library(mgcv);library(np)
library(fields); library(classInt); library(nlme); library(convoSPAT)

robust.gls<-function(formula, data, correlation, method="ML"){
  tryCatch(
    gls(formula, data=data, correlation=correlation, method=method, control=glsControl(opt="optim")),
    error = function(e){
      warning("GLS convergence fail: ", e$message, "\n Replaced by gls with no tolerance")
      gls(formula, data=data, correlation=correlation, method=method, 
          control=glsControl(opt="optim", maxIter=100, msMaxIter=100, msTol=1, tolerance=1))
    }
  )
}

iterative.test <- function(X, y, test.points, alpha = 0.05) {
  ##
  ## X : n x p matrix 
  ## y : response vector
  ## test.points : spatstat point pattern
  ## alpha : significance level
  ##
  
  # setup
  cov.names <- colnames(X)
  loc <- cbind(test.points$x, test.points$y)
  data <- data.frame(Y = y,X,x = loc[, 1],yloc = loc[, 2])
  iter <- list()
  t <- 1
  
  # rejection indicator
  rej <- matrix(0, nrow = 1, ncol = length(cov.names))
  colnames(rej) <- cov.names
  repeat {
    p <- length(cov.names)
    f <- as.formula(paste("Y ~", paste(cov.names, collapse = " + ")))
    mod <- robust.gls(f,data = data,correlation = corGaus(c(0.2, 0.1), ~ x + yloc, nugget = TRUE),method = "ML")
    sum.mod <- summary(mod)
    pval <- sum.mod$tTable[cov.names, "p-value"]
    pval <- matrix(pval, nrow = 1)
    colnames(pval) <- cov.names
    iter[[t]] <- pval
    if(max(pval) >= alpha){
      if(p == 1){
        message("No covariates are significant.")
        iter[[t+1]] <- rej
        return(iter)
      }
      drop.var <- (which.max(pval))
      cov.names <- cov.names[-drop.var]
      t<-t+1
    }else{
      rej[, cov.names] <- 1
      iter[[t+1]] <- rej
      return(iter)
    }
  }
}

n.iter<-1000;n.obs<-100;alpha<-0.05;R<-matrix(NA, nrow=n.iter, ncol=4)

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
  y <- beta0 + exp(C1) + C2^2 + C4^3 + W + nugget
  X<-as.matrix(X)
  
  res<-iterative.test(X=X,y=y,test.points=test.points)
  R[a,]<-res[[length(res)]]
}

R<-colMeans(R)
names(R)<-colnames(X)

save(R, file="s1_nl_cls.RData")