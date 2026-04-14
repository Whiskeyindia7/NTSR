##########################
#### classical method ####
##########################

load("Chlorophyll.RData")
set.seed(123)
## Packages
library(spatstat);library(mvtnorm);library(ggplot2);library(geoR);library(mgcv);library(np)
library(fields); library(classInt); library(nlme); library(convoSPAT)


### backward ###
iterative.test <- function(X, y, test.points, alpha = 0.05, nu = 2.5){
  
  cov.names <- colnames(X)
  loc <- cbind(test.points$x, test.points$y)
  
  data <- data.frame(
    Y = y,
    X,
    x = loc[,1],
    yloc = loc[,2]
  )
  
  iter <- list()
  t <- 1
  rej <- matrix(0, nrow = 1, ncol = length(cov.names))
  colnames(rej) <- cov.names
  repeat {
    p <- length(cov.names)
    f <- as.formula(
      paste("Y ~", paste(cov.names, collapse = " + "),
            "+ Matern(1 | x + yloc)"
      )
    )
    mod <- spaMM::fitme(
      f,
      data = data,
      family = gaussian(),
      method = "ML",
      fixed = list(nu = nu)
    )
    sm <- capture.output(
      tab <- as.data.frame(summary(mod, details = list(p_value = TRUE))$beta_table)
    )
    pval <- tab[cov.names, "p-value"]
    pval <- matrix(pval, nrow = 1)
    colnames(pval) <- cov.names
    iter[[t]] <- pval
    if(max(pval) >= alpha){
      if(p == 1){
        message("No covariates are significant.")
        iter[[t+1]] <- rej
        return(iter)
      }
      drop.var <- which.max(pval)
      cov.names <- cov.names[-drop.var]
      t <- t + 1
    } else {
      rej[, cov.names] <- 1
      iter[[t+1]] <- rej
      return(iter)
    }
  }
}

backward <- iterative.test(X=X, y=y, test.points=test.points, alpha=0.05, nu=2.5)
backward
