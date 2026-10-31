rm(list=ls())
## Packages
library(spatstat);library(mvtnorm);library(ggplot2);library(geoR);library(mgcv);library(np)
library(fields); library(classInt); library(nlme); library(convoSPAT);library(dcov)

## Functions
T2<-function(x,y){
  t<-mean(as.vector(abs(outer(x, x, "-"))))*mean(as.vector(abs(outer(y, y, "-"))))
  return(1/t)
}

make.nlgam.f <- function(X.reduce) {
  Xr <- as.data.frame(X.reduce)
  cn <- colnames(Xr)
  if (is.null(cn)) {
    cn <- paste0("V", seq_len(ncol(Xr)))
    colnames(Xr) <- cn
  }
  ucnt  <- sapply(Xr, function(z) length(unique(z)))
  terms <- sapply(seq_along(cn), function(j) {
    if (ucnt[j] <= 2) {
      sprintf("as.factor(%s)", cn[j])
    } else {
      k_val_internal <- min(10, ucnt[j] - 1)
      if(k_val_internal < 3) k_val_internal <- 3
      sprintf("s(%s, k=%d)", cn[j], k_val_internal)
    }
  })
  as.formula(paste("Y ~", paste(terms, collapse = " + "), "+ s(x, y)"))
}

iterative.test <- function(X, X.fields, y, test.points, reg.model,
                           test.method="var", test.statistic="cov", radius=0.5,
                           N.perm=c(99, 299, 499, 799, 999), alpha=0.05, seed=1234, verbose=FALSE){
  ##
  ## X = covariates, n x p matrix, colnames(X) must indicate the names of covariates
  ## X.fields = list of random fields of covariates
  ## y = response, n-dimensional vector
  ## test.points = spatstat object for location points
  ## reg.model = regression model. 1) lin.GAM, 2) NW or 3) nl.GAM can be used
  ## test.method = testing methods. 1) torus correction (tor), 2) variance correction (var)
  ## test.statistic = test statistic. 1) sample covariance (cov),
  ##                  2) Kendall's tau (tau), 3) distance covariance (dcov) 4
  ## radius = radius of uniform distribution for random shift vector
  ## N.perm = numbers of random shifts to be compared
  ## alpha = significance level. Default is 0.05
  ## seed = optional seed number
  ##

  N.perm.compare <- N.perm
  N.perm.max <- max(N.perm.compare)
  n.K <- length(N.perm.compare)

  n.obs <- length(y)
  Cnames <- colnames(X)
  loc <- cbind(test.points$x, test.points$y)

  if(!is.null(seed)){set.seed(seed)}

  ## Generate the maximum number of shifts only once.
  ## Every K uses the first K elements of this common shift sequence.
  shift.info <- vector("list", N.perm.max)
  shift.vectors <- matrix(NA, nrow=N.perm.max, ncol=2)
  colnames(shift.vectors) <- c("x", "y")

  if(test.method=="var"){
    for(k in 1:N.perm.max){
      jump <- runifdisc(1, radius=radius)
      shift.vectors[k,] <- c(jump$x, jump$y)

      test.points.shifted <- shift.ppp(test.points, c(jump$x, jump$y))
      W.reduced <- intersect.owin(test.points$window, test.points.shifted$window)
      test.points.reduced <- test.points[W.reduced]
      test.points.reduced.backshifted <- shift.ppp(test.points.reduced,
                                                   c(-jump$x, -jump$y))

      all.points <- cbind(test.points$x, test.points$y)
      reduced.points <- cbind(test.points.reduced$x, test.points.reduced$y)
      reduced.indices <- match(apply(reduced.points, 1, paste, collapse = ","),
                               apply(all.points, 1, paste, collapse = ","))

      shift.info[[k]] <- list(
        points.backshifted=test.points.reduced.backshifted,
        indices=reduced.indices,
        n=test.points.reduced$n
      )
    }
  } else if(test.method=="tor"){
    for(k in 1:N.perm.max){
      shift.info[[k]] <- rshift(test.points, edge="torus", radius=radius)
    }
  } else {
    stop("Unknown test.method")
  }

  iter.all <- vector("list", n.K)
  names(iter.all) <- as.character(N.perm.compare)

  rej.all <- matrix(0, nrow=n.K, ncol=length(Cnames))
  rownames(rej.all) <- as.character(N.perm.compare)
  colnames(rej.all) <- Cnames

  ## Each K has its own backward-selection path because different K values
  ## can remove different variables at a given step.
  for(l in 1:n.K){

    N.perm <- N.perm.compare[l]
    iter <- list()

    rej <- matrix(0, nrow=1, ncol=length(Cnames))
    colnames(rej) <- Cnames

    t <- 1
    selection.finished <- FALSE

    X.current <- as.matrix(X)
    X.fields.current <- X.fields
    names.current <- colnames(X.current)

    while(ncol(X.current) > 0 && !selection.finished){

      current_p <- ncol(X.current)
      pval_current <- matrix(Inf, nrow=1, ncol=current_p)
      colnames(pval_current) <- names.current

      for(i in 1:current_p){

        X.interest <- X.current[, i]
        Z.interest <- X.fields.current[[i]]

        if(current_p == 1){
          p <- 0
          res <- y
        } else {
          X.reduce <- as.matrix(X.current[, -i, drop=FALSE])
          colnames(X.reduce) <- names.current[-i]
          p <- ncol(X.reduce)
          data <- data.frame(Y=y, X.reduce, x=loc[,1], y=loc[,2])

          if(reg.model == "lin.GAM"){
            if(p == 1){
              lgammod <- gam(Y ~ X.reduce + s(x,y), data=data)
              res <- data$Y - cbind(1, X.reduce) %*% lgammod$coefficients[1:2]
            } else {
              terms <- paste0("X.reduce[,", seq_len(p), "]")
              f <- as.formula(paste("Y ~", paste(terms, collapse = " + "), "+ s(x,y)"))
              lgammod <- gam(f, data=data)
              res <- data$Y - cbind(1, X.reduce) %*% lgammod$coefficients[1:(p+1)]
            }
          } else if(reg.model == "NW"){
            bws <- npregbw(ydat=y, xdat=X.reduce, data=data,
                           regtype="lc", ckertype="epanechnikov")
            ksfit <- npreg(bws, residuals=T)
            res <- ksfit$resid
          } else if(reg.model == "nl.GAM"){
            f <- make.nlgam.f(X.reduce)
            nlgammod <- gam(f, data=data)
            nlgam.pred <- predict.gam(nlgammod, type="terms")
            res <- nlgammod$residuals + as.vector(nlgam.pred[, "s(x,y)"])
          }
        }

        if(test.method=="tor"){
          if(test.statistic %in% c("cov", "cov.uncentered")){ stat.obs <- cov(res, X.interest) }
          else if(test.statistic=="tau"){ stat.obs <- cor(res, X.interest, method="kendall") }
          else if(test.statistic=="dcov"){
            stat.obs <- dcov(res, X.interest, type="V") * n.obs * T2(res, X.interest)
          }

          simulated <- rep(NA, times=N.perm+1)
          simulated[1] <- stat.obs

          for(k in 1:N.perm){
            test.points.shift <- shift.info[[k]]
            cov.interest <- Z.interest[test.points.shift]

            if(test.statistic %in% c("cov", "cov.uncentered")){
              simulated[k+1] <- cov(res, cov.interest)
            } else if(test.statistic=="tau"){
              simulated[k+1] <- cor(cov.interest, res, method="kendall")
            } else if(test.statistic=="dcov"){
              simulated[k+1] <- dcov(res, cov.interest, type="V") * n.obs * T2(res, cov.interest)
            }
          }
        }
        else if(test.method=="var"){
          if(test.statistic %in% c("cov", "cov.uncentered")){ stat.obs <- cov(res, X.interest) }
          else if(test.statistic=="tau"){ stat.obs <- cor(res, X.interest, method="kendall") }
          else if(test.statistic=="dcov"){
            stat.obs <- dcov(res, X.interest, type="V") * T2(res, X.interest)
          }

          simulated <- rep(NA, times=N.perm+1)
          simulated[1] <- stat.obs
          n.simulated <- rep(NA, times=N.perm+1)
          n.simulated[1] <- n.obs

          for(k in 1:N.perm){
            test.points.reduced.backshifted <- shift.info[[k]]$points.backshifted
            reduced.indices <- shift.info[[k]]$indices

            res.reduced <- res[reduced.indices]
            vals2 <- Z.interest[test.points.reduced.backshifted]

            if(test.statistic %in% c("cov", "cov.uncentered")){
              simulated[k+1] <- cov(res.reduced, vals2)
            } else if(test.statistic=="tau"){
              simulated[k+1] <- cor(res.reduced, vals2, method="kendall")
            } else if(test.statistic=="dcov"){
              simulated[k+1] <- dcov(res.reduced, vals2, type="V") * T2(res.reduced, vals2)
            }
            n.simulated[k+1] <- shift.info[[k]]$n
          }

          if(test.statistic=="dcov"){
            simulated <- simulated * n.simulated
            stat.obs <- simulated[1]
          } else if(test.statistic=="cov.uncentered"){
            simulated <- simulated * sqrt(n.simulated)
          } else {
            simulated <- (simulated - mean(simulated)) * sqrt(n.simulated)
          }
        }

        if(test.statistic == "dcov"){
          pval_current[, i] <- sum(simulated >= stat.obs) / (N.perm + 1)
        } else {
          test.rank <- rank(simulated)[1]
          pval_current[, i] <- min(1, 2 * min(test.rank, N.perm + 2 - test.rank) / (N.perm + 1))
        }
      }

      iter[[t]] <- list(pval = pval_current)

      if(max(pval_current) >= alpha){

        if(current_p == 1){
          if(verbose){cat(sprintf("K=%d: No covariates are significant\n", N.perm))}
          iter[[t+1]] <- rej
          selection.finished <- TRUE
        } else {
          M <- which.max(pval_current)
          
          if(verbose){
            cat(sprintf(
              "K=%d, Step %d: Variable '%s' is removed (p-value: %.4f)\n",
              N.perm, t, names.current[M], pval_current[M]))
          }
          X.current <- as.matrix(X.current[, -M, drop=FALSE])
          X.fields.current <- X.fields.current[-M]
          names.current <- colnames(X.current)
          colnames(X.current) <- names.current
          t <- t + 1
        }
      } else {
        m <- match(names.current, colnames(rej))
        rej[, m] <- 1
        iter[[t+1]] <- rej
        if(verbose){
          cat(sprintf(
            "K=%d: Backward Selection Completed. Significant variables found.\n",
            N.perm))
        }
        selection.finished <- TRUE
      }
    }

    iter.all[[l]] <- iter
    rej.all[l,] <- rej
  }

  if(n.K == 1){
    return(iter.all[[1]])
  } else {
    return(list(
      N.perm=N.perm.compare,
      results=iter.all,
      selection=rej.all
    ))
  }
}
