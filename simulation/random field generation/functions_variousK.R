
robust.gls<-function(formula, data, correlation, method="ML"){
  tryCatch(
    gls(formula, data=data, correlation=correlation, method=method),
    error = function(e){
      warning("GLS convergence fail: ", e$message, "\n Replaced by gls with no tolerance")
      gls(formula, data=data, correlation=correlation, method=method, 
          control=glsControl(opt="optim", maxIter=100, msMaxIter=100, msTol=1, tolerance=1))
    }
  )
}

checkerboard <- function(coords, values, n, C) {
  if (!is.matrix(coords) || ncol(coords) != 2) {
    stop("`coords` must be a K×2 numeric matrix")
  }
  K <- nrow(coords)
  if (length(values) != K) {
    stop("Length of `values` must match number of rows in `coords`")
  }
  if (!is.numeric(n) || length(n) != 1 || n <= 0) {
    stop("`n` must be a positive integer")
  }
  if (!is.numeric(C) || length(C) != 1) {
    stop("`C` must be a single numeric value")
  }
  
  x <- coords[,1]
  y <- coords[,2]
  
  # Determine grid cell indices (1..n) for each point
  i <- pmin(floor(x * n) + 1, n)
  j <- pmin(floor(y * n) + 1, n)
  
  # Checkerboard sign: +1 for even-sum cells, -1 for odd-sum
  sign.pattern <- ifelse(((i + j) %% 2) == 0,  1, -1)
  
  # Apply adjustment
  adjusted.values <- values + sign.pattern * C
  return(adjusted.values)
}


T2<-function(x,y){
  t<-mean(as.vector(abs(outer(x, x, "-"))))*mean(as.vector(abs(outer(y, y, "-"))))
  return(1/t)
}

simulation.test.variousK<-function(n.trial=2000, N.perm=c(99, 299, 499, 799, 999), n.obs=100, test.statistic="cov", true.cor="L",
                                   window.bound=1, window=NULL, test.target="sig", scenario="SE1", 
                                   test.method="var", radius=window.bound/2, alpha=0.05, seed=123,
                                   max.shift.attempts=10000){
  #
  # INPUT DESCRIPTION
  # n.trial ... the number of repetition for whole test
  # N.perm ... the number of permutation (random shift)
  # n.obs ... the number of observation in window
  # test.statistic ... the test statistic. "cov", "tau" and "dcov" can be used
  # true.cor ... true correlation structure between y and C. "L" is linear, "S" is squared and "E" means exponential 
  # window.bound ... the length of side of window. default is 1 (unit square). useless if window is specified.
  # test.target ... power ("pow") or signinficance level ("sig"). default is significance test 
  # window ... specified window. if you want to specify the observational window.
  # scenario ... scenario. SE1, NS, LN, N, E1, SE4 can be used
  # test.method ... variance correction ("var") and torrus correction ("tor")
  # radius ... radius for shift vector
  # seed ... seed number for set.seed()
  # max.shift.attempts ... maximum number of attempts used to obtain a shift
  #                        with at least two retained observations
  #
  N.perm.compare <- N.perm
  N.perm <- max(N.perm.compare)
  n.K <- length(N.perm.compare)

  if(!is.numeric(max.shift.attempts) || length(max.shift.attempts) != 1L ||
     !is.finite(max.shift.attempts) || max.shift.attempts < 1){
    stop("`max.shift.attempts` must be a positive finite number.")
  }
  max.shift.attempts <- as.integer(max.shift.attempts)
  
  if(test.statistic=="dcov"){
    lgam.rej<-ks.rej<-nlgam.rej<-rep(0,n.K)
    lgam.p<-ks.p<-nlgam.p<-matrix(NA,nrow=n.trial,ncol=n.K)
    
    colnames(lgam.p)<-N.perm.compare
    colnames(ks.p)<-N.perm.compare
    colnames(nlgam.p)<-N.perm.compare
    
    names(lgam.rej)<-N.perm.compare
    names(ks.rej)<-N.perm.compare
    names(nlgam.rej)<-N.perm.compare
  }else{
    lgam.rej.uncentered<-ks.rej.uncentered<-nlgam.rej.uncentered<-rep(0,n.K)
    lgam.p.uncentered<-ks.p.uncentered<-nlgam.p.uncentered<-matrix(NA,nrow=n.trial,ncol=n.K)

    colnames(lgam.p.uncentered)<-N.perm.compare
    colnames(ks.p.uncentered)<-N.perm.compare
    colnames(nlgam.p.uncentered)<-N.perm.compare

    names(lgam.rej.uncentered)<-N.perm.compare
    names(ks.rej.uncentered)<-N.perm.compare
    names(nlgam.rej.uncentered)<-N.perm.compare
  }
  
  simulated.lgam.list <- vector("list", n.trial)
  simulated.ks.list <- vector("list", n.trial)
  simulated.nlgam.list <- vector("list", n.trial)
  n.simulated.list <- vector("list", n.trial)
  
  set.seed(seed)
  
  for(i in 1:n.trial){
    #
    # DATA DESCRIPTION
    # test.points ... the set of sampling locations, here represented as a point pattern
    #                 (the default observation window is the unit square)
    # Z1 .... (nuisance) one of the random fields, available in the whole observation window,
    #                 here represented as a pixel image
    # Z2 .... (interest) one of the random fields, available in the whole observation window,
    #                 here represented as a pixel image
    # y .... observational data, available only at test.points
    #
    
    if(is.null(window)){
      window <- owin(c(0,window.bound),c(0,window.bound))
    }else{window<-window}
    
    test.points <- runifpoint(n.obs, win=window)   ## Observed locations
    
    
    auxA <- attr(rLGCP("exp", mu=0, var=1, scale=0.2, saveLambda=TRUE, win=window),"Lambda") # Covariate 
    auxB <- attr(rLGCP("exp", mu=0, var=1, scale=0.2, saveLambda=TRUE, win=window),"Lambda") # Covariate 
    Z1 <- eval.im(log(auxA))
    Z2 <- eval.im(log(auxB))
    rm(auxA);rm(auxB);
    ## Observed covariates
    C1 <- Z1[test.points]; C2 <- Z2[test.points];
    
    ## Now, construct observed data y
    loc<-cbind(test.points$x,test.points$y)
    d<-rdist(loc)
    nugget<-rnorm(n.obs)
    beta0<- -0.5;
    
    if(scenario=="SE1"){
      auxC <- attr(rLGCP("gau", mu=0, var=1, scale=0.2, saveLambda=TRUE, win=window),"Lambda") # Random effect 
      Z3 <- eval.im(log(auxC))
      W <- Z3[test.points];rm(auxC)
    }else if(scenario=="NS"){
      conv.kernels<-f_mc_kernels(y.min=0, y.max=window.bound, x.min=0, x.max=window.bound, N.mc=4)
      conv.w<-NSconvo_sim(grid=F, y.min=0, y.max=window.bound, x.min=0, x.max=window.bound, N.obs=n.obs, 
                          sim.locations=loc, sigmasq=1, tausq=0, mc.kernels.obj=conv.kernels,kappa=.5,
                          beta.coefs=0, cov.model="matern")
      W <- conv.w$sim.data
    }else if(scenario=="LN"){
      auxC <- attr(rLGCP("gau", mu=0, var=1, scale=0.2, saveLambda=TRUE, win=window),"Lambda") # Random effect 
      Z3 <- eval.im(auxC)
      W <- (Z3[test.points]);rm(auxC)
    }
    else if(scenario=="N"){
      auxC <- attr(rLGCP("gau", mu=0, var=1, scale=0.2, saveLambda=TRUE, win=window),"Lambda") # Random effect 
      Z3 <- eval.im(log(auxC))
      W <- Z3[test.points];rm(auxC)
      W<-checkerboard(coords=loc, values=W, n=4, C=0.5)
    }
    else if(scenario=="E1"){
      auxC <- attr(rLGCP("exp", mu=0, var=1, scale=0.2, saveLambda=TRUE, win=window),"Lambda") # Random effect 
      Z3 <- eval.im(log(auxC))
      W <- Z3[test.points];rm(auxC)
    }
    else if(scenario=="SD"){
      Sigma.SD <- exp(-(d / 0.4)^2)
      W <- as.numeric(mvtnorm::rmvnorm(n=1, mean=rep(0, n.obs), sigma=Sigma.SD))
    }
    else if(scenario=="HT"){
      df<-3
      Sigma.SE1 <- exp(-(d / 0.2)^2)
      Sigma.t <- (df - 2) / df * Sigma.SE1
      W <- as.numeric(mvtnorm::rmvt(n=1, sigma=Sigma.t, df=df, delta=rep(0,n.obs)))
    }
    else if(scenario=="SE4"){
      auxC <- attr(rLGCP("gau", mu=0, var=4, scale=0.2, saveLambda=TRUE, win=window),"Lambda") # Random effect 
      Z3 <- eval.im(log(auxC))
      W <- Z3[test.points];rm(auxC)
    }else{stop("Unknown value for scenario :", scenario)}
    
    if(test.target=="sig"){
      if(true.cor=="L"){y <- beta0 +  C1 + W + nugget}
      else if(true.cor=="S"){y <- beta0 +  C1^2 + W + nugget}
      else{stop("Unknown value for true.cor :", true.cor)}
    }
    else if(test.target=="pow"){
      if(true.cor=="L"){y <- beta0 + C1 + C2 + W + nugget}
      else if(true.cor=="S"){y <- beta0 + C1^2 + C2^2 + W + nugget}
      else{stop("Unknown value for true.cor :", true.cor)}
    }
    
    data <- data.frame("response" = y, "covariate1" = C1, "covariate2" = C2, 
                       "x" = test.points$x, "y" = test.points$y)
    
    ############################################
    ## Now, regression moedel fitting: y ~ C1 ##
    ############################################
    
    ##################
    ## Linear Model ##
    ##################
    ## Method: Linear GAM (y = XB + f(s) + epsilon)
    lgammod <- gam(response ~ covariate1 + s(x,y), data=data)
    lgamres<- data$response - (lgammod$coefficients[1] + lgammod$coefficients[2]*data$covariate1)
    
    #####################
    ## Nonlinear model ##
    #####################
    
    ## Method: NW-estimation for mean trend
    bws<-npregbw(response ~ covariate1, data=data, regtype="lc", ckertype="epanechnikov", verbose = F)
    ksfit<-npreg(bws, residuals=T)
    ksres<-ksfit$resid
    
    ## Method: Nonlinear GAM
    nlgammod<-gam(response ~ s(covariate1) + s(x,y), data=data)
    nlgam.pred<-predict.gam(nlgammod, type="terms")
    nlgamres<-nlgammod$residuals + as.vector(nlgam.pred[,"s(x,y)"])
    
    if(test.method=="tor"){
      
      # test statistic
      if(test.statistic=="cov"){
        stat.lgam<-cov(lgamres, C2)
        stat.ks<-cov(ksres, C2)
        stat.nlgam<-cov(nlgamres, C2)
      }else if(test.statistic=="tau"){
        stat.lgam<-cor(lgamres, C2, method="kendall")
        stat.ks<-cor(ksres, C2, method="kendall")
        stat.nlgam<-cor(nlgamres, C2, method="kendall")
      }else if(test.statistic=="dcov"){
        stat.lgam<-dcov(lgamres, C2, type="V")*n.obs*T2(lgamres, C2)
        stat.ks<-dcov(ksres, C2, type="V")*n.obs*T2(ksres, C2)
        stat.nlgam<-dcov(nlgamres, C2, type="V")*n.obs*T2(nlgamres, C2)
      }
      
      
      simulated.lgam<-rep(NA, times=N.perm+1)
      simulated.ks<-rep(NA, times=N.perm+1)
      simulated.nlgam<-rep(NA, times=N.perm+1)
      
      simulated.lgam[1]<-stat.lgam
      simulated.ks[1]<-stat.ks
      simulated.nlgam[1]<-stat.nlgam
      
      
      for (k in 1:N.perm){
        test.points.shift <- rshift(test.points, edge="torus", radius=radius)
        cov.interest<-Z2[test.points.shift]
        if(test.statistic=="cov"){
          simulated.lgam[k+1] <- cov(cov.interest,lgamres)
          simulated.ks[k+1] <- cov(cov.interest,ksres)
          simulated.nlgam[k+1] <- cov(cov.interest,nlgamres)
        }else if(test.statistic=="tau"){
          simulated.lgam[k+1] <- cor(cov.interest, lgamres,method="kendall")
          simulated.ks[k+1] <- cor(cov.interest, ksres,method="kendall")
          simulated.nlgam[k+1] <- cor(cov.interest, nlgamres,method="kendall")
        }else if(test.statistic=="dcov"){
          simulated.lgam[k+1]<-dcov(lgamres, cov.interest, type="V")*n.obs*T2(lgamres, cov.interest)
          simulated.ks[k+1]<-dcov(ksres, cov.interest, type="V")*n.obs*T2(ksres, cov.interest)
          simulated.nlgam[k+1]<-dcov(nlgamres, cov.interest, type="V")*n.obs*T2(nlgamres, cov.interest)
        }
      }
      
    }else if(test.method=="var"){
      
      if(test.statistic=="cov"){
        stat.lgam<-cov(lgamres, C2)
        stat.ks<-cov(ksres, C2)
        stat.nlgam<-cov(nlgamres, C2)
      }else if(test.statistic=="tau"){
        stat.lgam<-cor(lgamres, C2, method="kendall")
        stat.ks<-cor(ksres, C2, method="kendall")
        stat.nlgam<-cor(nlgamres, C2, method="kendall")
      }else if(test.statistic=="dcov"){
        stat.lgam<-dcov(lgamres, C2, type="V")*T2(lgamres, C2)
        stat.ks<-dcov(ksres, C2, type="V")*T2(ksres, C2)
        stat.nlgam<-dcov(nlgamres, C2, type="V")*T2(nlgamres, C2)
      }
      
      
      simulated.lgam<-rep(NA, times=N.perm+1)
      simulated.ks<-rep(NA, times=N.perm+1)
      simulated.nlgam<-rep(NA, times=N.perm+1)
      
      simulated.lgam[1]<-stat.lgam
      simulated.ks[1]<-stat.ks
      simulated.nlgam[1]<-stat.nlgam
      n.simulated <- rep(NA, times=N.perm+1)
      n.simulated[1] <- test.points$n
      
      # Random shifts
      for (k in 1:N.perm){
        shift.attempt <- 0L

        repeat {
          shift.attempt <- shift.attempt + 1L
          jump <- runifdisc(1, radius=radius)
          test.points.shifted <- shift(test.points, c(jump$x,jump$y))
          W.reduced <- intersect.owin(test.points$window, test.points.shifted$window)
          test.points.reduced <- test.points[W.reduced]

          # cov(), Kendall's tau, and the normalized dCov statistic require
          # at least two retained observations. Invalid shifts are discarded.
          if(test.points.reduced$n >= 2L){
            break
          }

          if(shift.attempt >= max.shift.attempts){
            stop(sprintf(
              "Failed to obtain at least two retained observations for shift %d after %d attempts.",
              k, max.shift.attempts
            ))
          }
        }

        test.points.reduced.backshifted <- shift(test.points.reduced, c(-jump$x,-jump$y))
        
        all.points <- cbind(test.points$x, test.points$y)
        reduced.points <- cbind(test.points.reduced$x, test.points.reduced$y)
        reduced.indices <- match(apply(reduced.points, 1, paste, collapse = ","), apply(all.points, 1, paste, collapse = ","))
        
        lgamres.reduced<-lgamres[reduced.indices]
        ksres.reduced<-ksres[reduced.indices]
        nlgamres.reduced<-nlgamres[reduced.indices]
        
        cov.interest<-Z2[test.points.reduced.backshifted]
        
        if(test.statistic=="cov"){
          simulated.lgam[k+1] <- cov(lgamres.reduced, cov.interest)
          simulated.ks[k+1] <- cov(ksres.reduced, cov.interest)
          simulated.nlgam[k+1] <- cov(nlgamres.reduced, cov.interest)
          n.simulated[k+1] <- test.points.reduced$n
        }else if(test.statistic=="tau"){
          simulated.lgam[k+1] <- cor(lgamres.reduced, cov.interest, method="kendall")
          simulated.ks[k+1] <- cor(ksres.reduced, cov.interest, method="kendall")
          simulated.nlgam[k+1] <- cor(nlgamres.reduced, cov.interest, method="kendall")
          n.simulated[k+1] <- test.points.reduced$n
        }else if(test.statistic=="dcov"){
          simulated.lgam[k+1] <- dcov(lgamres.reduced, cov.interest, type="V")*T2(lgamres.reduced, cov.interest)
          simulated.ks[k+1] <- dcov(ksres.reduced, cov.interest, type="V")*T2(ksres.reduced, cov.interest)
          simulated.nlgam[k+1] <- dcov(nlgamres.reduced, cov.interest, type="V")*T2(nlgamres.reduced, cov.interest)
          n.simulated[k+1] <- test.points.reduced$n
        }
      }
      
      # Save raw test statistics
      simulated.lgam.list[[i]] <- simulated.lgam
      simulated.ks.list[[i]] <- simulated.ks
      simulated.nlgam.list[[i]] <- simulated.nlgam
      n.simulated.list[[i]] <- n.simulated
      
      # Variance correction for different numbers of random shifts
      for(j in 1:n.K){
        
        current.index <- 1:(N.perm.compare[j]+1)
        current.M <- N.perm.compare[j]+1
        current.n <- n.simulated[current.index]
        
        current.lgam <- simulated.lgam[current.index]
        current.ks <- simulated.ks[current.index]
        current.nlgam <- simulated.nlgam[current.index]
        if(test.statistic=="dcov"){
          normalized.lgam <-current.lgam*current.n
          normalized.ks <-current.ks*current.n
          normalized.nlgam <-current.nlgam*current.n
          
          ## Upper-tail p-values
          pval.lgam<-sum(normalized.lgam>=normalized.lgam[1])/current.M
          pval.ks<-sum(normalized.ks>=normalized.ks[1])/current.M
          pval.nlgam <-sum(normalized.nlgam>=normalized.nlgam[1])/current.M
          
          lgam.p[i,j] <- pval.lgam
          ks.p[i,j] <- pval.ks
          nlgam.p[i,j] <- pval.nlgam
          
          if(pval.lgam<=alpha){lgam.rej[j] <- lgam.rej[j]+1}
          if(pval.ks<=alpha){ks.rej[j] <- ks.rej[j]+1}
          if(pval.nlgam<=alpha){nlgam.rej[j] <- nlgam.rej[j]+1}
        }else{
          uncentered.lgam <-current.lgam*sqrt(current.n)
          uncentered.ks <-current.ks*sqrt(current.n)
          uncentered.nlgam <-current.nlgam*sqrt(current.n)

          ## p-values for uncentered statistics
          test.rank <- rank(uncentered.lgam)[1]
          pval.lgam.uncentered <- min(1,2*min(test.rank,current.M-test.rank+1)/current.M)
          
          test.rank <- rank(uncentered.ks)[1]
          pval.ks.uncentered <- min(1, 2*min(test.rank,current.M-test.rank+1)/current.M)
          
          test.rank <- rank(uncentered.nlgam)[1]
          pval.nlgam.uncentered <- min(1, 2*min(test.rank,current.M-test.rank+1)/current.M)
          
          lgam.p.uncentered[i,j] <-pval.lgam.uncentered
          ks.p.uncentered[i,j] <- pval.ks.uncentered
          nlgam.p.uncentered[i,j] <- pval.nlgam.uncentered
          
          if(pval.lgam.uncentered<=alpha){lgam.rej.uncentered[j] <-lgam.rej.uncentered[j]+1}
          if(pval.ks.uncentered<=alpha){ks.rej.uncentered[j] <- ks.rej.uncentered[j]+1}
          if(pval.nlgam.uncentered<=alpha){nlgam.rej.uncentered[j] <- nlgam.rej.uncentered[j]+1}
        }
      }
    }
  }
  ## Return results ##
  if(test.statistic=="dcov"){
    lgam.rej <- lgam.rej/n.trial
    ks.rej <- ks.rej/n.trial
    nlgam.rej <- nlgam.rej/n.trial
    
    result <- list(
      N.perm=N.perm.compare,
      lgam.rej=lgam.rej,
      ks.rej=ks.rej,
      nlgam.rej=nlgam.rej,
      lgam.p=lgam.p,
      ks.p=ks.p,
      nlgam.p=nlgam.p,
      simulated.lgam=simulated.lgam.list,
      simulated.ks=simulated.ks.list,
      simulated.nlgam=simulated.nlgam.list,
      n.simulated=n.simulated.list
    )
    
  }else{
    lgam.rej.uncentered <- lgam.rej.uncentered/n.trial
    ks.rej.uncentered <- ks.rej.uncentered/n.trial
    nlgam.rej.uncentered <- nlgam.rej.uncentered/n.trial
    
    result <- list(
      N.perm=N.perm.compare,

      lgam.rej.uncentered=lgam.rej.uncentered,
      ks.rej.uncentered=ks.rej.uncentered,
      nlgam.rej.uncentered=nlgam.rej.uncentered,

      lgam.p.uncentered=lgam.p.uncentered,
      ks.p.uncentered=ks.p.uncentered,
      nlgam.p.uncentered=nlgam.p.uncentered,
      
      simulated.lgam=simulated.lgam.list,
      simulated.ks=simulated.ks.list,
      simulated.nlgam=simulated.nlgam.list,
      n.simulated=n.simulated.list
    )
  }
  
  return(result)
}
