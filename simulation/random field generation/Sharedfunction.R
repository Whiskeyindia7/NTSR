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

simulation.test<-function(n.trial=2000, N.perm=499, n.obs=100, test.statistic="cov", true.cor="L",
                          window.bound=1, window=NULL, test.target="sig", dcov.type="V",
                          scenario="SE1", test.method="var", radius=window.bound/2, alpha=0.05){
  #
  # INPUT DESCRIPTION
  # n.trial ... the number of repetition for whole test
  # N.perm ... the number of permutation (random shift)
  # n.obs ... the number of observation in window
  # test.statistic ... the test statistic. "cov", "tau" and "dcov" can be used
  # true.cor ... true correlation structure between y and C. "L" is linear, "S" is squared and "E" means exponential 
  # window.bound ... the length of side of window. default is 1 (unit square). useless if window is specified.
  # test.target ... power ("pow") of significance level ("sig")? default is significance test 
  # window ... specified window. if you want to specify the observational window.
  # scenario ... scenario. SE1, NS, LN, N, E1, SE4 can be used
  # test.method ... variance correction ("var") and torrus correction ("tor")
  # radius ... radius for shift vector
  #
  lgam.rej<-ks.rej<-nlgam.rej<-0
  lgam.p<-ks.p<-nlgam.p<-c()
  num<-1
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
      test.points <- runifpoint(n.obs, win=window)   ## Observed locations
    }else{window<-window}
    
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
    else if(scenario=="SE4"){
      auxC <- attr(rLGCP("exp", mu=0, var=4, scale=0.2, saveLambda=TRUE, win=window),"Lambda") # Random effect 
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
    bws<-npregbw(response ~ covariate1, data=data, regtype="lc", ckertype="epanechnikov")
    ksfit<-npreg(bws, residuals=T)
    ksres<-ksfit$resid
    
    ## Method: Nonlinear GAM
    nlgammod<-gam(response ~ s(covariate1) + s(x,y), data=data)
    nlgam.pred<-predict.gam(nlgammod, type="terms")
    nlgamres<-nlgammod$residuals + as.vector(nlgam.pred[,"s(x,y)"])
    
    if(test.method=="tor"){
      
      # test statistic
      if(test.statistic=="cov"){
        stat.lgam<-cov(lgamres, data$covariate2)
        stat.ks<-cov(ksres, data$covariate2)
        stat.nlgam<-cov(nlgamres, data$covariate2)
      }else if(test.statistic=="tau"){
        stat.lgam<-cor(lgamres, data$covariate2, method="kendall")
        stat.ks<-cor(ksres, data$covariate2, method="kendall")
        stat.nlgam<-cor(nlgamres, data$covariate2, method="kendall")
      }else if(test.statistic=="dcov"){
        stat.lgam<-dcov(lgamres, data$covariate2, type=dcov.type)*n.obs*T2(lgamres, data$covariate2)
        stat.ks<-dcov(ksres, data$covariate2, type=dcov.type)*n.obs*T2(ksres, data$covariate2)
        stat.nlgam<-dcov(nlgamres, data$covariate2, type=dcov.type)*n.obs*T2(nlgamres, data$covariate2)
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
          simulated.lgam[k+1] <- cor(cov.interest,lgamres,method="kendall")
          simulated.ks[k+1] <- cor(cov.interest,ksres,method="kendall")
          simulated.nlgam[k+1] <- cor(cov.interest,nlgamres,method="kendall")
        }else if(test.statistic=="dcov"){
          simulated.lgam[k+1] <- dcov(cov.interest,lgamres,type=dcov.type)*n.obs*T2(lgamres, cov.interest)
          simulated.ks[k+1] <- dcov(cov.interest,ksres,type=dcov.type)*n.obs*T2(ksres, cov.interest)
          simulated.nlgam[k+1] <- dcov(cov.interest,nlgamres,type=dcov.type)*n.obs*T2(nlgamres, cov.interest)
        }
      }
      
    }else if(test.method=="var"){
      # test statistic
      if(test.statistic=="cov"){
        stat.lgam<-cov(lgamres, data$covariate2)
        stat.ks<-cov(ksres, data$covariate2)
        stat.nlgam<-cov(nlgamres, data$covariate2)
      }else if(test.statistic=="tau"){
        stat.lgam<-cor(lgamres, data$covariate2, method="kendall")
        stat.ks<-cor(ksres, data$covariate2, method="kendall")
        stat.nlgam<-cor(nlgamres, data$covariate2, method="kendall")
      }else if(test.statistic=="dcov"){
        stat.lgam<-dcov(lgamres, data$covariate2, type=dcov.type)*T2(lgamres, data$covariate2)
        stat.ks<-dcov(ksres, data$covariate2, type=dcov.type)*T2(ksres, data$covariate2)
        stat.nlgam<-dcov(nlgamres, data$covariate2, type=dcov.type)*T2(nlgamres, data$covariate2)
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
        jump <- runifdisc(1, radius=radius)
        test.points.shifted <- shift(test.points, c(jump$x,jump$y))
        W.reduced <- intersect.owin(test.points$window, test.points.shifted$window)
        test.points.reduced <- test.points[W.reduced]
        test.points.reduced.backshifted <- shift(test.points.reduced, c(-jump$x,-jump$y))
        
        all.points <- cbind(test.points$x, test.points$y)
        reduced.points <- cbind(test.points.reduced$x, test.points.reduced$y)
        reduced.indices <- match(apply(reduced.points, 1, paste, collapse = ","), apply(all.points, 1, paste, collapse = ","))
        
        lgamres.reduced<-lgamres[reduced.indices]
        ksres.reduced<-ksres[reduced.indices]
        nlgamres.reduced<-nlgamres[reduced.indices]
        
        vals2 <- Z2[test.points.reduced.backshifted]
        
        if(test.statistic=="cov"){
          simulated.lgam[k+1] <- cov(lgamres.reduced, vals2)
          simulated.ks[k+1] <- cov(ksres.reduced, vals2)
          simulated.nlgam[k+1] <- cov(nlgamres.reduced, vals2)
          n.simulated[k+1] <- test.points.reduced$n
        }else if(test.statistic=="tau"){
          simulated.lgam[k+1] <- cor(lgamres.reduced, vals2, method="kendall")
          simulated.ks[k+1] <- cor(ksres.reduced, vals2, method="kendall")
          simulated.nlgam[k+1] <- cor(nlgamres.reduced, vals2, method="kendall")
          n.simulated[k+1] <- test.points.reduced$n
        }else if(test.statistic=="dcov"){
          simulated.lgam[k+1] <- dcov(lgamres.reduced, vals2, type=dcov.type)*T2(lgamres.reduced, vals2)
          simulated.ks[k+1] <- dcov(ksres.reduced, vals2, type=dcov.type)*T2(ksres.reduced, vals2)
          simulated.nlgam[k+1] <- dcov(nlgamres.reduced, vals2, type=dcov.type)*T2(nlgamres.reduced, vals2)
          n.simulated[k+1] <- test.points.reduced$n
        }
      }
      
      # Variance correction with asymptotic variance
      if(test.statistic=="dcov"){
        simulated.lgam <- (simulated.lgam)*(n.simulated)
        simulated.ks <- (simulated.ks)*(n.simulated)
        simulated.nlgam <- (simulated.nlgam)*(n.simulated)
      }
      else{
        simulated.lgam <- (simulated.lgam - mean(simulated.lgam))*sqrt(n.simulated)
        simulated.ks <- (simulated.ks - mean(simulated.ks))*sqrt(n.simulated)
        simulated.nlgam <- (simulated.nlgam - mean(simulated.nlgam))*sqrt(n.simulated)
      }
      
      
    }
    
    ## p-value
    ## For linear GAM
    if(test.statistic=="dcov"){
      pval.lgam<-sum(simulated.lgam >= simulated.lgam[1]) / (N.perm + 1)
    }
    else{
      test.rank <- rank(simulated.lgam)[1]
      pval.lgam<- 2*min(test.rank, N.perm+1-test.rank)/(N.perm+1)
    }
    lgam.p[i]<-pval.lgam
    if(pval.lgam < alpha){lgam.rej <- lgam.rej+1}
    
    ## For Nadaraya-Watson estimator 
    if(test.statistic=="dcov"){
      pval.ks<-sum(simulated.ks >= simulated.ks[1]) / (N.perm + 1)
    }
    else{
      test.rank <- rank(simulated.ks)[1]
      pval.ks<- 2*min(test.rank, N.perm+1-test.rank)/(N.perm+1)
    }
    ks.p[i]<-pval.ks
    if(pval.ks < alpha){ks.rej <- ks.rej+1}
    
    ## For nonlinear GAM
    if(test.statistic=="dcov"){
      pval.nlgam<-sum(simulated.nlgam >= simulated.nlgam[1]) / (N.perm + 1)
    }
    else{
      test.rank <- rank(simulated.nlgam)[1]
      pval.nlgam<- 2*min(test.rank, N.perm+1-test.rank)/(N.perm+1)
    }
    nlgam.p[i]<-pval.nlgam
    if(pval.nlgam < alpha){nlgam.rej <- nlgam.rej+1}
    num<-num+1
  }
  lgam.rej<-lgam.rej/n.trial
  ks.rej<-ks.rej/n.trial
  nlgam.rej<-nlgam.rej/n.trial
  
  result<-list(lgam.rej=lgam.rej,
               ks.rej=ks.rej,
               nlgam.rej=nlgam.rej,
               lgam.p=lgam.p,
               ks.p=ks.p,
               nlgam.p=nlgam.p)
  
  return(result)
}
