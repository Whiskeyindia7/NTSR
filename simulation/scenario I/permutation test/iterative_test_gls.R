library(gstat)
gaussian.cov<-function(d, phi, partial.sill=1,nugget){
  C<-matrix(NA, ncol=dim(d)[1], nrow=dim(d)[1])
  idx1<-which(d==0)
  C[idx1]<-partial.sill+nugget
  idx2<-which(d>0)
  C[idx2]<-partial.sill*exp(-(d[idx2]/phi)^2)
  return(C)
}

iterative.gls <- function(X, y, test.points, cov.interest, tol=0.1, max.iter=10, is.full=TRUE){
  loc <- cbind(test.points$x, test.points$y)
  d <- as.matrix(dist(loc)); n.obs<-test.points$n  
  dist.vec<-c();b<-list();M<-list();M.inv<-list() 
  
  X<-cbind(1,X)
  lm.fit<-lm(y~X)
  lm.res<-lm.fit$residuals
  
  V<-variog(coords=loc, data=lm.res, message=F)
  S2<-mean((lm.res)^2)
  V.fit<-variofit(V, ini.cov.pars=c(2, 0.2), cov.model="gaussian", fix.nugget=F, nugget=1,
                  messages=F, limits=pars.limits(phi=c(1e-8,0.5), sigmasq=c(1e-8,Inf)))
  Vmat<-gaussian.cov(d=d, phi=V.fit$cov.pars[2], partial.sill=V.fit$cov.pars[1], nugget=V.fit$nugget)
  
  V.inv<-chol2inv(chol(Vmat))
  b.gls<-solve(t(X)%*%V.inv%*%X)%*%t(X)%*%V.inv%*%y
  gls.res<-as.numeric(y-X%*%b.gls)
  D <- matrix(NA, nrow=n.obs, ncol=max.iter)
  D[,1] <- gls.res; b[[1]] <- b.gls; dist.vec[1]<-sum(abs(lm.res - gls.res)); M[[1]]<-Vmat; M.inv[[1]]<-V.inv
  if(mean(abs(lm.res - D[,1]))<tol){
    R<-list(beta=b.gls, Vmat=Vmat, V.inv=V.inv, residuals=gls.res)
    return(R)
  }
  # Again...
  S2<-mean((gls.res)^2)
  V<-variog(coords=loc, data=gls.res, message=F)
  V.fit<-variofit(V, ini.cov.pars=c(2, 0.2), cov.model="gaussian", fix.nugget=F, nugget=1,
                  messages=F, limits=pars.limits(phi=c(1e-8,0.5), sigmasq=c(1e-8,Inf)))
  Vmat<-gaussian.cov(d=d, phi=V.fit$cov.pars[2], partial.sill=V.fit$cov.pars[1], nugget=V.fit$nugget)
  V.inv<-chol2inv(chol(Vmat)); b.gls<-solve(t(X)%*%V.inv%*%X)%*%t(X)%*%V.inv%*%y
  gls.res<-as.numeric(y-X%*%b.gls)
  D[,2]<-gls.res; dist.vec[2]<-sum(abs(D[,1] - D[,2])); b[[2]] <- b.gls; M[[2]]<-Vmat; M.inv[[2]]<-V.inv

  if(mean(abs(D[,1] - D[,2]))<tol){
    R<-list(beta=b.gls, Vmat=Vmat, V.inv=V.inv, residuals=gls.res)
    return(R)
  }else{
    i <- 2
    while(mean(abs(D[,i] - D[,i-1])) > tol){
      S2<-mean((D[,i])^2)
      V<-variog(coords=loc, data=D[,i], message=F)
      V.fit<-variofit(V, ini.cov.pars=c(2, 0.2), cov.model="gaussian", fix.nugget=F, nugget=1,
                      messages=F, limits=pars.limits(phi=c(1e-8,0.5), sigmasq=c(1e-8,Inf)))
      Vmat<-gaussian.cov(d=d, phi=V.fit$cov.pars[2], partial.sill=V.fit$cov.pars[1], nugget=V.fit$nugget)
      V.inv<-chol2inv(chol(Vmat)); b.gls<-solve(t(X)%*%V.inv%*%X)%*%t(X)%*%V.inv%*%y
      gls.res<-as.numeric(y-X%*%b.gls)
      D[,i+1] <- gls.res; dist.vec[i+1]<-sum(abs(D[,i] - D[,i+1])); b[[i+1]] <- b.gls; M[[i+1]]<-Vmat; M.inv[[i+1]]<-V.inv
      i <- i + 1
      if(i==max.iter){
        min.idx <- which.min(dist.vec)
        R<-list(beta=b[[min.idx]], Vmat=M[[min.idx]], V.inv=M.inv[[min.idx]], residuals=D[,min.idx])
        return(R)
      }
    }
    R<-list(beta=b.gls, Vmat=Vmat, V.inv=V.inv, residuals=gls.res)
    return(R)
  }
}

iterative.test<-function(X, y, test.points, N.perm=499, alpha=0.05){
  # constants
  iter<-list()
  n.obs<-length(y)
  pval<-matrix(Inf,nrow=1,ncol=ncol(X))
  rej<-matrix(0,nrow=1, ncol = ncol(X))
  Cnames<-names.of.cov<-colnames(X)
  colnames(rej)<-colnames(pval)<-colnames(X)
  loc<-cbind(test.points$x, test.points$y)
  t<-1
  # Testing
  while(any(pval > alpha)){
    for(i in 1:ncol(X)){
      X.reduce<-cbind(X[, setdiff(colnames(X), names.of.cov[i]), drop = FALSE])
      p<-ncol(as.matrix(X.reduce))
      iter.fit<-iterative.gls(X=X, y=y, test.points=test.points, is.full=TRUE)
      b<-iter.fit$beta; Vmat<-iter.fit$Vmat; res<-iter.fit$residuals; V.inv<-iter.fit$V.inv
      
      ## 3. Obtain T0
      A.iter<-(b[names.of.cov[i],])^2
      B.iter<-solve(t(X) %*% V.inv %*% X)[i,i]
      
      simulated <- numeric(N.perm + 1)
      simulated[1]<-(A.iter/B.iter)

      # Fit reduced model
      # Iterative GLS 
      iter.fit<-iterative.gls(X=X.reduce, y=y, test.points=test.points, cov.interest=names.of.cov[i])
      b<-iter.fit$beta; Vmat<-iter.fit$Vmat; res<-iter.fit$residuals; V.inv<-iter.fit$V.inv
      
      ## Obtain filtered residuals
      filter<-sqrtm(V.inv)%*%res
      refilter<-sqrtm(Vmat)
      
      for (k in 1:N.perm) {
        ## permute the filtered residuals
        permute <- sample(n.obs)
        perm.filter<-filter[permute]
        
        ## Reconstruct Y* using permuted filtered residuals and trend of nuisance
        perm.y<- cbind(1,X.reduce) %*% b + as.numeric(refilter%*%perm.filter)
        ## Fit the full model again : Y* ~ C1 + C2
        if(k==N.perm){print("finish")}
        iter.fit<-iterative.gls(X=X, y=perm.y, test.points=test.points, is.full=T)
        b.f<-iter.fit$beta; Vmat<-iter.fit$Vmat; res<-iter.fit$residuals; V.inv<-iter.fit$V.inv
        
        A.iter<-(b.f[names.of.cov[i],])^2
        B.iter<-solve(t(X) %*% V.inv %*% X)[i,i]
        
        simulated[k+1]<-(A.iter / B.iter)
      }  # end for permutation

      # p-value
      test.rank <- rank(simulated)[1]
      pval[,i]<- 2*min(test.rank, N.perm+1-test.rank)/(N.perm+1)
    }
    #out<-list(pval=pval, aic=AIC, bic=BIC)
    out<-list(pval=pval)
    iter[[t]]<-out
    
    if(max(pval) >= alpha){
      if(p==1){
        print("No covariates are significant")
        rej<-matrix(0, nrow=1, ncol = length(Cnames))
        colnames(rej)<-Cnames
        iter[[t+1]]<-rej
        return(iter)
      }
      else{
        M<-which.max(pval)
        X<-as.matrix(X[,-M]);names.of.cov<-names.of.cov[-M]
        pval<-matrix(Inf,nrow=1,ncol=ncol(X));
        colnames(pval)<-names.of.cov
        t<-t+1
      }
    }
    else{
      m<- match(colnames(pval), colnames(rej))
      rej[,m]<-1
      iter[[t+1]]<-rej
      return(iter)
    }
  }
  m<- match(colnames(pval), colnames(rej))
  rej[,m]<-1;q<-m[pval>=alpha];rej[,q]<-0
  iter[[t+1]]<-rej
  return(iter)
}
