rm(list=ls())
source('iterative_rs.R')
set.seed(1234)

load("Chlorophyll.RData")
X<-as.matrix(X)


# Constant 
n.obs<-test.points$n;
N.perm<-999;alpha<-0.05;radius<-0.30;



## Examples with covariance + variance correction
# Examples: linear GAM
reg.model<-"lin.GAM";test.method<-"var";test.statistic<-"cov"
test.lin <- iterative.test(X=X, X.fields=X.fields, y=y, reg.model=reg.model, test.method=test.method, radius=radius,N.perm=N.perm,
                           test.statistic=test.statistic, alpha=alpha, test.points=test.points)
# Examples: Nadaraya-Watson
reg.model<-"NW";test.method<-"var";test.statistic<-"cov"
test.nw <- iterative.test(X=X, X.fields=X.fields, y=y, reg.model=reg.model, test.method=test.method, radius=radius,N.perm=N.perm,
                          test.statistic=test.statistic, alpha=alpha, test.points=test.points)
# Examples: nonlinear GAM
reg.model<-"nl.GAM";test.method<-"var";test.statistic<-"cov"
test.nl <- iterative.test(X=X, X.fields=X.fields, y=y, reg.model=reg.model, test.method=test.method, radius=radius,N.perm=N.perm,
                          test.statistic=test.statistic, alpha=alpha, test.points=test.points)

# Examples: linear GAM
reg.model<-"lin.GAM";test.method<-"var";test.statistic<-"dcov"
test.lin.dcov <- iterative.test(X=X, X.fields=X.fields, y=y, reg.model=reg.model, test.method=test.method, radius=radius,N.perm=N.perm,
                                test.statistic=test.statistic, alpha=alpha, test.points=test.points)
# Examples: Nadaraya-Watson
reg.model<-"NW";test.method<-"var";test.statistic<-"dcov"
test.nw.dcov <- iterative.test(X=X, X.fields=X.fields, y=y, reg.model=reg.model, test.method=test.method, radius=radius,N.perm=N.perm,
                               test.statistic=test.statistic, alpha=alpha, test.points=test.points)
# Examples: nonlinear GAM
reg.model<-"nl.GAM";test.method<-"var";test.statistic<-"dcov"
test.nl.dcov <- iterative.test(X=X, X.fields=X.fields, y=y, reg.model=reg.model, test.method=test.method, radius=radius,N.perm=N.perm,
                               test.statistic=test.statistic, alpha=alpha, test.points=test.points)

test.lin;test.nw;test.nl
test.lin.dcov;test.nw.dcov;test.nl.dcov

save(test.lin,test.nw,test.nl,
     test.lin.dcov, test.nw.dcov, test.nl.dcov,
     file="rs_b.RData")
