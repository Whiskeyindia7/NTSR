rm(list=ls())
setwd('/Users/kanghyunwi/Desktop/revised_rs2')
source('iterative_rs.R')
#set.seed(1234)
n.seed<-1234
load("/Users/kanghyunwi/Downloads/Chlorophyll2.RData")
X<-as.matrix(X)

# Constant 
n.obs<-test.points$n;
N.perm<-999;
alpha<-0.05;radius<-inradius(test.points$window);

# Examples: linear GAM
reg.model<-"lin.GAM";test.method<-"var";test.statistic<-"dcov"
test.lin.dcov <- iterative.test(X=X, X.fields=X.fields, y=y, reg.model=reg.model, test.method=test.method, radius=radius,
                                test.statistic=test.statistic, alpha=alpha, test.points=test.points, seed=n.seed, N.perm=N.perm)
# Examples: Nadaraya-Watson
reg.model<-"NW";test.method<-"var";test.statistic<-"dcov"
test.nw.dcov <- iterative.test(X=X, X.fields=X.fields, y=y, reg.model=reg.model, test.method=test.method, radius=radius,
                               test.statistic=test.statistic, alpha=alpha, test.points=test.points, seed=n.seed, N.perm=N.perm)
# Examples: nonlinear GAM
reg.model<-"nl.GAM";test.method<-"var";test.statistic<-"dcov"
test.nl.dcov <- iterative.test(X=X, X.fields=X.fields, y=y, reg.model=reg.model, test.method=test.method, radius=radius,
                               test.statistic=test.statistic, alpha=alpha, test.points=test.points, seed=n.seed, N.perm=N.perm)

# Examples: linear GAM
reg.model<-"lin.GAM";test.method<-"var";test.statistic<-"cov.uncentered"
test.lin.uncentered <- iterative.test(X=X, X.fields=X.fields, y=y, reg.model=reg.model, test.method=test.method, radius=radius,
                                test.statistic=test.statistic, alpha=alpha, test.points=test.points, seed=n.seed, N.perm=N.perm)
# Examples: Nadaraya-Watson
reg.model<-"NW";test.method<-"var";test.statistic<-"cov.uncentered"
test.nw.uncentered  <- iterative.test(X=X, X.fields=X.fields, y=y, reg.model=reg.model, test.method=test.method, radius=radius,
                               test.statistic=test.statistic, alpha=alpha, test.points=test.points, seed=n.seed, N.perm=N.perm)
# Examples: nonlinear GAM
reg.model<-"nl.GAM";test.method<-"var";test.statistic<-"cov.uncentered"
test.nl.uncentered <- iterative.test(X=X, X.fields=X.fields, y=y, reg.model=reg.model, test.method=test.method, radius=radius,
                               test.statistic=test.statistic, alpha=alpha, test.points=test.points, seed=n.seed, N.perm=N.perm)

test.lin.dcov[[length(test.lin.dcov)]]; test.nw.dcov[[length(test.nw.dcov)]]; test.nl.dcov[[length(test.nl.dcov)]]
test.lin.uncentered[[length(test.lin.uncentered)]]; test.nw.uncentered[[length(test.nw.uncentered)]]; test.nl.uncentered[[length(test.nl.uncentered)]]


if(F){
  save(test.lin.dcov, test.nw.dcov, test.nl.dcov,
       test.lin.uncentered, test.nw.uncentered, test.nl.uncentered,
       file="rs_b.RData") 
}
