## Set environment##
rm(list = ls())
setwd('/home/whiskeyindia7/revised_rs')
source('functions_variousK.R')

## Packages
library(spatstat);library(fields);library(mvtnorm);library(ggplot2);library(geoR);library(mgcv);library(np)
library(fields); library(classInt); library(nlme); library(convoSPAT); library(dcov);


n.trial<-2000; n.obs<-100; seed<-123; N.perm<-c(99,299,499,799,999);scenario<-"HT";


## variance correction + sample covariance
true.cor<-"L";test.target<-"sig";test.method<-"var";test.statistic<-'cov'
HT.L.sig.var.cov<-simulation.test.variousK(n.trial=n.trial, n.obs=n.obs, test.statistic=test.statistic, N.perm=N.perm,
                               true.cor=true.cor, test.target=test.target, scenario=scenario, test.method=test.method, seed=seed)

true.cor<-"S";test.target<-"sig";test.method<-"var";test.statistic<-'cov'
HT.S.sig.var.cov<-simulation.test.variousK(n.trial=n.trial, n.obs=n.obs, test.statistic=test.statistic, N.perm=N.perm,
                               true.cor=true.cor, test.target=test.target, scenario=scenario, test.method=test.method, seed=seed)

true.cor<-"L";test.target<-"sig";test.method<-"var";test.statistic<-'dcov'
HT.L.sig.var.dcov<-simulation.test.variousK(n.trial=n.trial, n.obs=n.obs, test.statistic=test.statistic, N.perm=N.perm,
                               true.cor=true.cor, test.target=test.target, scenario=scenario, test.method=test.method, seed=seed)

true.cor<-"S";test.target<-"sig";test.method<-"var";test.statistic<-'dcov'
HT.S.sig.var.dcov<-simulation.test.variousK(n.trial=n.trial, n.obs=n.obs, test.statistic=test.statistic, N.perm=N.perm,
                               true.cor=true.cor, test.target=test.target, scenario=scenario, test.method=test.method, seed=seed)

save(HT.L.sig.var.cov,HT.S.sig.var.cov,
     HT.L.sig.var.dcov,HT.S.sig.var.dcov,
     file=paste0("var_", scenario, ".RData"))
