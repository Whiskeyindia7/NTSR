## Set environment##
rm(list = ls())
set.seed(123)
source('Sharedfunction.R')

## Packages
library(spatstat);library(fields);library(mvtnorm);library(ggplot2);library(geoR);library(mgcv);library(np)
library(fields); library(classInt); library(nlme); library(convoSPAT); library(dcov);


n.trial<-1000; N.perm<-499; n.obs<-100; dcov.type<-"V"
test.statistic<-'dcov'

# Torrus correction
true.cor<-"L";test.target<-"sig";test.method<-"tor";scenario<-"SE1"
SE1.L.sig.tor<-simulation.test(n.trial=n.trial, N.perm=N.perm, n.obs=n.obs, 
                               test.statistic=test.statistic, true.cor=true.cor,
                               test.target=test.target, dcov.type=dcov.type, 
                               scenario=scenario, test.method=test.method)

true.cor<-"S";test.target<-"sig";test.method<-"tor";scenario<-"SE1"
SE1.S.sig.tor<-simulation.test(n.trial=n.trial, N.perm=N.perm, n.obs=n.obs, 
                               test.statistic=test.statistic, true.cor=true.cor,
                               test.target=test.target, dcov.type=dcov.type, 
                               scenario=scenario, test.method=test.method)

true.cor<-"E";test.target<-"sig";test.method<-"tor";scenario<-"SE1"
SE1.E.sig.tor<-simulation.test(n.trial=n.trial, N.perm=N.perm, n.obs=n.obs, 
                               test.statistic=test.statistic, true.cor=true.cor,
                               test.target=test.target, dcov.type=dcov.type, 
                               scenario=scenario, test.method=test.method)


true.cor<-"L";test.target<-"pow";test.method<-"tor";scenario<-"SE1"
SE1.L.pow.tor<-simulation.test(n.trial=n.trial, N.perm=N.perm, n.obs=n.obs, 
                               test.statistic=test.statistic, true.cor=true.cor,
                               test.target=test.target, dcov.type=dcov.type, 
                               scenario=scenario, test.method=test.method)

true.cor<-"S";test.target<-"pow";test.method<-"tor";scenario<-"SE1"
SE1.S.pow.tor<-simulation.test(n.trial=n.trial, N.perm=N.perm, n.obs=n.obs, 
                               test.statistic=test.statistic, true.cor=true.cor,
                               test.target=test.target, dcov.type=dcov.type, 
                               scenario=scenario, test.method=test.method)

true.cor<-"E";test.target<-"pow";test.method<-"tor";scenario<-"SE1"
SE1.E.pow.tor<-simulation.test(n.trial=n.trial, N.perm=N.perm, n.obs=n.obs, 
                               test.statistic=test.statistic, true.cor=true.cor,
                               test.target=test.target, dcov.type=dcov.type, 
                               scenario=scenario, test.method=test.method)


## variance correction
true.cor<-"L";test.target<-"sig";test.method<-"var";scenario<-"SE1"
SE1.L.sig.var<-simulation.test(n.trial=n.trial, N.perm=N.perm, n.obs=n.obs, 
                               test.statistic=test.statistic, true.cor=true.cor,
                               test.target=test.target, dcov.type=dcov.type, 
                               scenario=scenario, test.method=test.method)

true.cor<-"S";test.target<-"sig";test.method<-"var";scenario<-"SE1"
SE1.S.sig.var<-simulation.test(n.trial=n.trial, N.perm=N.perm, n.obs=n.obs, 
                               test.statistic=test.statistic, true.cor=true.cor,
                               test.target=test.target, dcov.type=dcov.type, 
                               scenario=scenario, test.method=test.method)

true.cor<-"E";test.target<-"sig";test.method<-"var";scenario<-"SE1"
SE1.E.sig.var<-simulation.test(n.trial=n.trial, N.perm=N.perm, n.obs=n.obs, 
                               test.statistic=test.statistic, true.cor=true.cor,
                               test.target=test.target, dcov.type=dcov.type, 
                               scenario=scenario, test.method=test.method)


true.cor<-"L";test.target<-"pow";test.method<-"var";scenario<-"SE1"
SE1.L.pow.var<-simulation.test(n.trial=n.trial, N.perm=N.perm, n.obs=n.obs, 
                               test.statistic=test.statistic, true.cor=true.cor,
                               test.target=test.target, dcov.type=dcov.type, 
                               scenario=scenario, test.method=test.method)

true.cor<-"S";test.target<-"pow";test.method<-"var";scenario<-"SE1"
SE1.S.pow.var<-simulation.test(n.trial=n.trial, N.perm=N.perm, n.obs=n.obs, 
                               test.statistic=test.statistic, true.cor=true.cor,
                               test.target=test.target, dcov.type=dcov.type, 
                               scenario=scenario, test.method=test.method)

true.cor<-"E";test.target<-"pow";test.method<-"var";scenario<-"SE1"
SE1.E.pow.var<-simulation.test(n.trial=n.trial, N.perm=N.perm, n.obs=n.obs, 
                               test.statistic=test.statistic, true.cor=true.cor,
                               test.target=test.target, dcov.type=dcov.type, 
                               scenario=scenario, test.method=test.method)


save(SE1.L.sig.tor,SE1.S.sig.tor,SE1.E.sig.tor,
     SE1.L.pow.tor,SE1.S.pow.tor,SE1.E.pow.tor, 
     file="SE1_tor2.RData")

save(SE1.L.sig.var,SE1.S.sig.var,SE1.E.sig.var,
     SE1.L.pow.var,SE1.S.pow.var,SE1.E.pow.var, 
     file="SE1_var2.RData")