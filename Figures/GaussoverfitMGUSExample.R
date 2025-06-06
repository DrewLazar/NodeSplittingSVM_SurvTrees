rm(list=ls())
# Clear console
cat("\014")
# Force-set directory to this file in RStudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load the needed survsim scripts  
source("simple.surv.sim.int.R")
source("simple.ev.sim.R")
source("bootstrapPruning2.R")
source("DipolarSurvivalTree_PolyKernel.R")
source("DipolarSurvivalTree_GaussKernel.R") 
source("Predictionfunction.R")
source("BrierScoreFun.R")
source("kappacrossvalidationB.R")
source("UnivariateSurvivalTree.R")
library(ggplot2)

#mgus plots 
mgus = read.csv("./Data/mgus.csv")
names(mgus)[names(mgus) == 'time'] <- 'surv'
names(mgus)[names(mgus) == 'event'] <- 'status'
mgus$fac_sex<-as.integer(as.character(mgus$fac_sex)=="female")

alldata = mgus;
alldata=alldata[-1]
alldata=subset(alldata, select=-c(num_pctime,num_alb,num_creat))
alldata=na.omit(alldata)

covariates = c("num_age","num_hgb"); time = "surv"; censor = "status";
alldatas=alldata[alldata$status==1,]
alldatac=alldata[alldata$status==0,]

g0<- ggplot() + 
  geom_point(data=alldatas, aes(x = num_age, y = num_hgb, size = surv, color=surv),alpha=0.4) +
  geom_point(data=alldatac, aes(x = num_age, y = num_hgb,size=surv),show.legend = FALSE,color='green',alpha=0.4) 

X = alldata[, covariates]
distX = c(dist(X))
epsilon = quantile(distX[distX != 0], probs = c(0.23))
quantiles = c(.25,.75); tolerance = 10^-2; 
nsize = 10; n=nrow(alldata)
eta=-3; 
n = nrow(alldata)
fullsubset<-1:n

ksigma=intvarfun(alldata,covariates)
Dipolar.model <- DipolarSurvivalTree_GaussKernel$new(
  alldata, time, censor, covariates,
  quantiles, tolerance, epsilon, kappa=exp(eta), nsize,
  pureweight=1, mixedweight=1, Ksigma =ksigma/100^2 
)
Dipolar.model2 <- DipolarSurvivalTree_GaussKernel$new(
  alldata, time, censor, covariates,
  quantiles, tolerance, epsilon, kappa=exp(eta), nsize,
  pureweight=1, mixedweight=1, Ksigma =ksigma/100
)
Dipolar.model3 <- DipolarSurvivalTree_GaussKernel$new(
  alldata, time, censor, covariates,
  quantiles, tolerance, epsilon, kappa=exp(eta), nsize,
  pureweight=1, mixedweight=1, Ksigma =ksigma
)


dipolargtree<-Dipolar.model$createtree(fullsubset)
dipolargtree1<-Dipolar.model2$createtree(fullsubset)
dipolargtree2<-Dipolar.model3$createtree(fullsubset)

#model 1
X = as.matrix(dipolargtree$data)
mupmdiff = dipolargtree$opt_w0_mupmdiff$mupmdiff
w0 = dipolargtree$opt_w0_mupmdiff$w0

Kern.X = Dipolar.model$K

cc1 <- emdbook::curve3d(t(mupmdiff) %*%  Kern.X(X,c(x,y)) + w0,n = c(250,250), 
                        xlim=c(17.2,87.8), ylim=c(2.9,23.5), sys3d="none")
dimnames(cc1$z) <- list(cc1$x,cc1$y)
mm1 <- reshape2::melt(cc1$z) #eta1

#model 2
mupmdiff = dipolargtree1$opt_w0_mupmdiff$mupmdiff
w0 = dipolargtree1$opt_w0_mupmdiff$w0

Kern.X = Dipolar.model2$K

cc1 <- emdbook::curve3d(t(mupmdiff) %*%  Kern.X(X,c(x,y)) + w0,n = c(250,250), 
                        xlim=c(17.2,87.8), ylim=c(2.9,23.5), sys3d="none")
dimnames(cc1$z) <- list(cc1$x,cc1$y)
mm2 <- reshape2::melt(cc1$z) #eta1

#model 3
mupmdiff = dipolargtree2$opt_w0_mupmdiff$mupmdiff
w0 = dipolargtree2$opt_w0_mupmdiff$w0

Kern.X = Dipolar.model3$K

cc1 <- emdbook::curve3d(t(mupmdiff) %*%  Kern.X(X,c(x,y)) + w0,n = c(250,250), 
                        xlim=c(17.2,87.8), ylim=c(2.9,23.5), sys3d="none")
dimnames(cc1$z) <- list(cc1$x,cc1$y)
mm3 <- reshape2::melt(cc1$z) #eta1


g1<-g0 + geom_contour(data=mm1,aes(x=Var1,y=Var2,z=value),breaks=0,colour="red",linewidth=.70)

g2<-g1 + geom_contour(data=mm2,aes(x=Var1,y=Var2,z=value),breaks=0,colour="navyblue",linewidth=.70)

g3<- g2 + geom_contour(data=mm3,aes(x=Var1,y=Var2,z=value),breaks=0,colour="black",linewidth=.70)+
  xlab("Age") + ylab("HGB Level")