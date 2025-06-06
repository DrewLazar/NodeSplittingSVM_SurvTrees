# Remove preexisting variables from environment
rm(list=ls())
# Clear console
cat("\014")
# Force-set directory to this file in RStudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# Load the alldata data
load("../Data/Remission.rda")
# Load the DipolarSurvivalTree class; all necessary libraries included in the class file

# Load the needed scripts  
source("../ExtraFiles/simple.surv.sim.int.R")
source("../ExtraFiles/simple.ev.sim.R")
source("../ExtraFiles/Predictionfunction.R")
source("../ExtraFiles/kappacrossvalidation.R")
source("../ExtraFiles/PruningScript.R")
source("../Models/DipolarSurvivalTree_PolyKernel.R")
source("../Models/DipolarSurvivalTree_GaussKernel.R")
source("../Models/UnivariateSurvivalTree.R")
library(dplyr)
library(ggplot2)
library(emdbook)

set.seed(26)

# Define some variables for instantiating our model classes 
alldata = Remission; covariates = c("logWBC","TR"); time = "survt"; censor = "status";

X = alldata[, covariates]
distX = c(dist(X))
epsilon = quantile(distX[distX != 0], probs = c(0.10))
quantiles = c(.15,.55); tolerance = 10^-2; kappa = 1
nsize = 10; eta=3; alpha=2.2

n=nrow(alldata)


#Instantiate our model classes 
Dipolar.model <- DipolarSurvivalTree_PolyKernel$new(
  alldata, time, censor, covariates,
  quantiles, tolerance, epsilon, kappa=exp(eta), nsize,
  pureweight=1, mixedweight=1,Kconstant=0,Kpoly_order =1
)

Dipolarq.model <- DipolarSurvivalTree_PolyKernel$new(
  alldata, time, censor, covariates,
  quantiles, tolerance, epsilon, kappa=exp(eta), nsize,
  pureweight=1, mixedweight=1,Kconstant=1,Kpoly_order =2
)
Dipolarg.model <- DipolarSurvivalTree_GaussKernel$new(
  alldata, time, censor, covariates,
  quantiles, tolerance, epsilon, kappa=exp(eta), nsize,
  pureweight=1, mixedweight=1, Ksigma =1
)
Univariate.model <- UnivariateSurvivalTree$new(alldata, time, censor, covariates,nsize)

#Create tree objects 
dipolartree<-Dipolar.model$createtree(1:n)
dipolarqtree<-Dipolarq.model$createtree(1:n)
dipolargtree<-Dipolarg.model$createtree(1:n)
univariatetree<-Univariate.model$createtree(1:n)

#print full trees and log rank statistics  
print(dipolartree,"lrstat")
print(dipolarqtree,"lrstat")
print(dipolargtree,"lrstat")
print(univariatetree,"lrstat")

dipolartree.prune<-bootstrapPruning(dipolartree,Dipolar.model,alpha)[[4]]
dipolarqtree.prune<-bootstrapPruning(dipolarqtree,Dipolarq.model,alpha)[[4]]
dipolargtree.prune<-bootstrapPruning(dipolargtree,Dipolarg.model,alpha)[[4]]
univariatetree.prune<-bootstrapPruning(univariatetree,Univariate.model,alpha)[[4]]


#Plots for Linear Tree ----
rlindex = as.numeric(rownames(dipolartree$Node0r$Node0rl$data))
rrindex = as.numeric(rownames(dipolartree$Node0r$Node0rr$data))
llindex = as.numeric(rownames(dipolartree$Node0l$Node0ll$data))
lrindex = as.numeric(rownames(dipolartree$Node0l$Node0lr$data))
for (i in 1:42){
  if (i %in% rlindex) {
    alldata$termnode[i] = 'rl' 
  } else if (i %in% rrindex) {
    alldata$termnode[i] = 'rr' 
  } else if (i %in% llindex) {
    alldata$termnode[i] = 'll' 
  } else {
    alldata$termnode[i] = 'lr' 
  }
}
g0<- ggplot() + 
  geom_point(data=alldata, aes(x = logWBC, y = TR, size = survt, color=termnode),alpha=0.4)+
  scale_color_manual(values=c("blue","red", "cyan3","green"))+labs(color = "Node")
#rootnode 
NodeCurr = dipolartree
X = as.matrix(NodeCurr$data)
mupmdiff = NodeCurr$opt_w0_mupmdiff$mupmdiff
w0 = NodeCurr$opt_w0_mupmdiff$w0

Kern.X = Dipolar.model$K

cc1 <- emdbook::curve3d(t(mupmdiff) %*%  Kern.X(X,c(x,y)) + w0,n = c(250,250), 
                        xlim=c(0,5.0), ylim=c(-.015,1.015), sys3d="none")
dimnames(cc1$z) <- list(cc1$x,cc1$y)
mm1 <- reshape2::melt(cc1$z) #eta1

g1<-g0 + geom_contour(data=mm1,aes(x=Var1,y=Var2,z=value),breaks=0,colour="black")


#first left node  
NodeCurr2 = dipolartree$Node0l
X = as.matrix(NodeCurr2$data)
mupmdiff2 = NodeCurr2$opt_w0_mupmdiff$mupmdiff
w02 = NodeCurr2$opt_w0_mupmdiff$w0

Kern.X = Dipolar.model$K

cc2 <- emdbook::curve3d(t(mupmdiff2) %*%  Kern.X(X,c(x,y)) + w02,n = c(250,250), 
                        xlim=c(0,5),ylim=c(-.015,1.015), sys3d="none")
dimnames(cc2$z) <- list(cc2$x,cc2$y)
mm2 <- reshape2::melt(cc2$z) #eta1


g2<-g1 + geom_contour(data=mm2,aes(x=Var1,y=Var2,z=value),breaks=0,colour="black")

#first right node
NodeCurr3 = dipolartree$Node0r
X = as.matrix(NodeCurr3$data)
mupmdiff3 = NodeCurr3$opt_w0_mupmdiff$mupmdiff
w03 = NodeCurr3$opt_w0_mupmdiff$w0

Kern.X = Dipolar.model$K

cc3 <- emdbook::curve3d(t(mupmdiff3) %*%  Kern.X(X,c(x,y)) + w03,n = c(250,250), 
                        xlim=c(0,5.1),ylim=c(-.015,1.015), sys3d="none")
dimnames(cc3$z) <- list(cc3$x,cc3$y)
mm3 <- reshape2::melt(cc3$z) #eta1


g3lt<-g2 + geom_contour(data=mm3,aes(x=Var1,y=Var2,z=value),breaks=0,colour="black")


#Plots for Quadratic Tree ----
rlindex = as.numeric(rownames(dipolarqtree$Node0r$Node0rl$data))
rrindex = as.numeric(rownames(dipolarqtree$Node0r$Node0rr$data))
llindex = as.numeric(rownames(dipolarqtree$Node0l$Node0ll$data))
lrindex = as.numeric(rownames(dipolarqtree$Node0l$Node0lr$data))
for (i in 1:42){
  if (i %in% rlindex) {
    alldata$termnode[i] = 'rl' 
  } else if (i %in% rrindex) {
    alldata$termnode[i] = 'rr' 
  } else if (i %in% llindex) {
    alldata$termnode[i] = 'll' 
  } else {
    alldata$termnode[i] = 'lr' 
  }
}
g0<- ggplot() + 
  geom_point(data=alldata, aes(x = logWBC, y = TR, size = survt, color=termnode),alpha=0.4)+
  scale_color_manual(values=c("blue","red", "cyan3","green"))+labs(color = "Node")

#rootnode 
NodeCurr = dipolarqtree
X = as.matrix(NodeCurr$data)
mupmdiff = NodeCurr$opt_w0_mupmdiff$mupmdiff
w0 = NodeCurr$opt_w0_mupmdiff$w0

Kern.X = Dipolarq.model$K

cc1 <- emdbook::curve3d(t(mupmdiff) %*%  Kern.X(X,c(x,y)) + w0,n = c(250,250), 
                        xlim=c(0,5), ylim=c(-.015,1.015), sys3d="none")
dimnames(cc1$z) <- list(cc1$x,cc1$y)
mm1 <- reshape2::melt(cc1$z) #eta1

g1<-g0 + geom_contour(data=mm1,aes(x=Var1,y=Var2,z=value),breaks=0,colour="black")


#first left node  
NodeCurr2 = dipolarqtree$Node0l
X = as.matrix(NodeCurr2$data)
mupmdiff2 = NodeCurr2$opt_w0_mupmdiff$mupmdiff
w02 = NodeCurr2$opt_w0_mupmdiff$w0

Kern.X = Dipolarq.model$K

cc2 <- emdbook::curve3d(t(mupmdiff2) %*%  Kern.X(X,c(x,y)) + w02,n = c(250,250), 
                        xlim=c(0,5), ylim=c(-.015,1.015), sys3d="none")
dimnames(cc2$z) <- list(cc2$x,cc2$y)
mm2 <- reshape2::melt(cc2$z) #eta1


g2<-g1 + geom_contour(data=mm2,aes(x=Var1,y=Var2,z=value),breaks=0,colour="black")

#first right node
NodeCurr3 = dipolarqtree$Node0r
X = as.matrix(NodeCurr3$data)
mupmdiff3 = NodeCurr3$opt_w0_mupmdiff$mupmdiff
w03 = NodeCurr3$opt_w0_mupmdiff$w0

Kern.X = Dipolarq.model$K

cc3 <- emdbook::curve3d(t(mupmdiff3) %*%  Kern.X(X,c(x,y)) + w03,n = c(250,250), 
                        xlim=c(3,5), ylim=c(-.015,1.015), sys3d="none")
dimnames(cc3$z) <- list(cc3$x,cc3$y)
mm3 <- reshape2::melt(cc3$z) #eta1


g3qt<-g2 + geom_contour(data=mm3,aes(x=Var1,y=Var2,z=value),breaks=0,colour="black")


#Plots for Gaussian Tree ----
rlindex = as.numeric(rownames(dipolargtree$Node0r$Node0rl$data))
rrindex = as.numeric(rownames(dipolargtree$Node0r$Node0rr$data))
llindex = as.numeric(rownames(dipolargtree$Node0l$Node0ll$data))
lrindex = as.numeric(rownames(dipolargtree$Node0l$Node0lr$data))
for (i in 1:42){
  if (i %in% rlindex) {
    alldata$termnode[i] = 'rl' 
  } else if (i %in% rrindex) {
    alldata$termnode[i] = 'rr' 
  } else if (i %in% llindex) {
    alldata$termnode[i] = 'll' 
  } else {
    alldata$termnode[i] = 'lr' 
  }
}
g0<- ggplot() + 
  geom_point(data=alldata, aes(x = logWBC, y = TR, size = survt, color=termnode),alpha=0.4)+
  scale_color_manual(values=c("blue","red", "cyan3","green"))+labs(color = "Node")
#rootnode 
NodeCurr = dipolargtree
X = as.matrix(NodeCurr$data)
mupmdiff = NodeCurr$opt_w0_mupmdiff$mupmdiff
w0 = NodeCurr$opt_w0_mupmdiff$w0

Kern.X = Dipolarg.model$K

cc1 <- emdbook::curve3d(t(mupmdiff) %*%  Kern.X(X,c(x,y)) + w0,n = c(250,250), 
                        xlim=c(0,5.0), ylim=c(-.015,1.015), sys3d="none")
dimnames(cc1$z) <- list(cc1$x,cc1$y)
mm1 <- reshape2::melt(cc1$z) #eta1

g1<-g0 + geom_contour(data=mm1,aes(x=Var1,y=Var2,z=value),breaks=0,colour="black")


#first left node  
NodeCurr2 = dipolargtree$Node0l
X = as.matrix(NodeCurr2$data)
mupmdiff2 = NodeCurr2$opt_w0_mupmdiff$mupmdiff
w02 = NodeCurr2$opt_w0_mupmdiff$w0

Kern.X = Dipolarg.model$K

cc2 <- emdbook::curve3d(t(mupmdiff2) %*%  Kern.X(X,c(x,y)) + w02,n = c(250,250), 
                        xlim=c(1,5),ylim=c(-.015,1.015), sys3d="none")
dimnames(cc2$z) <- list(cc2$x,cc2$y)
mm2 <- reshape2::melt(cc2$z) #eta1


g2<-g1 + geom_contour(data=mm2,aes(x=Var1,y=Var2,z=value),breaks=0,colour="black")

#first right node
NodeCurr3 = dipolargtree$Node0r
X = as.matrix(NodeCurr3$data)
mupmdiff3 = NodeCurr3$opt_w0_mupmdiff$mupmdiff
w03 = NodeCurr3$opt_w0_mupmdiff$w0

Kern.X = Dipolarg.model$K

cc3 <- emdbook::curve3d(t(mupmdiff3) %*%  Kern.X(X,c(x,y)) + w03,n = c(250,250), 
                        xlim=c(0,5.1),ylim=c(-.015,1.015), sys3d="none")
dimnames(cc3$z) <- list(cc3$x,cc3$y)
mm3 <- reshape2::melt(cc3$z) #eta1


g3gt<-g2 + geom_contour(data=mm3,aes(x=Var1,y=Var2,z=value),breaks=0,colour="black")


#Plots for Univariate Tree ----
rlindex = as.numeric(rownames(dipolargtree$Node0r$Node0rl$data))
rrindex = as.numeric(rownames(dipolargtree$Node0r$Node0rr$data))
llindex = as.numeric(rownames(dipolargtree$Node0l$Node0ll$data))
lrindex = as.numeric(rownames(dipolargtree$Node0l$Node0lr$data))
for (i in 1:42){
  if (i %in% rlindex) {
    alldata$termnode[i] = 'rl' 
  } else if (i %in% rrindex) {
    alldata$termnode[i] = 'rr' 
  } else if (i %in% llindex) {
    alldata$termnode[i] = 'll' 
  } else {
    alldata$termnode[i] = 'lr' 
  }
}






r.index = as.numeric(rownames(univariatetree$Node0r$data))
lr.index = as.numeric(rownames(univariatetree$Node0l$Node0lr$data))
llr.index = as.numeric(rownames(univariatetree$Node0l$Node0ll$Node0llr$data))
llll.index = as.numeric(rownames(univariatetree$Node0l$Node0ll$Node0lll$Node0llll$data))
lllr.index = as.numeric(rownames(univariatetree$Node0l$Node0ll$Node0lll$Node0lllr$data))
for (i in 1:42){
  if (i %in% r.index) {
    alldata$termnode[i] = 'r' 
  } else if (i %in% lr.index) {
    alldata$termnode[i] = 'lr' 
  } else if (i %in% llr.index) {
    alldata$termnode[i] = 'llr' 
  } else if (i %in% llll.index) {
    alldata$termnode[i] = 'llll' 
  } else {
    alldata$termnode[i] = 'lllr' 
  }
}
g0<- ggplot() + 
  geom_point(data=alldata, aes(x = logWBC, y = TR, size = survt, color=termnode),alpha=0.4)+
  scale_color_manual(values=c("blue","red", "cyan3","green","darkmagenta"))+labs(color = "Node")

#rootnode 
split1=-univariatetree$optv[1]
g1<-g0 + geom_vline(xintercept = split1, linetype="solid", color = "black")

#left node  
split2=-univariatetree$Node0l$optv[1]
g2<-g1 + geom_vline(xintercept =split2, linetype="solid", color = "black")

#left left node  
univariatetree$Node0l$Node0ll$optv
g3<- g2 + geom_segment(aes(x=1.45,xend=split2,y=0.5,yend=0.5))

#lelf left left node
split3=-univariatetree$Node0l$Node0ll$Node0lll$optv[1]
g4<-g3 + geom_segment(aes(x=split3,xend=split3,y=0,yend=0.5))









