library(dplyr)
library(ggplot2)
library(osqp)
# Remove preexisting variables from environment
rm(list=ls())
# Clear console
cat("\014")
# Force-set directory to this file in RStudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load the needed scripts  
source("../ExtraFiles/simple.surv.sim.int.R")
source("../ExtraFiles/simple.ev.sim.R")
source("../ExtraFiles/kappacrossvalidation.R")
source("../Models/DipolarSurvivalTree_PolyKernel.R")
source("../Models/DipolarSurvivalTree_GaussKernel.R")
source("../Models/UnivariateSurvivalTree.R")

#To suppress warnings in ggplot for points outside of plotted area 
options( warn = -1 )

#Linear Hazard, Linear Kernel ----
#LINEAR HAZARD 
#Generate some data with linear hazard 
#Define parameters for Weibull distribution for censoring distribution
set.seed(1234)
pc=1; b0c =-5.3;
anc.cens = pc; beta0.cens = -b0c/pc; 
#Define parameters for Weibull distribution for time-to-event distribution
pt = 1; b=1.2; b0t =-b*6; betat = list(-4*b,4*b)
anc.ev=pt; beta0.ev = -b0t/pt; beta =  lapply(betat,"*",-1/pt)

sim.data <- simple.surv.sim.int(n=180, foltime=100, dist.ev=c('weibull'),
                                anc.ev,beta0.ev,dist.cens=c('weibull'),anc.cens,
                               beta0.cens,z=NULL, beta, x=list(c("normal", 1.2,.4),c("normal", 2.6, .3)))
prop_censored=1-sum(sim.data$status)/nrow(sim.data)
print(prop_censored)

#Create Hazard column 
sim.data['hr']=b0t + betat[[1]]*sim.data['x']+ betat[[2]]*sim.data['x.1']
#Correlation between Hazard and survival times
cor(sim.data['stop'],sim.data['hr'])

#Fit a splitting surface

alldata = sim.data; covariates = c("x","x.1"); time = "stop"; censor = "status";
alldatas=alldata[alldata$status==1,]
alldatac=alldata[alldata$status==0,]

alldata = sim.data; covariates = c("x","x.1"); time = "stop"; censor = "status"
eta=3;

X = alldata[, covariates]
distX = c(dist(X))
epsilon = quantile(distX[distX != 0], probs = c(0.40))
quantiles = c(.35,.65); tolerance = 10^-2; kappa = exp(eta)
nsize = 10

n=nrow(alldata)

#Instantiate our model classes 
Dipolar.model <- DipolarSurvivalTree_PolyKernel$new(
  alldata, time, censor, covariates,
  quantiles, tolerance, epsilon, kappa, nsize,
  pureweight=1, mixedweight=1,Kconstant=0,Kpoly_order =1
)

#Create tree objects 
fullsubset<-1:n
dipolarqtree<-Dipolar.model$createtree(fullsubset)

#print full trees and log rank statistics  
print(dipolarqtree,"lrstat")

X = as.matrix(dipolarqtree$data)
mupmdiff = dipolarqtree$opt_w0_mupmdiff$mupmdiff
w0 = dipolarqtree$opt_w0_mupmdiff$w0

Kern.X = function(X,x) {
  (X%*%x)
}

g0<- ggplot() + 
  geom_point(data=alldatas, aes(x = x, y = x.1, size = stop, color=stop),alpha=0.4) +
  geom_point(data=alldatac, aes(x = x, y = x.1,size=stop),show.legend = FALSE,color='green',alpha=0.4) 

cc1 <- emdbook::curve3d(t(mupmdiff) %*%  Kern.X(X,c(x,y)) + w0,n = c(250,250), 
                       xlim=c(-.1,2.2), ylim=c(1.5,4), sys3d="none")
dimnames(cc1$z) <- list(cc1$x,cc1$y)
mm1 <- reshape2::melt(cc1$z)

#cc2 <- emdbook::curve3d(t(mupmdiff) %*%  Kern.X(X,c(x,y)) + w0,n = c(250,250), 
#                        xlim=c(-.1,2.2), ylim=c(1.5,4), sys3d="none")
#dimnames(cc2$z) <- list(cc2$x,cc2$y)
#mm2 <- reshape2::melt(cc2$z)

#cc3 <- emdbook::curve3d(t(mupmdiff) %*%  Kern.X(X,c(x,y)) + w0,n = c(250,250), 
#                        xlim=c(-.1,2.2), ylim=c(1.5,4), sys3d="none")
#dimnames(cc3$z) <- list(cc3$x,cc3$y)
#mm3 <- reshape2::melt(cc3$z)

g1<-g0 +   geom_abline(slope=1.5, intercept=1,color='black')+ 
  geom_contour(data=mm1,aes(x=Var1,y=Var2,z=value),breaks=0,colour="red") #+
# geom_contour(data=mm2,aes(x=Var1,y=Var2,z=value),breaks=0,colour="deeppink4")+
#  geom_contour(data=mm3, aes(x=Var1,y=Var2,z=value),breaks=0,colour="darkred") 
g1
#ggsave("plot1.jpeg",g1,width=5,height=5)


#Parabolic Hazard, Quadratic Kernel ----
#NONLINEAR HAZARD - Downward turning parabola 
set.seed(1234)
#Generate some data with parabolic hazard y = h(x) = ax^2-2ax+a+1
#Define parameters for Weibull distribution for censoring distribution
pc=1; b0c =-1.8;
anc.cens = pc; beta0.cens = -b0c/pc; 
#Define parameters for Weibull distribution for time-to-event distribution
pt = 1; a=-2;b=3.5; b0t =b*(-a-1); betat = list(b*2*a,b*1,b*(-a))
anc.ev=pt; beta0.ev = -b0t/pt; beta =  lapply(betat,"*",-1/pt)
#Generate data 
sim.data <- simple.surv.sim.int(n=180, foltime=10, dist.ev=c('weibull'),
                                anc.ev,beta0.ev,dist.cens=c('weibull'),anc.cens,
                                beta0.cens,z=NULL, beta, x=list(c("normal", 1,.4),c("normal", 1, .3),c("quad",1)))
prop_censored=1-sum(sim.data$status)/nrow(sim.data)
print(prop_censored)
#Create Hazard column 
sim.data['hr']=b0t + betat[[1]]*sim.data['x']+ betat[[2]]*sim.data['x.1']+betat[[3]]*sim.data['x.2']
#Correlation between Hazard and survival times
cor(sim.data['stop'],sim.data['hr'])

#Fit a splitting surface
alldata = sim.data; covariates = c("x","x.1"); time = "stop"; censor = "status";
alldatas=alldata[alldata$status==1,]
alldatac=alldata[alldata$status==0,]

X = alldata[, covariates]
distX = c(dist(X))
epsilon = quantile(distX[distX != 0], probs = c(0.20))
quantiles = c(.35,.65); tolerance = 10^-2; eta =3 
nsize = 10

n=nrow(alldata)

#Instantiate our model classes 
Dipolar.model <- DipolarSurvivalTree_PolyKernel$new(
  alldata, time, censor, covariates,
  quantiles, tolerance, epsilon, kappa=exp(eta), nsize,
  pureweight=1, mixedweight=1,Kconstant=1,Kpoly_order =2
)

#Create tree objects 
fullsubset<-1:n
dipolarqtree<-Dipolar.model$createtree(fullsubset)

#print full trees and log rank statistics  
print(dipolarqtree,"lrstat")

X = as.matrix(dipolarqtree$data)
mupmdiff = dipolarqtree$opt_w0_mupmdiff$mupmdiff
w0 = dipolarqtree$opt_w0_mupmdiff$w0

Kern.X = function(X,x) {
  (X%*%x+1)^2
}

#eta1
cc1 <- emdbook::curve3d(t(mupmdiff) %*%  Kern.X(X,c(x,y)) + w0,n = c(250,250), 
                       xlim=c(-.1,2.2), ylim=c(0,2), sys3d="none")
dimnames(cc1$z) <- list(cc1$x,cc1$y)
mm1 <- reshape2::melt(cc1$z)
#eta2
#cc2 <- emdbook::curve3d(t(mupmdiff) %*%  Kern.X(X,c(x,y)) + w0,n = c(250,250), 
#                      xlim=c(-.1,2.2), ylim=c(0,2), sys3d="none")
#dimnames(cc2$z) <- list(cc2$x,cc2$y)
#mm2 <- reshape2::melt(cc2$z)
#eta3
#cc3 <- emdbook::curve3d(t(mupmdiff) %*%  Kern.X(X,c(x,y)) + w0,n = c(250,250), 
#                       xlim=c(-.1,2.2), ylim=c(0,1.8), sys3d="none")
#dimnames(cc3$z) <- list(cc3$x,cc3$y)
#mm3 <- reshape2::melt(cc3$z)

g0<- ggplot() +   geom_point(data=alldatac,aes(x=x,y=x.1,size=stop,color=as.character(status)),show.legend = FALSE,color='green',alpha=.4,position="jitter")+
  geom_point(data=alldatas, aes(x = x, y = x.1, size = stop,color=stop),alpha=.4,position="jitter")+xlim(-.1,2.2)+ylim(.3,1.6)+
  scale_size_continuous(name = "Surv\nTime")+scale_color_continuous(name = "Surv\nTime")+labs(x="x1",y="x2")+
  stat_function(data=sim.data,aes(x=x),fun=function(x) (a+1)-2*a*x+a*x^2,colour="black")
g1<-g0 + geom_contour(data=mm1, aes(x=Var1,y=Var2,z=value),breaks=0,colour="red") #eta1
#geom_contour(data=mm2, aes(x=Var1,y=Var2,z=value),breaks=0,colour="deeppink4") + #eta2 
#geom_contour(data=mm3, aes(x=Var1,y=Var2,z=value),breaks=0,colour="red") # + #eta3
#ggsave("plot3.jpeg",g1,width=5,height=5)



#Elliptical Hazard, Quadratric Kernel ----
#NONLINEAR HAZARD - Ellipse
#Generate some data with parabolic hazard y = h(x) = ax^2-2ax+a+1
#Define parameters for Weibull distribution for censoring distribution
set.seed(123)
pc=1; b0c =-3.5;
anc.cens = pc; beta0.cens = -b0c/pc; 
#Define parameters for Weibull distribution for time-to-event distribution
pt = 1; a=3.5; b=2; c=3; b0t =c*(1/a^2+1/b^2-1);  betat = list(-2/a^2,-2/b^2,1/a^2,1/b^2)
anc.ev=pt; beta0.ev = -b0t/pt; beta =  lapply(betat,"*",-c/pt)
#Generate data 
sim.data <- simple.surv.sim.int(n=180, foltime=25, dist.ev=c('weibull'),
                                anc.ev,beta0.ev,dist.cens=c('weibull'),anc.cens,
                                beta0.cens,z=NULL, beta, x=list(c("normal", 1,1),c("normal", 1,2),c("quad",1),c("quad",2)))
prop_censored=1-sum(sim.data$status)/nrow(sim.data)
print(prop_censored)
#Create Hazard column 
sim.data['hr']=b0t + betat[[1]]*sim.data['x']+ betat[[2]]*sim.data['x.1']+betat[[3]]*sim.data['x.2']+betat[[4]]*sim.data['x.3']
#Correlation between Hazard and survival times
cor(sim.data['stop'],sim.data['hr'])

#Scatter Plot and Level Hazard with coded points 

#Fit a model 
#Fit a splitting surface
alldata = sim.data; covariates = c("x","x.1"); time = "stop"; censor = "status";
alldatas=alldata[alldata$status==1,]
alldatac=alldata[alldata$status==0,]

alldata = sim.data; covariates = c("x","x.1"); time = "stop"; censor = "status";

X = alldata[, covariates]
distX = c(dist(X))
epsilon = quantile(distX[distX != 0], probs = c(0.23))
quantiles = c(.35,.65); tolerance = 10^-2; eta=3;  
nsize = 10

n=nrow(alldata)

#Instantiate our model classes 
Dipolar.model <- DipolarSurvivalTree_PolyKernel$new(
  alldata, time, censor, covariates,
  quantiles, tolerance, epsilon, kappa=exp(eta), nsize,
  pureweight=1, mixedweight=1,Kconstant=1,Kpoly_order =2
)

#Create tree objects 
fullsubset<-1:n
dipolarqtree<-Dipolar.model$createtree(fullsubset)

#print full trees and log rank statistics  
print(dipolarqtree,"lrstat")

X = as.matrix(dipolarqtree$data)
mupmdiff = dipolarqtree$opt_w0_mupmdiff$mupmdiff
w0 = dipolarqtree$opt_w0_mupmdiff$w0

Kern.X = function(X,x) {
  (X%*%x+1)^2
}

#Plot data plus hazard   
g0<- ggplot() + 
  geom_point(data=alldatas, aes(x = x, y = x.1, size = stop, color=stop),alpha=0.4)+
  geom_point(data=alldatac, aes(x = x, y = x.1,size=stop),show.legend = FALSE,color='green',alpha=0.4)

cc <- emdbook::curve3d((x-1)^2/a^2+(y-1)^2/b^2,n = c(200,200),
                       xlim=c(-3,3), ylim=c(-5,5), sys3d="none")
dimnames(cc$z) <- list(cc$x,cc$y)
mm <- reshape2::melt(cc$z)

g1<-g0 + geom_contour(data=mm,aes(x=Var1,y=Var2,z=value),breaks=.32,colour="black") +
  scale_size_continuous(name = "Surv\nTime")+scale_color_continuous(name = "Surv\nTime")+
  labs(x="x1",y="x2")

cc1 <- emdbook::curve3d(t(mupmdiff) %*%  Kern.X(X,c(x,y)) + w0,n = c(250,250), 
                       xlim=c(-3,7), ylim=c(-2,6), sys3d="none")
dimnames(cc1$z) <- list(cc1$x,cc1$y)
mm1 <- reshape2::melt(cc1$z)

#cc2<- emdbook::curve3d(t(mupmdiff) %*%  Kern.X(X,c(x,y)) + w0,n = c(250,250), 
#                        xlim=c(-3,7), ylim=c(-2,6), sys3d="none")
#dimnames(cc2$z) <- list(cc2$x,cc2$y)
#mm2 <- reshape2::melt(cc2$z)

#cc3<- emdbook::curve3d(t(mupmdiff) %*%  Kern.X(X,c(x,y)) + w0,n = c(250,250), 
#                       xlim=c(-3,7), ylim=c(-2,6), sys3d="none")
#dimnames(cc3$z) <- list(cc3$x,cc3$y)
#mm3 <- reshape2::melt(cc3$z)

g2<-g1+ labs(x="x1",y="x2") + ylim(-3,4.5)+geom_contour(data=mm1, aes(x=Var1,y=Var2,z=value),breaks=0,colour='red') #eta1+
#geom_contour(data=mm2, aes(x=Var1,y=Var2,z=value),breaks=0,colour='deeppink4') #eta2+
#geom_contour(data=mm3, aes(x=Var1,y=Var2,z=value),breaks=0,colour="darkred") #eta3
#ggsave("plot2.jpeg",g0,width=5,height=5)

#Hyperbolic Hazard, Quadratic Kernel ----
#NONLINEAR HAZARD - Hyperbola
#Generate some data with parabolic hazard y = h(x) = ax^2-2ax+a+1
#Define parameters for Weibull distribution for censoring distribution
set.seed(123)
pc=1; b0c =-6.4;
anc.cens = pc; beta0.cens = -b0c/pc; 
#Define parameters for Weibull distribution for time-to-event distribution
pt = 1; a=3.5; b=2; c=3; b0t =c*(1/a^2-1/b^2-1);  betat = list(2/a^2,-2/b^2,-1/a^2,1/b^2)
anc.ev=pt; beta0.ev = -b0t/pt; beta =  lapply(betat,"*",-c/pt)
#Generate data 
sim.data <- simple.surv.sim.int(n=180, foltime=100, dist.ev=c('weibull'),
                                anc.ev,beta0.ev,dist.cens=c('weibull'),anc.cens,
                                beta0.cens,z=NULL, beta, x=list(c("normal", 1,1),c("normal", 1,2),c("quad",1),c("quad",2)))
prop_censored=1-sum(sim.data$status)/nrow(sim.data)
print(prop_censored)
#Create Hazard column 
sim.data['hr']=b0t + betat[[1]]*sim.data['x']+ betat[[2]]*sim.data['x.1']+betat[[3]]*sim.data['x.2']+betat[[4]]*sim.data['x.3']
#Correlation between Hazard and survival times
cor(sim.data['stop'],sim.data['hr'])

#Fit a model 
alldata = sim.data; covariates = c("x","x.1"); time = "stop"; censor = "status";
alldatas=alldata[alldata$status==1,]
alldatac=alldata[alldata$status==0,]

X = alldata[, covariates]
distX = c(dist(X))
epsilon = quantile(distX[distX != 0], probs = c(0.23))
quantiles = c(.35,.65); tolerance = 10^-2; eta=3;  
nsize = 10

n=nrow(alldata)

#Instantiate our model classes 
Dipolar.model <- DipolarSurvivalTree_PolyKernel$new(
  alldata, time, censor, covariates,
  quantiles, tolerance, epsilon, kappa=exp(eta), nsize,
  pureweight=1, mixedweight=1,Kconstant=1,Kpoly_order =2
)

#Create tree objects 
fullsubset<-1:n
dipolarqtree<-Dipolar.model$createtree(fullsubset)

#print full trees and log rank statistics  
print(dipolarqtree,"lrstat")

X = as.matrix(dipolarqtree$data)
mupmdiff = dipolarqtree$opt_w0_mupmdiff$mupmdiff
w0 = dipolarqtree$opt_w0_mupmdiff$w0

Kern.X = function(X,x) {
  (X%*%x+1)^2
}

#Scatter Plot and Level Hazard with coded points 
g0<- ggplot() + 
  geom_point(data=alldatas, aes(x = x, y = x.1, size = stop, color=stop),alpha=0.4)+
  geom_point(data=alldatac, aes(x = x, y = x.1,size=stop),show.legend = FALSE,color='green',alpha=0.4)

cc <- emdbook::curve3d(-(x-1)^2/a^2+(y-1)^2/b^2,n = c(200,200),
                       xlim=c(-2,4), ylim=c(-3,6), sys3d="none")
dimnames(cc$z) <- list(cc$x,cc$y)
mm <- reshape2::melt(cc$z)

g1<-g0 + geom_contour(data=mm,
                      aes(x=Var1,y=Var2,z=value),breaks=.25,
                      colour="black")
cc1 <- emdbook::curve3d(t(mupmdiff) %*%  Kern.X(X,c(x,y)) + w0,n = c(250,250), 
                        xlim=c(-2,4), ylim=c(-3,6), sys3d="none")
dimnames(cc1$z) <- list(cc1$x,cc1$y)
mm1 <- reshape2::melt(cc1$z)

#cc2<- emdbook::curve3d(t(mupmdiff) %*%  Kern.X(X,c(x,y)) + w0,n = c(250,250), 
#                       xlim=c(-2,4), ylim=c(-3,6), sys3d="none")
#dimnames(cc2$z) <- list(cc2$x,cc2$y)
#mm2 <- reshape2::melt(cc2$z)

#cc3<- emdbook::curve3d(t(mupmdiff) %*%  Kern.X(X,c(x,y)) + w0,n = c(250,250), 
#                       xlim=c(-2,4), ylim=c(-3,6), sys3d="none")
#dimnames(cc3$z) <- list(cc3$x,cc3$y)
#mm3 <- reshape2::melt(cc3$z)

g2<-g1 +  scale_size_continuous(name = "Surv\nTime")+scale_color_continuous(name = "Surv\nTime")+ labs(x="x1",y="x2") +
  ylim(-3,4.5)+
  geom_contour(data=mm1, aes(x=Var1,y=Var2,z=value),breaks=0,colour='red') #+
 # geom_contour(data=mm2, aes(x=Var1,y=Var2,z=value),breaks=.3,colour='deeppink4') +
#  geom_contour(data=mm3, aes(x=Var1,y=Var2,z=value),breaks=0,colour="darkred")
#ggsave("plot2.jpeg",g2,width=5,height=5)


#Parabolic Hazard, Gaussian Kernel ----
#NONLINEAR HAZARD - Downward turning parabola 
set.seed(123)
#Generate some data with parabolic hazard y = h(x) = ax^2-2ax+a+1
#Define parameters for Weibull distribution for censoring distribution
pc=1; b0c =-1.8;
anc.cens = pc; beta0.cens = -b0c/pc; 
#Define parameters for Weibull distribution for time-to-event distribution
pt = 1; a=-2;b=3.5; b0t =b*(-a-1); betat = list(b*2*a,b*1,b*(-a))
anc.ev=pt; beta0.ev = -b0t/pt; beta =  lapply(betat,"*",-1/pt)
#Generate data 
sim.data <- simple.surv.sim.int(n=180, foltime=10, dist.ev=c('weibull'),
                                anc.ev,beta0.ev,dist.cens=c('weibull'),anc.cens,
                                beta0.cens,z=NULL, beta, x=list(c("normal", 1,.4),c("normal", 1, .3),c("quad",1)))
prop_censored=1-sum(sim.data$status)/nrow(sim.data)
print(prop_censored)
#Create Hazard column 
sim.data['hr']=b0t + betat[[1]]*sim.data['x']+ betat[[2]]*sim.data['x.1']+betat[[3]]*sim.data['x.2']
#Correlation between Hazard and survival times
cor(sim.data['stop'],sim.data['hr'])

#Fit a splitting surface
alldata = sim.data; covariates = c("x","x.1"); time = "stop"; censor = "status";
alldatas=alldata[alldata$status==1,]
alldatac=alldata[alldata$status==0,]

X = alldata[, covariates]
distX = c(dist(X))
epsilon = quantile(distX[distX != 0], probs = c(0.20))
quantiles = c(.35,.65); tolerance = 10^-2; eta =3 
nsize = 10

n=nrow(alldata)

ksigma=intvarfun(alldata,covariates)
#Instantiate our model classes
Dipolar.model1 <- DipolarSurvivalTree_GaussKernel$new(
  alldata, time, censor, covariates,
  quantiles, tolerance, epsilon, kappa=exp(eta), nsize,
  pureweight=1, mixedweight=1, Ksigma = ksigma
)

#Create tree objects
fullsubset<-1:n
dipolarqtree1<-Dipolar.model1$createtree(fullsubset)

#print full trees and log rank statistics
print(dipolarqtree1,"lrstat")


NodeCurr = dipolarqtree1
X = as.matrix(NodeCurr$data)
mupmdiff = NodeCurr$opt_w0_mupmdiff$mupmdiff
w0 = NodeCurr$opt_w0_mupmdiff$w0

Kern.X = Dipolar.model1$K

g0<- ggplot() +   geom_point(data=alldatac,aes(x=x,y=x.1,size=stop,color=as.character(status)),show.legend = FALSE,color='green',alpha=.4,position="jitter")+
  geom_point(data=alldatas, aes(x = x, y = x.1, size = stop,color=stop),alpha=.4,position="jitter")+xlim(-.1,2.2)+ylim(.3,1.6)+
  scale_size_continuous(name = "Surv\nTime")+scale_color_continuous(name = "Surv\nTime")+labs(x="x1",y="x2")

g1<-g0+stat_function(data=sim.data,aes(x=x),fun=function(x) (a+1)-2*a*x+a*x^2,colour="black")

cc1 <- emdbook::curve3d(t(mupmdiff) %*%  Kern.X(X,c(x,y)) + w0,n = c(250,250), 
                        xlim=c(0,2.1), ylim=c(-2,4), sys3d="none")
dimnames(cc1$z) <- list(cc1$x,cc1$y)
mm1 <- reshape2::melt(cc1$z) #eta1

#cc2 <- emdbook::curve3d(t(mupmdiff) %*%  Kern.X(X,c(x,y)) + w0,n = c(250,250), 
#                        xlim=c(-2,3.2), ylim=c(-2,4), sys3d="none")
#dimnames(cc2$z) <- list(cc2$x,cc2$y)
#mm2 <- reshape2::melt(cc2$z) #eta2

#cc3 <- emdbook::curve3d(t(mupmdiff) %*%  Kern.X(X,c(x,y)) + w0,n = c(250,250), 
#                         xlim=c(-2,3.2), ylim=c(-2,4), sys3d="none")
# dimnames(cc3$z) <- list(cc3$x,cc3$y)
# mm3 <- reshape2::melt(cc3$z)

g2<-g1 + geom_contour(data=mm1, aes(x=Var1,y=Var2,z=value),breaks=0,colour="red")#eta1 +
#geom_contour(data=mm2, aes(x=Var1,y=Var2,z=value),breaks=0,colour='deeppink4')#eta2 +
#geom_contour(data=mm3, aes(x=Var1,y=Var2,z=value),breaks=0,colour="darkred") #eta3 +

 
#Ellipse Hazard, Gaussian Kernel----
set.seed(123)
pc=1; b0c =-4.5;
anc.cens = pc; beta0.cens = -b0c/pc; 
#Define parameters for Weibull distribution for time-to-event distribution
pt = 1; a=1; b=2; c=5.5; b0t =c*(1/a^2+1/b^2-1);  betat = list(-2/a^2,-2/b^2,1/a^2,1/b^2)
anc.ev=pt; beta0.ev = -b0t/pt; beta =  lapply(betat,"*",-c/pt)

#Generate data 
sim.data <- simple.surv.sim.int(n=180, foltime=80, dist.ev=c('weibull'),
                                anc.ev,beta0.ev,dist.cens=c('weibull'),anc.cens,
                                beta0.cens,z=NULL, beta, x=list(c("normal", 1,1),c("normal", 1,2),c("quad",1),c("quad",2)))
prop_censored=1-sum(sim.data$status)/nrow(sim.data)
print(prop_censored)
#Create Hazard column 
sim.data['hr']=b0t + betat[[1]]*sim.data['x']+ betat[[2]]*sim.data['x.1']+betat[[3]]*sim.data['x.2']+betat[[4]]*sim.data['x.3']
#Correlation between Hazard and survival times
cor(sim.data['stop'],sim.data['hr'])

#Scatter Plot and Level Hazard with coded points 
# ggplot() + 
#   geom_point(data=sim.data, aes(x = x, y = x.1, size = stop, color=stop,alpha=0.1))

#Fit a splitting surface
alldata = sim.data; covariates = c("x","x.1"); time = "stop"; censor = "status";
alldatas=alldata[alldata$status==1,]
alldatac=alldata[alldata$status==0,]

X = alldata[, covariates]
distX = c(dist(X))
epsilon = quantile(distX[distX != 0], probs = c(0.25))
quantiles = c(.35,.65); tolerance = 10^-2; eta =6
nsize = 20

n=nrow(alldata)

ksigma=intvarfun(alldata,covariates)
#Instantiate our model classes
Dipolar.model1 <- DipolarSurvivalTree_GaussKernel$new(
  alldata, time, censor, covariates,
  quantiles, tolerance, epsilon, kappa=exp(eta), nsize,
  pureweight=1, mixedweight=1, Ksigma = ksigma
)


#Create tree objects
fullsubset<-1:n
dipolarqtree1<-Dipolar.model1$createtree(fullsubset)

#print full trees and log rank statistics
print(dipolarqtree1,"lrstat")


NodeCurr = dipolarqtree1
X = as.matrix(NodeCurr$data)
mupmdiff = NodeCurr$opt_w0_mupmdiff$mupmdiff
w0 = NodeCurr$opt_w0_mupmdiff$w0

Kern.X = Dipolar.model1$K


g0<- ggplot() + 
  geom_point(data=alldatas, aes(x = x, y = x.1, size = stop, color=stop),alpha=0.4, position="jitter")+
  geom_point(data=alldatac, aes(x = x, y = x.1,size=stop),show.legend = FALSE,color='green',alpha=0.4,position="jitter")

cc <- emdbook::curve3d((x-1)^2/a^2+(y-1)^2/b^2,n = c(200,200),
                       xlim=c(-3,3), ylim=c(-5,5), sys3d="none")
dimnames(cc$z) <- list(cc$x,cc$y)
mm <- reshape2::melt(cc$z)

g1<-g0 + geom_contour(data=mm,aes(x=Var1,y=Var2,z=value),breaks=.70,colour="black") +
scale_color_continuous(name = "Surv\nTime")+  scale_size_continuous(name = "Surv\nTime")+
  labs(x="x1",y="x2")+ylim(-3,4.5) +xlim(-1,3)


cc1 <- emdbook::curve3d(t(mupmdiff) %*%  Kern.X(X,c(x,y)) + w0,n = c(250,250), 
                 xlim=c(0,2.1), ylim=c(-2,4), sys3d="none")
dimnames(cc1$z) <- list(cc1$x,cc1$y)
mm1 <- reshape2::melt(cc1$z) #eta1
 
#cc2 <- emdbook::curve3d(t(mupmdiff) %*%  Kern.X(X,c(x,y)) + w0,n = c(250,250), 
#                        xlim=c(-2,3.2), ylim=c(-2,4), sys3d="none")
#dimnames(cc2$z) <- list(cc2$x,cc2$y)
#mm2 <- reshape2::melt(cc2$z) #eta2
 
#cc3 <- emdbook::curve3d(t(mupmdiff) %*%  Kern.X(X,c(x,y)) + w0,n = c(250,250), 
#                         xlim=c(-2,3.2), ylim=c(-2,4), sys3d="none")
# dimnames(cc3$z) <- list(cc3$x,cc3$y)
# mm3 <- reshape2::melt(cc3$z)
 
 g2<-g1 + geom_contour(data=mm1, aes(x=Var1,y=Var2,z=value),breaks=0,colour="red")#eta1 +
   #geom_contour(data=mm2, aes(x=Var1,y=Var2,z=value),breaks=0,colour='deeppink4')#eta2 +
   #geom_contour(data=mm3, aes(x=Var1,y=Var2,z=value),breaks=0,colour="darkred") #eta3 +
  

 
 
 

 
 
 
 