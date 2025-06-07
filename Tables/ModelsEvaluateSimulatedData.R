#CLEAR, SET WORKSPACE, LOAD LIBRARIES ----
# Remove preexisting variables from environment
rm(list=ls())
# Clear console
cat("\014")
# Force-set directory to this file in RStudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# Load the needed survsim scripts  
source("../ExtraFiles/simple.surv.sim.int.R")
source("../ExtraFiles/simple.ev.sim.R")
source("../ExtraFiles/PruningScript.R")
source("../ExtraFiles/Predictionfunction.R")
source("../ExtraFiles/kappacrossvalidation.R")
source("../ExtraFiles/BrierScoreFun.R")
source("../Models/DipolarSurvivalTree_PolyKernel.R")
source("../Models/DipolarSurvivalTree_GaussKernel.R") 
source("../Models/UnivariateSurvivaltree.R")



#SIMULATION PLANAR (2vars) ----

pc=1; b0c =-3.3;
anc.cens = pc; beta0.cens = -b0c/pc; 
#Define parameters for Weibull distribution for time-to-event distribution
covariates = c("x","x.1"); time = "stop"; censor = "status";
pt = 1; b=1.2; b0t =-b*6; betat = list(-4*b,4*b)
anc.ev=pt; beta0.ev = -b0t/pt; beta =  lapply(betat,"*",-1/pt)

r1=1000;  #Number of loops 
Num_nodes_dt.l<-c(); Num_nodes_dtp.l<-c(); CIndex.l<-c(); BScores.l<-c();
Num_nodes_dt.q<-c(); Num_nodes_dtp.q<-c(); CIndex.q<-c(); BScores.q<-c();
Num_nodes_dt.g<-c(); Num_nodes_dtp.g<-c(); CIndex.g<-c(); BScores.g<-c();
Num_nodes_dt.u<-c(); Num_nodes_dtp.u<-c(); CIndex.u<-c(); BScores.u<-c();
for (j in 1:r1){
  ##Generate a training and testing set 
  sim.data <- simple.surv.sim.int(n=500, foltime=70, dist.ev=c('weibull'),
                                  anc.ev,beta0.ev,dist.cens=c('weibull'),anc.cens,
                                  beta0.cens,z=NULL, beta, x=list(c("normal", 1.2,.4),c("normal", 2.6, .3)))
  prop_censored=1-sum(sim.data$status)/nrow(sim.data);alldata=sim.data; 
  alldata.train=alldata[1:250,]; alldata.test=alldata[251:500,]
  X = alldata.train[, covariates]
  distX = c(dist(X))
  epsilon = quantile(distX[distX != 0], probs = c(0.23))
  quantiles = c(.25,.75); tolerance = 10^-2;  
  nsize = 20; n1=nrow(alldata.train);n2=nrow(alldata.test)

##Linear Model 
  #tune for eta linear model 
  etalist=seq(-4,4,by=0.5)
  etachoice=choosekappa(alldata.train,etalist,model="l")
  eta=etachoice[[3]]
  #Instantiate model 
  Dipolar.model.l <- DipolarSurvivalTree_PolyKernel$new(
    alldata.train, time, censor, covariates,
    quantiles, tolerance, epsilon, kappa=exp(eta), nsize,
    pureweight=1, mixedweight=1,Kconstant=0,Kpoly_order =1
  )
  #validate linear model
  dipolartree.l<-Dipolar.model.l$createtree(1:n1)
  Num_nodes_dt.l[j]<-dipolartree.l$totalCount
  dipolartree.prune.l<-bootstrapPruning(dipolartree.l,Dipolar.model.l,2.8)[[4]]
  Num_nodes_dtp.l[j]<-dipolartree.prune.l$totalCount
  IBSrange.l = gmsfun(alldata.test$stop,alldata.test$status,5)
  BScores.l[j]=Int_Brierscore(1:n2,dipolartree.prune.l,Dipolar.model.l,IBSrange.l,0)
  actualtime.l = alldata.test$stop
  predictedtime2.l = Dipolar.model.l$predicttime(alldata.test,dipolartree.prune.l)
  CIndex.l[j] = Dipolar.model.l$cindex(actualtime.l,predictedtime2.l,alldata.test$status)

##Quadratic Model 
#tune for eta quadratic model
  etalist=seq(-4,4,by=0.5)
  etachoice=choosekappa(alldata.train,etalist,model="q")
  eta=etachoice[[3]]
  #Instantiate quadratic model 
  Dipolar.model.q <- DipolarSurvivalTree_PolyKernel$new(
    alldata.train, time, censor, covariates,
    quantiles, tolerance, epsilon, kappa=exp(eta), nsize,
    pureweight=1, mixedweight=1,Kconstant=1,Kpoly_order =2
  )
  #validate quadratic model
  dipolartree.q<-Dipolar.model.q$createtree(1:n1)
  Num_nodes_dt.q[j]<-dipolartree.q$totalCount
  dipolartree.prune.q<-bootstrapPruning(dipolartree.q,Dipolar.model.q,2.8)[[4]]
  Num_nodes_dtp.q[j]<-dipolartree.prune.q$totalCount
  IBSrange.q = gmsfun(alldata.test$stop,alldata.test$status,5)
  BScores.q[j]=Int_Brierscore(1:n2,dipolartree.prune.q,Dipolar.model.q,IBSrange.q,0)
  actualtime.q = alldata.test$stop
  predictedtime2.q = Dipolar.model.q$predicttime(alldata.test,dipolartree.prune.q)
  CIndex.q[j] = Dipolar.model.q$cindex(actualtime.q,predictedtime2.q,alldata.test$status)

##Gaussian Model 
  #tune for eta Gaussian model 
  etalist=seq(-2,4,by=0.5)
  etachoice=choosekappa(alldata.train,etalist,model="g")
  eta=etachoice[[3]]
  #Instantiate Gaussian model 
  ksigma=intvarfun(alldata.train,covariates)
  Dipolar.model.g <- DipolarSurvivalTree_GaussKernel$new(
    alldata.train, time, censor, covariates,
    quantiles, tolerance, epsilon, kappa=exp(eta), nsize,
    pureweight=1, mixedweight=1, Ksigma = ksigma
  )
  #Validadate Gaussian Model 
  dipolartree.g<-Dipolar.model.g$createtree(1:n1)
  Num_nodes_dt.g[j]<-dipolartree.g$totalCount
  dipolartree.prune.g<-bootstrapPruning(dipolartree.g,Dipolar.model.g,2.8)[[4]]
  Num_nodes_dtp.g[j]<-dipolartree.prune.g$totalCount
  IBSrange.g = gmsfun(alldata.test$stop,alldata.test$status,5)
  BScores.g[j]=Int_Brierscore(1:n2,dipolartree.prune.g,Dipolar.model.g,IBSrange.g,0)
  actualtime.g = alldata.test$stop
  predictedtime2.g = Dipolar.model.g$predicttime(alldata.test,dipolartree.prune.g)
  CIndex.g[j] = Dipolar.model.g$cindex(actualtime.g,predictedtime2.g,alldata.test$status)

  #Instantiate model 
  Univariate.model <- UnivariateSurvivalTree$new(alldata.train, time, censor, covariates,nsize)  
  #validate univariate model
  univariatetree<-Univariate.model$createtree(1:n1)
  Num_nodes_dt.u[j]<-univariatetree$totalCount
  univariatetree.prune<-bootstrapPruning(univariatetree,Univariate.model,2.8)[[4]]
  Num_nodes_dtp.u[j]<-univariatetree.prune$totalCount
  IBSrange.u = gmsfun(alldata.test$stop,alldata.test$status,5)
  BScores.u[j]=Int_Brierscore(1:n2,univariatetree.prune,Univariate.model,IBSrange.u,1)
  actualtime.u = alldata.test$stop
  predictedtime2.u = Univariate.model$predicttime(alldata.test,univariatetree.prune)
  CIndex.u[j] = Univariate.model$cindex(actualtime.u,predictedtime2.u,alldata.test$status)
}
mean(BScores.l)
mean(CIndex.l)
mean(Num_nodes_dt.l)
mean(Num_nodes_dtp.l)


mean(BScores.q)
mean(CIndex.q)
mean(Num_nodes_dt.q)
mean(Num_nodes_dtp.q)


mean(BScores.g)
mean(CIndex.g)
mean(Num_nodes_dt.g)
mean(Num_nodes_dtp.g)


mean(BScores.u)
mean(CIndex.u)
mean(Num_nodes_dt.u)
mean(Num_nodes_dtp.u)


#SIMULATION PARABOLIC (2vars)----

#specify variable information 
#specify parameters for simulation
pc=1; b0c =-2;
anc.cens = pc; beta0.cens = -b0c/pc; 
covariates = c("x","x.1"); time = "stop"; censor = "status";
#Define parameters for Weibull distribution for time-to-event distribution
pt = 1; a=-2;b=3.6; b0t =b*(-a-1); betat = list(b*2*a,b*1,b*(-a))
anc.ev=pt; beta0.ev = -b0t/pt; beta =  lapply(betat,"*",-1/pt)

r2=1;  
Num_nodes_dt.l<-c(); Num_nodes_dtp.l<-c(); CIndex.l<-c(); BScores.l<-c();
Num_nodes_dt.q<-c(); Num_nodes_dtp.q<-c(); CIndex.q<-c(); BScores.q<-c();
Num_nodes_dt.g<-c(); Num_nodes_dtp.g<-c(); CIndex.g<-c(); BScores.g<-c();
Num_nodes_dt.u<-c(); Num_nodes_dtp.u<-c(); CIndex.u<-c(); BScores.u<-c();
for (j in 1:r2){
  ##Generate a training and test set  
  sim.data <- simple.surv.sim.int(n=500, foltime=15, dist.ev=c('weibull'),
                                  anc.ev,beta0.ev,dist.cens=c('weibull'),anc.cens,
                                  beta0.cens,z=NULL, beta, x=list(c("normal", 1,.4),c("normal", 1, .3),c("quad",1)))
  prop_censored=1-sum(sim.data$status)/nrow(sim.data);alldata=sim.data; 
  alldata.train=alldata[1:250,]; alldata.test=alldata[251:500,]
  X = alldata.train[, covariates]
  distX = c(dist(X))
  epsilon = quantile(distX[distX != 0], probs = c(0.23))
  quantiles = c(.25,.75); tolerance = 10^-2;   
  nsize = 20; n1=nrow(alldata.train);n2=nrow(alldata.test)
  
  ##Linear Model 
  #tune for eta linear model 
  etalist=seq(-4,4,by=0.5)
  etachoice=choosekappa(alldata.train,etalist,model="l")
  eta=etachoice[[3]]
  
  #Instantiate model 
  Dipolar.model.l <- DipolarSurvivalTree_PolyKernel$new(
    alldata.train, time, censor, covariates,
    quantiles, tolerance, epsilon, kappa=exp(eta), nsize,
    pureweight=1, mixedweight=1,Kconstant=0,Kpoly_order =1
  )
  #validate linear model
  dipolartree.l<-Dipolar.model.l$createtree(1:n1)
  Num_nodes_dt.l[j]<-dipolartree.l$totalCount
  dipolartree.prune.l<-bootstrapPruning(dipolartree.l,Dipolar.model.l,2.8)[[4]]
  Num_nodes_dtp.l[j]<-dipolartree.prune.l$totalCount
  IBSrange.l = gmsfun(alldata.test$stop,alldata.test$status,5)
  BScores.l[j]=Int_Brierscore(1:n2,dipolartree.prune.l,Dipolar.model.l,IBSrange.l,0)
  actualtime.l = alldata.test$stop
  predictedtime2.l = Dipolar.model.l$predicttime(alldata.test,dipolartree.prune.l)
  CIndex.l[j] = Dipolar.model.l$cindex(actualtime.l,predictedtime2.l,alldata.test$status)
  
  ##Quadratic Model 
  #tune for eta quadratic model
  etalist=seq(-4,4,by=0.5)
  etachoice=choosekappa(alldata.train,etalist,model="q")
  eta=etachoice[[3]]
  
  #Instantiate quadratic model 
  Dipolar.model.q <- DipolarSurvivalTree_PolyKernel$new(
    alldata.train, time, censor, covariates,
    quantiles, tolerance, epsilon, kappa=exp(eta), nsize,
    pureweight=1, mixedweight=1,Kconstant=1,Kpoly_order =2
  )
  #validate quadratic model
  dipolartree.q<-Dipolar.model.q$createtree(1:n1)
  Num_nodes_dt.q[j]<-dipolartree.q$totalCount
  dipolartree.prune.q<-bootstrapPruning(dipolartree.q,Dipolar.model.q,2.8)[[4]]
  Num_nodes_dtp.q[j]<-dipolartree.prune.q$totalCount
  IBSrange.q = gmsfun(alldata.test$stop,alldata.test$status,5)
  BScores.q[j]=Int_Brierscore(1:n2,dipolartree.prune.q,Dipolar.model.q,IBSrange.q,0)
  actualtime.q = alldata.test$stop
  predictedtime2.q = Dipolar.model.q$predicttime(alldata.test,dipolartree.prune.q)
  CIndex.q[j] = Dipolar.model.q$cindex(actualtime.q,predictedtime2.q,alldata.test$status)
  
  ##Gaussian Model 
  #tune for eta Gaussian model 
  etalist=seq(-2,4,by=0.5)
  etachoice=choosekappa(alldata.train,etalist,model="g")
  eta=etachoice[[3]]
  
  ksigma=intvarfun(alldata.train,covariates)
  #Instantiate Gaussian model 
  Dipolar.model.g <- DipolarSurvivalTree_GaussKernel$new(
    alldata.train, time, censor, covariates,
    quantiles, tolerance, epsilon, kappa=exp(eta), nsize,
    pureweight=1, mixedweight=1, Ksigma = ksigma
  )
  #Validadate Gaussian Model 
  dipolartree.g<-Dipolar.model.g$createtree(1:n1)
  Num_nodes_dt.g[j]<-dipolartree.g$totalCount
  dipolartree.prune.g<-bootstrapPruning(dipolartree.g,Dipolar.model.g,2.8)[[4]]
  Num_nodes_dtp.g[j]<-dipolartree.prune.g$totalCount
  IBSrange.g = gmsfun(alldata.test$stop,alldata.test$status,5)
  BScores.g[j]=Int_Brierscore(1:n2,dipolartree.prune.g,Dipolar.model.g,IBSrange.g,0)
  actualtime.g = alldata.test$stop
  predictedtime2.g = Dipolar.model.g$predicttime(alldata.test,dipolartree.prune.g)
  CIndex.g[j] = Dipolar.model.g$cindex(actualtime.g,predictedtime2.g,alldata.test$status)

  ##Univariate Model
  #Instantiate model 
  Univariate.model <- UnivariateSurvivalTree$new(alldata.train, time, censor, covariates,nsize)  
  #validate univariate model
  univariatetree<-Univariate.model$createtree(1:n1)
  Num_nodes_dt.u[j]<-univariatetree$totalCount
  univariatetree.prune<-bootstrapPruning(univariatetree,Univariate.model,2.8)[[4]]
  Num_nodes_dtp.u[j]<-univariatetree.prune$totalCount
  IBSrange.u = gmsfun(alldata.test$stop,alldata.test$status,5)
  BScores.u[j]=Int_Brierscore(1:n2,univariatetree.prune,Univariate.model,IBSrange.u,1)
  actualtime.u = alldata.test$stop
  predictedtime2.u = Univariate.model$predicttime(alldata.test,univariatetree.prune)
  CIndex.u[j] = Univariate.model$cindex(actualtime.u,predictedtime2.u,alldata.test$status)
  }

mean(BScores.l)
mean(CIndex.l)
mean(Num_nodes_dt.l)
mean(Num_nodes_dtp.l)

mean(BScores.q)
mean(CIndex.q)
mean(Num_nodes_dt.q)
mean(Num_nodes_dtp.q)

mean(BScores.g)
mean(CIndex.g)
mean(Num_nodes_dt.g)
mean(Num_nodes_dtp.g)

mean(BScores.u)
mean(CIndex.u)
mean(Num_nodes_dt.u)
mean(Num_nodes_dtp.u)



#SIMULATION ELLIPSE (2vars) ----

#specify variable information 
#specify parameters for simulation
covariates = c("x","x.1"); time = "stop"; censor = "status";
pc=1; b0c =-5.5;
anc.cens = pc; beta0.cens = -b0c/pc; 
#Define parameters for Weibull distribution for time-to-event distribution
pt = 1; a=3.5; b=2; c=1; b0t =c*(1/a^2+1/b^2-1);  betat = list(-2/a^2,-2/b^2,1/a^2,1/b^2)
anc.ev=pt; beta0.ev = -b0t/pt; beta =  lapply(betat,"*",-c/pt)

r1=1000;  
Num_nodes_dt.l<-c(); Num_nodes_dtp.l<-c(); CIndex.l<-c(); BScores.l<-c();
Num_nodes_dt.q<-c(); Num_nodes_dtp.q<-c(); CIndex.q<-c(); BScores.q<-c();
Num_nodes_dt.g<-c(); Num_nodes_dtp.g<-c(); CIndex.g<-c(); BScores.g<-c();
Num_nodes_dt.u<-c(); Num_nodes_dtp.u<-c(); CIndex.u<-c(); BScores.u<-c();
for (j in 1:r1){
  ##Generate a training set 
  sim.data <- simple.surv.sim.int(n=500, foltime=25, dist.ev=c('weibull'),
                                        anc.ev,beta0.ev,dist.cens=c('weibull'),anc.cens,
                                        beta0.cens,z=NULL, beta, x=list(c("normal", 1,1.5),c("normal", 1,2.5),c("quad",1),c("quad",2)))
  prop_censored=1-sum(sim.data$status)/nrow(sim.data);alldata=sim.data; 
  alldata.train=alldata[1:250,]; alldata.test=alldata[251:500,]
  X = alldata.train[, covariates]
  distX = c(dist(X))
  epsilon = quantile(distX[distX != 0], probs = c(0.23))
  quantiles = c(.25,.75); tolerance = 10^-2;  
  nsize = 20; n1=nrow(alldata.train);n2=nrow(alldata.test)
  
  ##Linear Model 
  etalist=seq(-4,4,by=0.5)
  etachoice=choosekappa(alldata.train,etalist,model="l")
  eta=etachoice[[3]]
  
  #Instantiate model 
  Dipolar.model.l <- DipolarSurvivalTree_PolyKernel$new(
    alldata.train, time, censor, covariates,
    quantiles, tolerance, epsilon, kappa=exp(eta), nsize,
    pureweight=1, mixedweight=1,Kconstant=0,Kpoly_order =1
  )
  #validate linear model
  dipolartree.l<-Dipolar.model.l$createtree(1:n1)
  Num_nodes_dt.l[j]<-dipolartree.l$totalCount
  dipolartree.prune.l<-bootstrapPruning(dipolartree.l,Dipolar.model.l,2.8)[[4]]
  Num_nodes_dtp.l[j]<-dipolartree.prune.l$totalCount
  IBSrange.l = gmsfun(alldata.test$stop,alldata.test$status,5)
  BScores.l[j]=Int_Brierscore(1:n2,dipolartree.prune.l,Dipolar.model.l,IBSrange.l,0)
  actualtime.l = alldata.test$stop
  predictedtime2.l = Dipolar.model.l$predicttime(alldata.test,dipolartree.prune.l)
  CIndex.l[j] = Dipolar.model.l$cindex(actualtime.l,predictedtime2.l,alldata.test$status)
  
  ##Quadratic Model 
  #tune for eta quadratic model
  etalist=seq(-4,4,by=0.5)
  etachoice=choosekappa(alldata.train,etalist,model="q")
  eta=etachoice[[3]]
  
  #Instantiate quadratic model 
  Dipolar.model.q <- DipolarSurvivalTree_PolyKernel$new(
    alldata.train, time, censor, covariates,
    quantiles, tolerance, epsilon, kappa=exp(eta), nsize,
    pureweight=1, mixedweight=1,Kconstant=1,Kpoly_order =2
  )
  #validate quadratic model
  dipolartree.q<-Dipolar.model.q$createtree(1:n1)
  Num_nodes_dt.q[j]<-dipolartree.q$totalCount
  dipolartree.prune.q<-bootstrapPruning(dipolartree.q,Dipolar.model.q,2.8)[[4]]
  Num_nodes_dtp.q[j]<-dipolartree.prune.q$totalCount
  IBSrange.q = gmsfun(alldata.test$stop,alldata.test$status,5)
  BScores.q[j]=Int_Brierscore(1:n2,dipolartree.prune.q,Dipolar.model.q,IBSrange.q,0)
  actualtime.q = alldata.test$stop
  predictedtime2.q = Dipolar.model.q$predicttime(alldata.test,dipolartree.prune.q)
  CIndex.q[j] = Dipolar.model.q$cindex(actualtime.q,predictedtime2.q,alldata.test$status)
  
  ##Gaussian Model 
  #tune for eta Gaussian model 
  etalist=seq(-2,4,by=0.5)
  etachoice=choosekappa(alldata.train,etalist,model="g")
  eta=etachoice[[3]]

  #Instantiate Gaussian model 
  ksigma=intvarfun(alldata.train,covariates)
  Dipolar.model.g <- DipolarSurvivalTree_GaussKernel$new(
    alldata.train, time, censor, covariates,
    quantiles, tolerance, epsilon, kappa=exp(eta), nsize,
    pureweight=1, mixedweight=1, Ksigma = ksigma
  )
  #Validadate Gaussian Model 
  dipolartree.g<-Dipolar.model.g$createtree(1:n1)
  Num_nodes_dt.g[j]<-dipolartree.g$totalCount
  dipolartree.prune.g<-bootstrapPruning(dipolartree.g,Dipolar.model.g,2.8)[[4]]
  Num_nodes_dtp.g[j]<-dipolartree.prune.g$totalCount
  IBSrange.g = gmsfun(alldata.test$stop,alldata.test$status,5)
  BScores.g[j]=Int_Brierscore(1:n2,dipolartree.prune.g,Dipolar.model.g,IBSrange.g,0)
  actualtime.g = alldata.test$stop
  predictedtime2.g = Dipolar.model.g$predicttime(alldata.test,dipolartree.prune.g)
  CIndex.g[j] = Dipolar.model.g$cindex(actualtime.g,predictedtime2.g,alldata.test$status)

  ##Univariate Model 
  #Instantiate model 
  Univariate.model <- UnivariateSurvivalTree$new(alldata.train, time, censor, covariates,nsize)  
  #validate univariate model
  univariatetree<-Univariate.model$createtree(1:n1)
  Num_nodes_dt.u[j]<-univariatetree$totalCount
  univariatetree.prune<-bootstrapPruning(univariatetree,Univariate.model,2.8)[[4]]
  Num_nodes_dtp.u[j]<-univariatetree.prune$totalCount
  IBSrange.u = gmsfun(alldata.test$stop,alldata.test$status,5)
  BScores.u[j]=Int_Brierscore(1:n2,univariatetree.prune,Univariate.model,IBSrange.u,1)
  actualtime.u = alldata.test$stop
  predictedtime2.u = Univariate.model$predicttime(alldata.test,univariatetree.prune)
  CIndex.u[j] = Univariate.model$cindex(actualtime.u,predictedtime2.u,alldata.test$status)
  }

mean(BScores.l)
mean(CIndex.l)
mean(Num_nodes_dt.l)
mean(Num_nodes_dtp.l)

mean(BScores.q)
mean(CIndex.q)
mean(Num_nodes_dt.q)
mean(Num_nodes_dtp.q)

mean(BScores.g)
mean(CIndex.g)
mean(Num_nodes_dt.g)
mean(Num_nodes_dtp.g)

mean(BScores.u)
mean(CIndex.u)
mean(Num_nodes_dt.u)
mean(Num_nodes_dtp.u)


#SIMULATION PLANAR Weibull (2vars) ----

#specify variable information 
#specify parameters for simulation
pc=1; b0c =-3.3;
anc.cens = pc; beta0.cens = -b0c/pc; 
#Define parameters for Weibull distribution for time-to-event distribution
covariates = c("x","x.1"); time = "stop"; censor = "status";
pt = 2; b=1.2; b0t =-b*6; betat = list(-4*b,4*b)
anc.ev=pt; beta0.ev = -b0t/pt; beta =  lapply(betat,"*",-1/pt)

r1=1000;  
Num_nodes_dt.l<-c(); Num_nodes_dtp.l<-c();  CIndex.l<-c(); BScores.l<-c();
Num_nodes_dt.q<-c(); Num_nodes_dtp.q<-c();  CIndex.q<-c(); BScores.q<-c();
Num_nodes_dt.g<-c(); Num_nodes_dtp.g<-c();  CIndex.g<-c(); BScores.g<-c();
for (j in 1:r1){
  #Generate a training set 
  sim.data <- simple.surv.sim.int(n=500, foltime=70, dist.ev=c('weibull'),
                                  anc.ev,beta0.ev,dist.cens=c('weibull'),anc.cens,
                                  beta0.cens,z=NULL, beta, x=list(c("normal", 1.2,.4),c("normal", 2.6, .3)))
  prop_censored=1-sum(sim.data$status)/nrow(sim.data);alldata=sim.data; 
  alldata.train=alldata[1:250,]; alldata.test=alldata[251:500,]
  X = alldata.train[, covariates]
  distX = c(dist(X))
  epsilon = quantile(distX[distX != 0], probs = c(0.23))
  quantiles = c(.25,.75); tolerance = 10^-2;  
  nsize = 20; n1=nrow(alldata.train);n2=nrow(alldata.test)
  
  #LINEAR MODEL
  #tune for eta linear model 
  etalist=seq(-4,4,by=0.5)
  etachoice=choosekappa(alldata.train,etalist,model="l")
  eta=etachoice[[3]]
  
  
  #Instantiate model 
  Dipolar.model.l <- DipolarSurvivalTree_PolyKernel$new(
    alldata.train, time, censor, covariates,
    quantiles, tolerance, epsilon, kappa=exp(eta), nsize,
    pureweight=1, mixedweight=1,Kconstant=0,Kpoly_order =1
  )
  #validate linear model
  dipolartree.l<-Dipolar.model.l$createtree(1:n1)
  Num_nodes_dt.l[j]<-dipolartree.l$totalCount
  dipolartree.prune.l<-bootstrapPruning(dipolartree.l,Dipolar.model.l,2.8)[[4]]
  Num_nodes_dtp.l[j]<-dipolartree.prune.l$totalCount
  IBSrange.l = gmsfun(alldata.test$stop,alldata.test$status,5)
  BScores.l[j]=Int_Brierscore(1:n2,dipolartree.prune.l,Dipolar.model.l,IBSrange.l,0)
  actualtime.l = alldata.test$stop
  predictedtime2.l = Dipolar.model.l$predicttime(alldata.test,dipolartree.prune.l)
  CIndex.l[j] = Dipolar.model.l$cindex(actualtime.l,predictedtime2.l,alldata.test$status)
  
  #QUADRATIC MODEL 
  #tune for eta quadratic model
  etalist=seq(-4,4,by=0.5)
  etachoice=choosekappa(alldata.train,etalist,model="q")
  eta=etachoice[[3]]
  
  
  #Instantiate quadratic model 
  Dipolar.model.q <- DipolarSurvivalTree_PolyKernel$new(
    alldata.train, time, censor, covariates,
    quantiles, tolerance, epsilon, kappa=exp(eta), nsize,
    pureweight=1, mixedweight=1,Kconstant=1,Kpoly_order =2
  )
  #validate quadratic model
  dipolartree.q<-Dipolar.model.q$createtree(1:n1)
  Num_nodes_dt.q[j]<-dipolartree.q$totalCount
  dipolartree.prune.q<-bootstrapPruning(dipolartree.q,Dipolar.model.q,2.8)[[4]]
  Num_nodes_dtp.q[j]<-dipolartree.prune.q$totalCount
  IBSrange.q = gmsfun(alldata.test$stop,alldata.test$status,5)
  BScores.q[j]=Int_Brierscore(1:n2,dipolartree.prune.q,Dipolar.model.q,IBSrange.q,0)
  actualtime.q = alldata.test$stop
  predictedtime2.q = Dipolar.model.q$predicttime(alldata.test,dipolartree.prune.q)
  CIndex.q[j] = Dipolar.model.q$cindex(actualtime.q,predictedtime2.q,alldata.test$status)
  
  #GAUSSIAN MODEL 
  #tune for eta Gaussian model 
  etalist=seq(-2,4,by=0.5)
  etachoice=choosekappa(alldata.train,etalist,model="g")
  eta=etachoice[[3]]
  
  
  #Instantiate Gaussian model 
  ksigma=intvarfun(alldata.train,covariates)
  Dipolar.model.g <- DipolarSurvivalTree_GaussKernel$new(
    alldata.train, time, censor, covariates,
    quantiles, tolerance, epsilon, kappa=exp(eta), nsize,
    pureweight=1, mixedweight=1, Ksigma =ksigma
  )
  #Validadate Gaussian Model 
  dipolartree.g<-Dipolar.model.g$createtree(1:n1)
  Num_nodes_dt.g[j]<-dipolartree.g$totalCount
  dipolartree.prune.g<-bootstrapPruning(dipolartree.g,Dipolar.model.g,2.8)[[4]]
  Num_nodes_dtp.g[j]<-dipolartree.prune.g$totalCount
  IBSrange.g = gmsfun(alldata.test$stop,alldata.test$status,5)
  BScores.g[j]=Int_Brierscore(1:n2,dipolartree.prune.g,Dipolar.model.g,IBSrange.g,0)
  actualtime.g = alldata.test$stop
  predictedtime2.g = Dipolar.model.g$predicttime(alldata.test,dipolartree.prune.g)
  CIndex.g[j] = Dipolar.model.g$cindex(actualtime.g,predictedtime2.g,alldata.test$status)
}

mean(BScores.l)
mean(CIndex.l)
mean(Num_nodes_dt.l)
mean(Num_nodes_dtp.l)

mean(BScores.q)
mean(CIndex.q)
mean(Num_nodes_dt.q)
mean(Num_nodes_dtp.q)

mean(BScores.g)
mean(CIndex.g)
mean(Num_nodes_dt.g)
mean(Num_nodes_dtp.g)


#SIMULATION ELLIPSE Weibull (2vars)----

#specify variable information 
#specify parameters for simulation
covariates = c("x","x.1"); time = "stop"; censor = "status";
pc=2; b0c =-5.5;
anc.cens = pc; beta0.cens = -b0c/pc; 
#Define parameters for Weibull distribution for time-to-event distribution
pt = 1; a=3.5; b=2; c=1; b0t =c*(1/a^2+1/b^2-1);  betat = list(-2/a^2,-2/b^2,1/a^2,1/b^2)
anc.ev=pt; beta0.ev = -b0t/pt; beta =  lapply(betat,"*",-c/pt)

r1=1000;  
Num_nodes_dt.l<-c(); Num_nodes_dtp.l<-c(); CIndex.l<-c(); BScores.l<-c();
Num_nodes_dt.q<-c(); Num_nodes_dtp.q<-c();  CIndex.q<-c(); BScores.q<-c();
Num_nodes_dt.g<-c(); Num_nodes_dtp.g<-c();  CIndex.g<-c(); BScores.g<-c();
for (j in 1:r1){
  #Generate a training set 
  sim.data <- simple.surv.sim.int(n=500, foltime=25, dist.ev=c('weibull'),
                                  anc.ev,beta0.ev,dist.cens=c('weibull'),anc.cens,
                                  beta0.cens,z=NULL, beta, x=list(c("normal", 1,1.5),c("normal", 1,2.5),c("quad",1),c("quad",2)))
  prop_censored=1-sum(sim.data$status)/nrow(sim.data);alldata=sim.data; 
  alldata.train=alldata[1:250,]; alldata.test=alldata[251:500,]
  X = alldata.train[, covariates]
  distX = c(dist(X))
  epsilon = quantile(distX[distX != 0], probs = c(0.23))
  quantiles = c(.25,.75); tolerance = 10^-2; eta=3;  
  nsize = 20; n1=nrow(alldata.train);n2=nrow(alldata.test)
  
  #LINEAR MODEL
  #tune for eta linear model 
  etalist=seq(-4,4,by=0.5)
  etachoice=choosekappa(alldata.train,etalist,model="l")
  eta=etachoice[[3]]
  
  #Instantiate model 
  Dipolar.model.l <- DipolarSurvivalTree_PolyKernel$new(
    alldata.train, time, censor, covariates,
    quantiles, tolerance, epsilon, kappa=exp(eta), nsize,
    pureweight=1, mixedweight=1,Kconstant=0,Kpoly_order =1
  )
  #validate linear model
  dipolartree.l<-Dipolar.model.l$createtree(1:n1)
  Num_nodes_dt.l[j]<-dipolartree.l$totalCount
  dipolartree.prune.l<-bootstrapPruning(dipolartree.l,Dipolar.model.l,2.8)[[4]]
  Num_nodes_dtp.l[j]<-dipolartree.prune.l$totalCount
  IBSrange.l = gmsfun(alldata.test$stop,alldata.test$status,5)
  BScores.l[j]=Int_Brierscore(1:n2,dipolartree.prune.l,Dipolar.model.l,IBSrange.l,0)
  actualtime.l = alldata.test$stop
  predictedtime2.l = Dipolar.model.l$predicttime(alldata.test,dipolartree.prune.l)
  CIndex.l[j] = Dipolar.model.l$cindex(actualtime.l,predictedtime2.l,alldata.test$status)
  
  #QUADRATIC MODEL 
  #tune for eta quadratic model
  etalist=seq(-4,4,by=0.5)
  etachoice=choosekappa(alldata.train,etalist,model="q")
  eta=etachoice[[3]]
  
  #Instantiate quadratic model 
  Dipolar.model.q <- DipolarSurvivalTree_PolyKernel$new(
    alldata.train, time, censor, covariates,
    quantiles, tolerance, epsilon, kappa=exp(eta), nsize,
    pureweight=1, mixedweight=1,Kconstant=1,Kpoly_order =2
  )
  #validate quadratic model
  dipolartree.q<-Dipolar.model.q$createtree(1:n1)
  Num_nodes_dt.q[j]<-dipolartree.q$totalCount
  dipolartree.prune.q<-bootstrapPruning(dipolartree.q,Dipolar.model.q,2.8)[[4]]
  Num_nodes_dtp.q[j]<-dipolartree.prune.q$totalCount
  IBSrange.q = gmsfun(alldata.test$stop,alldata.test$status,5)
  BScores.q[j]=Int_Brierscore(1:n2,dipolartree.prune.q,Dipolar.model.q,IBSrange.q,0)
  actualtime.q = alldata.test$stop
  predictedtime2.q = Dipolar.model.q$predicttime(alldata.test,dipolartree.prune.q)
  CIndex.q[j] = Dipolar.model.q$cindex(actualtime.q,predictedtime2.q,alldata.test$status)
  
  #GAUSSIAN MODEL 
  etalist=seq(-2,4,by=0.5)
  etachoice=choosekappa(alldata.train,etalist,model="g")
  eta=etachoice[[3]]
  
  #Instantiate Gaussian model 
  ksigma=intvarfun(alldata.train,covariates)
  Dipolar.model.g <- DipolarSurvivalTree_GaussKernel$new(
    alldata.train, time, censor, covariates,
    quantiles, tolerance, epsilon, kappa=exp(eta), nsize,
    pureweight=1, mixedweight=1, Ksigma = ksigma
  )
  #Validadate Gaussian Model 
  dipolartree.g<-Dipolar.model.g$createtree(1:n1)
  Num_nodes_dt.g[j]<-dipolartree.g$totalCount
  dipolartree.prune.g<-bootstrapPruning(dipolartree.g,Dipolar.model.g,2.8)[[4]]
  Num_nodes_dtp.g[j]<-dipolartree.prune.g$totalCount
  IBSrange.g = gmsfun(alldata.test$stop,alldata.test$status,5)
  BScores.g[j]=Int_Brierscore(1:n2,dipolartree.prune.g,Dipolar.model.g,IBSrange.g,0)
  actualtime.g = alldata.test$stop
  predictedtime2.g = Dipolar.model.g$predicttime(alldata.test,dipolartree.prune.g)
  CIndex.g[j] = Dipolar.model.g$cindex(actualtime.g,predictedtime2.g,alldata.test$status)
}

mean(BScores.l)
mean(CIndex.l)
mean(Num_nodes_dt.l)
mean(Num_nodes_dtp.l)

mean(BScores.q)
mean(CIndex.q)
mean(Num_nodes_dt.q)
mean(Num_nodes_dtp.q)


mean(BScores.g)
mean(CIndex.g)
mean(Num_nodes_dt.g)
mean(Num_nodes_dtp.g)

















#SIMULATION HYPERPLANE (4d) ----

#specify variable information 
#specify parameters for simulation
pc=1; b0c=-7.3;
anc.cens = pc; beta0.cens = -b0c/pc; 
#Define parameters for Weibull distribution for time-to-event distribution
covariates = c("x","x.1","x.2","x.3"); time = "stop"; censor = "status";
pt = 1; b=.3; b0t =-b; betat = list(-2*b,-2*b,-1*b,-1*b)
anc.ev=pt; beta0.ev = -b0t/pt; beta =  lapply(betat,"*",-1/pt)

r1=1000;  
Num_nodes_dt.l<-c(); Num_nodes_dtp.l<-c();  CIndex.l<-c(); BScores.l<-c();
Num_nodes_dt.q<-c(); Num_nodes_dtp.q<-c();  CIndex.q<-c(); BScores.q<-c();
Num_nodes_dt.g<-c(); Num_nodes_dtp.g<-c();  CIndex.g<-c(); BScores.g<-c();
Num_nodes_dt.u<-c(); Num_nodes_dtp.u<-c(); CIndex.u<-c(); BScores.u<-c();
for (j in 1:r1){
  #Generate a training set 
  sim.data <- simple.surv.sim.int(n=500, foltime=150, dist.ev=c('weibull'),
                                  anc.ev,beta0.ev,dist.cens=c('weibull'),anc.cens,
                                  beta0.cens,z=NULL, beta, x=list(c("normal", 1.2,1.2),c("normal", 2.6,1),c("normal", 1.2,1.41),c("normal", 1.2,1.41)))
  prop_censored=1-sum(sim.data$status)/nrow(sim.data);alldata=sim.data; 
  alldata.train=alldata[1:250,]; alldata.test=alldata[251:500,]
  X = alldata.train[, covariates]
  distX = c(dist(X))
  epsilon = quantile(distX[distX != 0], probs = c(0.23))
  quantiles = c(.25,.75); tolerance = 10^-2;  
  nsize = 20; n1=nrow(alldata.train);n2=nrow(alldata.test)
  
  #LINEAR MODEL
  #tune for eta linear model 
  etalist=seq(-4,4,by=0.5)
  etachoice=choosekappa(alldata.train,etalist,model="l")
  eta=etachoice[[3]]
  
  
  #Instantiate model 
  Dipolar.model.l <- DipolarSurvivalTree_PolyKernel$new(
    alldata.train, time, censor, covariates,
    quantiles, tolerance, epsilon, kappa=exp(eta), nsize,
    pureweight=1, mixedweight=1,Kconstant=0,Kpoly_order =1
  )
  #validate linear model
  dipolartree.l<-Dipolar.model.l$createtree(1:n1)
  Num_nodes_dt.l[j]<-dipolartree.l$totalCount
  dipolartree.prune.l<-bootstrapPruning(dipolartree.l,Dipolar.model.l,2.8)[[4]]
  Num_nodes_dtp.l[j]<-dipolartree.prune.l$totalCount
  IBSrange.l = gmsfun(alldata.test$stop,alldata.test$status,5)
  BScores.l[j]=Int_Brierscore(1:n2,dipolartree.prune.l,Dipolar.model.l,IBSrange.l,0)
  actualtime.l = alldata.test$stop
  predictedtime2.l = Dipolar.model.l$predicttime(alldata.test,dipolartree.l)
  CIndex.l[j] = Dipolar.model.l$cindex(actualtime.l,predictedtime2.l,alldata.test$status)
  
  #QUADRATIC MODEL 
  #tune for eta quadratic model
  etalist=seq(-4,4,by=0.5)
  etachoice=choosekappa(alldata.train,etalist,model="q")
  eta=etachoice[[3]]
  
  #Instantiate quadratic model 
  Dipolar.model.q <- DipolarSurvivalTree_PolyKernel$new(
    alldata.train, time, censor, covariates,
    quantiles, tolerance, epsilon, kappa=exp(eta), nsize,
    pureweight=1, mixedweight=1,Kconstant=1,Kpoly_order =2
  )
  #validate quadratic model
  dipolartree.q<-Dipolar.model.q$createtree(1:n1)
  Num_nodes_dt.q[j]<-dipolartree.q$totalCount
  dipolartree.prune.q<-bootstrapPruning(dipolartree.q,Dipolar.model.q,2.8)[[4]]
  Num_nodes_dtp.q[j]<-dipolartree.prune.q$totalCount
  IBSrange.q = gmsfun(alldata.test$stop,alldata.test$status,5)
  BScores.q[j]=Int_Brierscore(1:n2,dipolartree.prune.q,Dipolar.model.q,IBSrange.q,0)
  actualtime.q = alldata.test$stop
  predictedtime2.q = Dipolar.model.q$predicttime(alldata.test,dipolartree.prune.q)
  CIndex.q[j] = Dipolar.model.q$cindex(actualtime.q,predictedtime2.q,alldata.test$status)
  
  #GAUSSIAN MODEL 
  #tune for eta Gaussian model 
  etalist=seq(-2,4,by=0.5)
  etachoice=choosekappa(alldata.train,etalist,model="g")
  eta=etachoice[[3]]
  
  ksigma=intvarfun(alldata.train,covariates)
  #Instantiate Gaussian model 
  Dipolar.model.g <- DipolarSurvivalTree_GaussKernel$new(
    alldata.train, time, censor, covariates,
    quantiles, tolerance, epsilon, kappa=exp(eta), nsize,
    pureweight=1, mixedweight=1, Ksigma = ksigma
  )
  #Validadate Gaussian Model 
  dipolartree.g<-Dipolar.model.g$createtree(1:n1)
  Num_nodes_dt.g[j]<-dipolartree.g$totalCount
  dipolartree.prune.g<-bootstrapPruning(dipolartree.g,Dipolar.model.g,2.8)[[4]]
  Num_nodes_dtp.g[j]<-dipolartree.prune.g$totalCount
  IBSrange.g = gmsfun(alldata.test$stop,alldata.test$status,5)
  BScores.g[j]=Int_Brierscore(1:n2,dipolartree.prune.g,Dipolar.model.g,IBSrange.g,0)
  actualtime.g = alldata.test$stop
  predictedtime2.g = Dipolar.model.g$predicttime(alldata.test,dipolartree.prune.g)
  CIndex.g[j] = Dipolar.model.g$cindex(actualtime.g,predictedtime2.g,alldata.test$status)

  ##Univariate Model 
  #Instantiate model 
  Univariate.model <- UnivariateSurvivalTree$new(alldata.train, time, censor, covariates,nsize)  
  #validate univariate model
  univariatetree<-Univariate.model$createtree(1:n1)
  Num_nodes_dt.u[j]<-univariatetree$totalCount
  univariatetree.prune<-bootstrapPruning(univariatetree,Univariate.model,2.8)[[4]]
  Num_nodes_dtp.u[j]<-univariatetree.prune$totalCount
  IBSrange.u = gmsfun(alldata.test$stop,alldata.test$status,5)
  BScores.u[j]=Int_Brierscore(1:n2,univariatetree.prune,Univariate.model,IBSrange.u,1)
  actualtime.u = alldata.test$stop
  predictedtime2.u = Univariate.model$predicttime(alldata.test,univariatetree.prune)
  CIndex.u[j] = Univariate.model$cindex(actualtime.u,predictedtime2.u,alldata.test$status)
}

mean(BScores.l)
mean(CIndex.l)
mean(Num_nodes_dt.l)
mean(Num_nodes_dtp.l)

mean(BScores.q)
mean(CIndex.q)
mean(Num_nodes_dt.q)
mean(Num_nodes_dtp.q)

mean(BScores.g)
mean(CIndex.g)
mean(Num_nodes_dt.g)
mean(Num_nodes_dtp.g)

mean(BScores.u)
mean(CIndex.u)
mean(Num_nodes_dt.u)
mean(Num_nodes_dtp.u)


#SIMULATION HYPERELLIPSE (4d)----

#specify variable information 
#specify parameters for simulation
pc=1; b0c=-6.3;
anc.cens = pc; beta0.cens = -b0c/pc; 
#Define parameters for Weibull distribution for time-to-event distribution
covariates = c("x","x.1","x.2"); time = "stop"; censor = "status";
pt = 1; b=.3; b0t =-b; betat = list(0,0,0,0,-b,-b,-b,-b)
anc.ev=pt; beta0.ev = -b0t/pt; beta =  lapply(betat,"*",-1/pt)

r1=1000;  
Num_nodes_dt.l<-c(); Num_nodes_dtp.l<-c();  CIndex.l<-c(); BScores.l<-c();
Num_nodes_dt.q<-c(); Num_nodes_dtp.q<-c();  CIndex.q<-c(); BScores.q<-c();
Num_nodes_dt.g<-c(); Num_nodes_dtp.g<-c();  CIndex.g<-c(); BScores.g<-c();
for (j in 1:r1){
  #Generate a training set 
  sim.data <- simple.surv.sim.int(n=500, foltime=150, dist.ev=c('weibull'),
                                  anc.ev,beta0.ev,dist.cens=c('weibull'),anc.cens,
                                  beta0.cens,z=NULL, beta, x=list(c("normal", 0.5,1.2),c("normal", 1,1),c("normal", -.5,1.41),c("normal", -.5,1.41),
                                                                  c("quad",1),c("quad",2),c("quad",3),c("quad",4)))
  prop_censored=1-sum(sim.data$status)/nrow(sim.data);alldata=sim.data; 
  alldata.train=alldata[1:250,]; alldata.test=alldata[251:500,]
  X = alldata.train[, covariates]
  distX = c(dist(X))
  epsilon = quantile(distX[distX != 0], probs = c(0.23))
  quantiles = c(.25,.75); tolerance = 10^-2;  
  nsize = 20; n1=nrow(alldata.train);n2=nrow(alldata.test)
  
  #LINEAR MODEL
  #tune for eta linear model 
  etalist=seq(-4,4,by=0.5)
  etachoice=choosekappa(alldata.train,etalist,model="l")
  eta=etachoice[[3]]
  
  #Instantiate model 
  Dipolar.model.l <- DipolarSurvivalTree_PolyKernel$new(
    alldata.train, time, censor, covariates,
    quantiles, tolerance, epsilon, kappa=exp(eta), nsize,
    pureweight=1, mixedweight=1,Kconstant=0,Kpoly_order =1
  )
  #validate linear model
  dipolartree.l<-Dipolar.model.l$createtree(1:n1)
  Num_nodes_dt.l[j]<-dipolartree.l$totalCount
  dipolartree.prune.l<-bootstrapPruning(dipolartree.l,Dipolar.model.l,2.8)[[4]]
  Num_nodes_dtp.l[j]<-dipolartree.prune.l$totalCount
  IBSrange.l = gmsfun(alldata.test$stop,alldata.test$status,5)
  BScores.l[j]=Int_Brierscore(1:n2,dipolartree.prune.l,Dipolar.model.l,IBSrange.l,0)
  actualtime.l = alldata.test$stop
  predictedtime2.l = Dipolar.model.l$predicttime(alldata.test,dipolartree.l)
  CIndex.l[j] = Dipolar.model.l$cindex(actualtime.l,predictedtime2.l,alldata.test$status)
  
  #QUADRATIC MODEL 
  #tune for eta quadratic model
  etalist=seq(-4,4,by=0.5)
  etachoice=choosekappa(alldata.train,etalist,model="q")
  eta=etachoice[[3]]
  
  #Instantiate quadratic model 
  Dipolar.model.q <- DipolarSurvivalTree_PolyKernel$new(
    alldata.train, time, censor, covariates,
    quantiles, tolerance, epsilon, kappa=exp(eta), nsize,
    pureweight=1, mixedweight=1,Kconstant=1,Kpoly_order =2
  )
  #validate quadratic model
  dipolartree.q<-Dipolar.model.q$createtree(1:n1)
  Num_nodes_dt.q[j]<-dipolartree.q$totalCount
  dipolartree.prune.q<-bootstrapPruning(dipolartree.q,Dipolar.model.q,2.8)[[4]]
  Num_nodes_dtp.q[j]<-dipolartree.prune.q$totalCount
  IBSrange.q = gmsfun(alldata.test$stop,alldata.test$status,5)
  BScores.q[j]=Int_Brierscore(1:n2,dipolartree.prune.q,Dipolar.model.q,IBSrange.q,0)
  actualtime.q = alldata.test$stop
  predictedtime2.q = Dipolar.model.q$predicttime(alldata.test,dipolartree.prune.q)
  CIndex.q[j] = Dipolar.model.q$cindex(actualtime.q,predictedtime2.q,alldata.test$status)
  
  #GAUSSIAN MODEL 
  #tune for eta Gaussian model 
  etalist=seq(-2,4,by=0.5)
  etachoice=choosekappa(alldata.train,etalist,model="g")
  eta=etachoice[[3]]
  
  #Instantiate Gaussian model 
  ksigma=intvarfun(alldata.train,covariates)
  Dipolar.model.g <- DipolarSurvivalTree_GaussKernel$new(
    alldata.train, time, censor, covariates,
    quantiles, tolerance, epsilon, kappa=exp(eta), nsize,
    pureweight=1, mixedweight=1, Ksigma = ksigma
  )
  #Validadate Gaussian Model 
  dipolartree.g<-Dipolar.model.g$createtree(1:n1)
  Num_nodes_dt.g[j]<-dipolartree.g$totalCount
  dipolartree.prune.g<-bootstrapPruning(dipolartree.g,Dipolar.model.g,2.8)[[4]]
  Num_nodes_dtp.g[j]<-dipolartree.prune.g$totalCount
  IBSrange.g = gmsfun(alldata.test$stop,alldata.test$status,5)
  BScores.g[j]=Int_Brierscore(1:n2,dipolartree.prune.g,Dipolar.model.g,IBSrange.g,0)
  actualtime.g = alldata.test$stop
  predictedtime2.g = Dipolar.model.g$predicttime(alldata.test,dipolartree.prune.g)
  CIndex.g[j] = Dipolar.model.g$cindex(actualtime.g,predictedtime2.g,alldata.test$status)
}

mean(BScores.l)
mean(CIndex.l)
mean(Num_nodes_dt.l)
mean(Num_nodes_dtp.l)

mean(BScores.q)
mean(CIndex.q)
mean(Num_nodes_dt.q)
mean(Num_nodes_dtp.q)

mean(BScores.g)
mean(CIndex.g)
mean(Num_nodes_dt.g)
mean(Num_nodes_dtp.g)

#SIMULATION HYPERPLANE (7d) ----

#specify variable information 
#specify parameters for simulation
pc=1; b0c=-7.3;
anc.cens = pc; beta0.cens = -b0c/pc; 
#Define parameters for Weibull distribution for time-to-event distribution
covariates = c("x","x.1","x.2","x.3","x.4","x.5","x.6"); time = "stop"; censor = "status";
pt = 1; b=.2; b0t =-b; betat = list(-2*b,-2*b,-1*b,-1*b,-1.5*b,-1.5*b,-2.5*b)
anc.ev=pt; beta0.ev = -b0t/pt; beta =  lapply(betat,"*",-1/pt)

r1=1000;  
Num_nodes_dt.l<-c(); Num_nodes_dtp.l<-c();  CIndex.l<-c(); BScores.l<-c();
Num_nodes_dt.q<-c(); Num_nodes_dtp.q<-c();  CIndex.q<-c(); BScores.q<-c();
Num_nodes_dt.g<-c(); Num_nodes_dtp.g<-c();  CIndex.g<-c(); BScores.g<-c();
Num_nodes_dt.u<-c(); Num_nodes_dtp.u<-c();  CIndex.u<-c(); BScores.u<-c();
for (j in 1:r1){
  #Generate a training set 
  sim.data <- simple.surv.sim.int(n=500, foltime=150, dist.ev=c('weibull'),
                                  anc.ev,beta0.ev,dist.cens=c('weibull'),anc.cens,
                                  beta0.cens,z=NULL, beta, x=list(c("normal", 1.2,1.2),c("normal", 2.6,1),c("normal", 1.2,1.41),c("normal", 1.2,1.41),c("normal", 1.2,1.41),c("normal", 1.2,1.41),c("normal", 1.2,1.41)))
  prop_censored=1-sum(sim.data$status)/nrow(sim.data);alldata=sim.data; 
  alldata.train=alldata[1:250,]; alldata.test=alldata[251:500,]
  X = alldata.train[, covariates]
  distX = c(dist(X))
  epsilon = quantile(distX[distX != 0], probs = c(0.23))
  quantiles = c(.25,.75); tolerance = 10^-2;  
  nsize = 20; n1=nrow(alldata.train);n2=nrow(alldata.test)
  
  #LINEAR MODEL
  #tune for eta linear model 
  etalist=seq(-4,4,by=0.5)
  etachoice=choosekappa(alldata.train,etalist,model="l")
  eta=etachoice[[3]]
  
  #Instantiate model 
  Dipolar.model.l <- DipolarSurvivalTree_PolyKernel$new(
    alldata.train, time, censor, covariates,
    quantiles, tolerance, epsilon, kappa=exp(eta), nsize,
    pureweight=1, mixedweight=1,Kconstant=0,Kpoly_order =1
  )
  #validate linear model
  dipolartree.l<-Dipolar.model.l$createtree(1:n1)
  Num_nodes_dt.l[j]<-dipolartree.l$totalCount
  dipolartree.prune.l<-bootstrapPruning(dipolartree.l,Dipolar.model.l,2.8)[[4]]
  Num_nodes_dtp.l[j]<-dipolartree.prune.l$totalCount
  IBSrange.l = gmsfun(alldata.test$stop,alldata.test$status,5)
  BScores.l[j]=Int_Brierscore(1:n2,dipolartree.prune.l,Dipolar.model.l,IBSrange.l,0)
  actualtime.l = alldata.test$stop
  predictedtime2.l = Dipolar.model.l$predicttime(alldata.test,dipolartree.l)
  CIndex.l[j] = Dipolar.model.l$cindex(actualtime.l,predictedtime2.l,alldata.test$status)
  
  #QUADRATIC MODEL 
  #tune for eta quadratic model
  etalist=seq(-4,4,by=0.5)
  etachoice=choosekappa(alldata.train,etalist,model="q")
  eta=etachoice[[3]]
  
  
  #Instantiate quadratic model 
  Dipolar.model.q <- DipolarSurvivalTree_PolyKernel$new(
    alldata.train, time, censor, covariates,
    quantiles, tolerance, epsilon, kappa=exp(eta), nsize,
    pureweight=1, mixedweight=1,Kconstant=1,Kpoly_order =2
  )
  #validate quadratic model
  dipolartree.q<-Dipolar.model.q$createtree(1:n1)
  Num_nodes_dt.q[j]<-dipolartree.q$totalCount
  dipolartree.prune.q<-bootstrapPruning(dipolartree.q,Dipolar.model.q,2.8)[[4]]
  Num_nodes_dtp.q[j]<-dipolartree.prune.q$totalCount
  IBSrange.q = gmsfun(alldata.test$stop,alldata.test$status,5)
  BScores.q[j]=Int_Brierscore(1:n2,dipolartree.prune.q,Dipolar.model.q,IBSrange.q,0)
  actualtime.q = alldata.test$stop
  predictedtime2.q = Dipolar.model.q$predicttime(alldata.test,dipolartree.prune.q)
  CIndex.q[j] = Dipolar.model.q$cindex(actualtime.q,predictedtime2.q,alldata.test$status)
  
  #GAUSSIAN MODEL 
  #tune for eta Gaussian model 
  etalist=seq(-2,4,by=0.5)
  etachoice=choosekappa(alldata.train,etalist,model="g")
  eta=etachoice[[3]]
  
  
  ksigma=intvarfun(alldata.train,covariates)
  #Instantiate Gaussian model 
  Dipolar.model.g <- DipolarSurvivalTree_GaussKernel$new(
    alldata.train, time, censor, covariates,
    quantiles, tolerance, epsilon, kappa=exp(eta), nsize,
    pureweight=1, mixedweight=1, Ksigma = ksigma
  )
  #Validadate Gaussian Model 
  dipolartree.g<-Dipolar.model.g$createtree(1:n1)
  Num_nodes_dt.g[j]<-dipolartree.g$totalCount
  dipolartree.prune.g<-bootstrapPruning(dipolartree.g,Dipolar.model.g,2.8)[[4]]
  Num_nodes_dtp.g[j]<-dipolartree.prune.g$totalCount
  IBSrange.g = gmsfun(alldata.test$stop,alldata.test$status,5)
  BScores.g[j]=Int_Brierscore(1:n2,dipolartree.prune.g,Dipolar.model.g,IBSrange.g,0)
  actualtime.g = alldata.test$stop
  predictedtime2.g = Dipolar.model.g$predicttime(alldata.test,dipolartree.prune.g)
  CIndex.g[j] = Dipolar.model.g$cindex(actualtime.g,predictedtime2.g,alldata.test$status)
  
  ##Univariate Model
  #Instantiate model 
  Univariate.model <- UnivariateSurvivalTree$new(alldata.train, time, censor, covariates,nsize)  
  #validate univariate model
  univariatetree<-Univariate.model$createtree(1:n1)
  Num_nodes_dt.u[j]<-univariatetree$totalCount
  univariatetree.prune<-bootstrapPruning(univariatetree,Univariate.model,2.8)[[4]]
  Num_nodes_dtp.u[j]<-univariatetree.prune$totalCount
  IBSrange.u = gmsfun(alldata.test$stop,alldata.test$status,5)
  BScores.u[j]=Int_Brierscore(1:n2,univariatetree.prune,Univariate.model,IBSrange.u,1)
  actualtime.u = alldata.test$stop
  predictedtime2.u = Univariate.model$predicttime(alldata.test,univariatetree.prune)
  CIndex.u[j] = Univariate.model$cindex(actualtime.u,predictedtime2.u,alldata.test$status)
}

mean(BScores.l)
mean(CIndex.l)
mean(Num_nodes_dt.l)
mean(Num_nodes_dtp.l)

mean(BScores.q)
mean(CIndex.q)
mean(Num_nodes_dt.q)
mean(Num_nodes_dtp.q)

mean(BScores.g)
mean(CIndex.g)
mean(Num_nodes_dt.g)
mean(Num_nodes_dtp.g)

mean(BScores.u)
mean(CIndex.u)
mean(Num_nodes_dt.u)
mean(Num_nodes_dtp.u)


#SIMULATION HYPERELLIPSE (7d) ----

#specify variable information 
#specify parameters for simulation
pc=1; b0c=-7.3;
anc.cens = pc; beta0.cens = -b0c/pc; 
#Define parameters for Weibull distribution for time-to-event distribution
covariates = c("x","x.1","x.2","x.3"); time = "stop"; censor = "status";
pt = 1; b=.3; b0t =-b; betat = list(0,0,0,0,0,0,0,-b,-b,-b,-b,-b,-b,-b)
anc.ev=pt; beta0.ev = -b0t/pt; beta =  lapply(betat,"*",-1/pt)

r1=1000;  
Num_nodes_dt.l<-c(); Num_nodes_dtp.l<-c();  CIndex.l<-c(); BScores.l<-c();
Num_nodes_dt.q<-c(); Num_nodes_dtp.q<-c();  CIndex.q<-c(); BScores.q<-c();
Num_nodes_dt.g<-c(); Num_nodes_dtp.g<-c(); CIndex.g<-c(); BScores.g<-c();
for (j in 1:r1){
  #Generate a training set 
  sim.data <- simple.surv.sim.int(n=500, foltime=150, dist.ev=c('weibull'),
                                  anc.ev,beta0.ev,dist.cens=c('weibull'),anc.cens,
                                  beta0.cens,z=NULL, beta, 
                                  x=list(c("normal", .2,1.2),c("normal", .6,1),c("normal", .2,1.41),c("normal", -.2,1.41),
                                         c("normal", .2,1.2),c("normal", .2,1.2),c("normal", -.2,1.2),c("quad",1),c("quad",2),
                                         c("quad",3),c("quad",4),c("quad",5),c("quad",6),c("quad",7)))
                                         
  prop_censored=1-sum(sim.data$status)/nrow(sim.data);alldata=sim.data; 
  alldata.train=alldata[1:250,]; alldata.test=alldata[251:500,]
  X = alldata.train[, covariates]
  distX = c(dist(X))
  epsilon = quantile(distX[distX != 0], probs = c(0.23))
  quantiles = c(.25,.75); tolerance = 10^-2;  
  nsize = 20; n1=nrow(alldata.train);n2=nrow(alldata.test)
  
  #LINEAR MODEL
  #tune for eta linear model 
  etalist=seq(-4,4,by=0.5)
  etachoice=choosekappa(alldata.train,etalist,model="l")
  eta=etachoice[[3]]
  
  #Instantiate model 
  Dipolar.model.l <- DipolarSurvivalTree_PolyKernel$new(
    alldata.train, time, censor, covariates,
    quantiles, tolerance, epsilon, kappa=exp(eta), nsize,
    pureweight=1, mixedweight=1,Kconstant=0,Kpoly_order =1
  )
  #validate linear model
  dipolartree.l<-Dipolar.model.l$createtree(1:n1)
  Num_nodes_dt.l[j]<-dipolartree.l$totalCount
  dipolartree.prune.l<-bootstrapPruning(dipolartree.l,Dipolar.model.l,2.8)[[4]]
  Num_nodes_dtp.l[j]<-dipolartree.prune.l$totalCount
  IBSrange.l = gmsfun(alldata.test$stop,alldata.test$status,5)
  BScores.l[j]=Int_Brierscore(1:n2,dipolartree.prune.l,Dipolar.model.l,IBSrange.l,0)
  actualtime.l = alldata.test$stop
  predictedtime2.l = Dipolar.model.l$predicttime(alldata.test,dipolartree.l)
  CIndex.l[j] = Dipolar.model.l$cindex(actualtime.l,predictedtime2.l,alldata.test$status)
  
  #QUADRATIC MODEL 
  #tune for eta quadratic model
  etalist=seq(-4,4,by=0.5)
  etachoice=choosekappa(alldata.train,etalist,model="q")
  eta=etachoice[[3]]
  
  #Instantiate quadratic model 
  Dipolar.model.q <- DipolarSurvivalTree_PolyKernel$new(
    alldata.train, time, censor, covariates,
    quantiles, tolerance, epsilon, kappa=exp(eta), nsize,
    pureweight=1, mixedweight=1,Kconstant=1,Kpoly_order =2
  )
  #validate quadratic model
  dipolartree.q<-Dipolar.model.q$createtree(1:n1)
  Num_nodes_dt.q[j]<-dipolartree.q$totalCount
  dipolartree.prune.q<-bootstrapPruning(dipolartree.q,Dipolar.model.q,2.8)[[4]]
  Num_nodes_dtp.q[j]<-dipolartree.prune.q$totalCount
  IBSrange.q = gmsfun(alldata.test$stop,alldata.test$status,5)
  BScores.q[j]=Int_Brierscore(1:n2,dipolartree.prune.q,Dipolar.model.q,IBSrange.q,0)
  actualtime.q = alldata.test$stop
  predictedtime2.q = Dipolar.model.q$predicttime(alldata.test,dipolartree.prune.q)
  CIndex.q[j] = Dipolar.model.q$cindex(actualtime.q,predictedtime2.q,alldata.test$status)
  
  #GAUSSIAN MODEL 
  #tune for eta Gaussian model 
  etalist=seq(-2,4,by=0.5)
  etachoice=choosekappa(alldata.train,etalist,model="g")
  eta=etachoice[[3]]
  
  #Instantiate Gaussian model 
  ksigma=intvarfun(alldata.train,covariates)
  Dipolar.model.g <- DipolarSurvivalTree_GaussKernel$new(
    alldata.train, time, censor, covariates,
    quantiles, tolerance, epsilon, kappa=exp(eta), nsize,
    pureweight=1, mixedweight=1, Ksigma = ksigma
  )
  #Validadate Gaussian Model 
  dipolartree.g<-Dipolar.model.g$createtree(1:n1)
  Num_nodes_dt.g[j]<-dipolartree.g$totalCount
  dipolartree.prune.g<-bootstrapPruning(dipolartree.g,Dipolar.model.g,2.8)[[4]]
  Num_nodes_dtp.g[j]<-dipolartree.prune.g$totalCount
  IBSrange.g = gmsfun(alldata.test$stop,alldata.test$status,5)
  BScores.g[j]=Int_Brierscore(1:n2,dipolartree.prune.g,Dipolar.model.g,IBSrange.g,0)
  actualtime.g = alldata.test$stop
  predictedtime2.g = Dipolar.model.g$predicttime(alldata.test,dipolartree.prune.g)
  CIndex.g[j] = Dipolar.model.g$cindex(actualtime.g,predictedtime2.g,alldata.test$status)
}

mean(BScores.l)
mean(CIndex.l)
mean(Num_nodes_dt.l)
mean(Num_nodes_dtp.l)

mean(BScores.q)
mean(CIndex.q)
mean(Num_nodes_dt.q)
mean(Num_nodes_dtp.q)

mean(BScores.g)
mean(CIndex.g)
mean(Num_nodes_dt.g)
mean(Num_nodes_dtp.g)














