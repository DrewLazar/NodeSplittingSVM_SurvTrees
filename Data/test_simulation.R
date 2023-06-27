# Remove preexisting variables from environment
rm(list=ls())
# Clear console
cat("\014")
# Force-set directory to this file in RStudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# Load the needed survsim scripts  
source("../Data/simple.surv.sim.int.R")
source("../Data/simple.ev.sim.R")
# Load the DipolarSurvivalTree class; all necessary libraries included in the class file
source("../Models/UnivariateSurvivaltree.R")
source("../Models/dipolarsurvivaltree.R")
source("../Models/quadraticdipolarsurvivaltree.R")
source("../extrafunctions/bootstrapPruning.R")
source('../extrafunctions/BrierScoreFun.R')
source("../extrafunctions/predictionfunction.R")


set.seed(46) 


#Generate some data
#Define parameters for Weibull distribution for censoring distribution
pc=1.2; b0c = 5.8;
anc.cens = pc; beta0.cens = b0c/pc;
#Define parameters for Weibull distribution for time-to-event distribution
pt = 7.8; b0t = 8.4; betat = list(2.8,6,7)
anc.ev=pt; beta0.ev = b0t/pt; beta =  lapply(betat,"*",1/pt)

sim.data <- simple.surv.sim.int(n=160, foltime=100, dist.ev=c('weibull'),
                                anc.ev,beta0.ev,dist.cens=c('weibull'),anc.cens,
                                beta0.cens,z=list(c("unif", 2, 5.5)), beta, x=list(c("bern", 0.3),c("normal", .6, .3),c("int",1,2)))

names(sim.data)[names(sim.data) == 'stop'] <- 'survt'


alldata = sim.data; covariates = c("x","x.1"); time = "survt"; censor = "status"; covariatestoaugment = list(c("x","x.1"),c("x.1","x.1")); quantiles = c(.20,.80); tolerance = 10^-3; epsilon = 1;
nsize = 10; allpairs = FALSE


#Find the percentage of censoring
n=nrow(alldata)
p_censored = sum(alldata$status)/n
print(p_censored)


#Instantiate our model classes 
Dipolar.model <- DipolarSurvivalTree$new(
  alldata, time, censor, covariates,
  quantiles, tolerance, epsilon, nsize,
  pureweight=1, mixedweight=1
)
Dipolarq.model <- QuadraticDipolarSurvivalTree$new(
  alldata, time, censor, covariates, covariatestoaugment,
  quantiles, tolerance, epsilon, nsize,
  pureweight=1, mixedweight=1, allpairs=allpairs
)
Univariate.model <- UnivariateSurvivalTree$new(
  alldata, time, censor, covariates,nsize)


#Create tree objects 
#fullsubset<-1:n
#dipolartree<-Dipolar.model$createtree(fullsubset)
#dipolarqtree<-Dipolarq.model$createtree(fullsubset)
#univariatetree<-Univariate.model$createtree(fullsubset)

#print full trees and log rank statistics  
#print(dipolartree,"lrstat")
#print(dipolarqtree,"lrstat")
#print(univariatetree,"lrstat")

# prune our trees using 25, 80% bootstrap samples each 
#univariatetree_prun=bootstrapPruning(univariatetree,Univariate.model,c(25,.9),1,4)
#dipolartree_prun=bootstrapPruning(dipolartree,Dipolar.model,c(25,.8),0,4)
#dipolartreeq_prun=bootstrapPruning(dipolarqtree,Dipolarq.model,c(25,.8),0,4)


#print pruned trees 
#print(dipolartree_prun,"lrstat")
#print(dipolartreeq_prun,"lrstat")
#print(univariatetree_prun,"lrstat")

#plot our pruned trees 
#plot(univariatetree_prun)
#plot(dipolartree_prun )
#plot(dipolartreeq_prun)

#Function to create Integrated Briar Scores times for computing IBS. Returns the 
#arithmetic means of uncensored survival observations. 
gmsfun = function(cforgms){
  cforgms=sort(cforgms[!duplicated(cforgms)])
  n=length(cforgms); a=c()
  for (i in 1:(n-1)){
    a[i]=sqrt(cforgms[i]*cforgms[i+1])
  }
  return(a)
}

#Concordance and Briar Scores of models pruned and unpruned, evaluated on full data  

#Concordance dipolartreeu and dipolartreeu_prun
#actualtime = alldata$survt
#predictedtime = Univariate.model$predicttime(alldata, univariatetree)
#Univariate.model$cindex(actualtime,predictedtime,alldata$status)
#predictedtime_prun = Univariate.model$predicttime(alldata, univariatetree_prun)
#Univariate.model$cindex(actualtime,predictedtime_prun,alldata$status)


#Concordance dipolartree and dipolartree_prun
#actualtime = alldata$survt
#predictedtime = Dipolar.model$predicttime(alldata, dipolartree)
#Dipolar.model$cindex(actualtime,predictedtime,alldata$status)
#predictedtime_prun = Dipolar.model$predicttime(alldata, dipolartree_prun)
#Dipolar.model$cindex(actualtime,predictedtime_prun,alldata$status)

#Concordance dipolartreeq
#actualtime = alldata$survt
#predictedtime2 = Dipolarq.model$predicttime(alldata, dipolarqtree)
#Dipolar.model$cindex(actualtime,predictedtime2,alldata$status)
#predictedtime_prun = Dipolarq.model$predicttime(alldata, dipolartreeq_prun)
#Dipolar.model$cindex(actualtime,predictedtime_prun,alldata$status)

#IBSrange for Brier Scores
#alleventtimes = alldata$survt[alldata$status==1]; IBSrange = gmsfun(alleventtimes)


#BrierScore2 dipolartreeu
#IntBS_u = Int_Brierscore(1:n,univariatetree,Univariate.model,IBSrange,1)
#IntBS_up = Int_Brierscore(1:n,univariatetree_prun,Univariate.model,IBSrange,1)

#BrierScore2 dipolartree
#IntBS_dt = Int_Brierscore(1:n,dipolartree,Dipolar.model,IBSrange,0)
#IntBS_dtp = Int_Brierscore(1:n,dipolartree_prun,Dipolar.model,IBSrange,0)

#BrierScore2 dipolartreeq
#IntBS_dtq = Int_Brierscore(1:n,dipolarqtree,Dipolarq.model,IBSrange,0)
#IntBS_dtqp = Int_Brierscore(1:n,dipolartreeq_prun,Dipolarq.model,IBSrange,0)


#k-fold cross-validation of Briar Scores and Concordance Indices 

#Make k fold partition function 
foldsfun =function(data,k){
  n=nrow(data)
  a = rep(floor(n/k),k) 
  b = rep(0:1,times=c(k-n%%k,n%%k)) 
  c = a + b 
  folds<-split(1:n, sample(rep(1:k,times=c)))
  return(folds)
}

#Create indices of k folds 
k=5
folds = foldsfun(alldata,k); BScores = c(); CIndex = c()


#Cross-validation for dipolar tree 
BScores<-c(); CIndex<-c(); Num_nodes_dt<-c(); Num_nodes_dtp<-c()
for (i in 1:5){
  trainsubset<-setdiff(1:nrow(alldata),folds[[i]])
  dipolartree.kfold<-Dipolar.model$createtree(trainsubset)
  Num_nodes_dt[i]<-dipolartree.kfold$totalCount
  dipolartree.kfold.prune<-bootstrapPruning(dipolartree.kfold,Dipolar.model,c(25,.8),0,3)
  Num_nodes_dtp[i]<-dipolartree.kfold.prune$totalCount
  alldatatest=Dipolar.model$traindata[folds[[i]],]
  IBSrange = gmsfun(alldatatest$survt)
  BScores[i]=Int_Brierscore(folds[[i]],dipolartree.kfold.prune,Dipolar.model,IBSrange,0)
  actualtime = alldatatest$survt
  predictedtime2 = Dipolar.model$predicttime(alldatatest,dipolartree.kfold.prune)
  CIndex[i] = Dipolar.model$cindex(actualtime,predictedtime2,alldatatest$status)
}
mean(BScores)
mean(CIndex)
mean(Num_nodes_dt)
mean(Num_nodes_dtp)



#Cross-validation for quadratic tree 
BScoresq<-c(); CIndexq<-c(); Num_nodes_dtq<-c(); Num_nodes_dtpq<-c()
for (i in 1:k){
  trainsubset<-setdiff(1:nrow(alldata),folds[[i]])
  dipolartreeq.kfold<-Dipolarq.model$createtree(trainsubset)
  Num_nodes_dtq[i]<-dipolartreeq.kfold$totalCount
  dipolartreeq.kfold.prune<-bootstrapPruning(dipolartreeq.kfold,Dipolarq.model,c(25,.8),0,3)
  Num_nodes_dtpq[i]<-dipolartreeq.kfold.prune$totalCount
  alldatatest=Dipolarq.model$traindata[folds[[i]],]
  IBSrange = gmsfun(alldatatest$survt)
  BScoresq[i]=Int_Brierscore(1:n,dipolartreeq.kfold.prune,Dipolarq.model,IBSrange,0)
  actualtime = alldatatest$survt
  predictedtime = Dipolarq.model$predicttime(alldatatest,dipolartreeq.kfold.prune)
  CIndexq[i] = Dipolar.model$cindex(actualtime,predictedtime,alldatatest$status)
}
mean(BScoresq)
mean(CIndexq)
mean(Num_nodes_dtq)
mean(Num_nodes_dtpq)

#Cross-validation for univariate tree 

BScoresu<-c(); CIndexu<-c(); Num_nodes_dtu<-c(); Num_nodes_dtpu<-c()

for (i in 1:k){
  trainsubset<-setdiff(1:nrow(alldata),folds[[i]])
  univariate.kfold<-Univariate.model$createtree(trainsubset)
  Num_nodes_dtu[i]<- univariate.kfold$totalCount
  univariate.kfold.prune<-bootstrapPruning(univariate.kfold,Univariate.model,c(25,.8),1,3)
  Num_nodes_dtpu[i]<-univariate.kfold.prune$totalCount
  alldatatest=Univariate.model$traindata[folds[[i]],]
  IBSrange = gmsfun(alldatatest$survt)
  BScoresu[i]=Int_Brierscore(folds[[i]],univariate.kfold.prune,Univariate.model,IBSrange,1)
  actualtime = alldatatest$survt
  predictedtime = Univariate.model$predicttime(alldatatest,univariate.kfold.prune)
  CIndexu[i] = Univariate.model$cindex(actualtime,predictedtime,alldatatest$status)
}
mean(BScoresu)
mean(CIndexu)
mean(Num_nodes_dtu)
mean(Num_nodes_dtpu)

