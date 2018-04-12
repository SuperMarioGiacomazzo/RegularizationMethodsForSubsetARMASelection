options(scipen=999)

#Necessary Packages
library(bayesreg)
library(MCMCpack)
library(doParallel)
library(foreach)
library(forecast)

#Set Working Directory
setwd("D:/Mario Documents/Graduate School/Github Repositories/RegularizationMethodsForSubsetARMASelection")

source("subsetARMACode.R")

#Simulation Specs
S=500  #Number of Replications of Simulated Time Series
N=c(120,240,360)   #Lengths of Simulated Time Series
noise=c(0.5,1,1.5) #Inherent Noise in Simulated Time Series
maxP=14 #Maximum Autoregressive Order
maxQ=14 #Maximum Moving Average Order
true.ar=rep(0,14)
true.ma=c(0.8,rep(0,4),0.7,0.56,rep(0,7))
coef.true=c(true.ar,true.ma)

#Methods Considered
methods=c("ADLASSO AIC","ADLASSO AIC/BIC","ADLASSO BIC")

#Empty Dataframe to Save Results
RESULTS=data.frame(matrix(NA,1,8+maxP+maxQ))
names(RESULTS)=c("Method","Length","Noise","Sim",
                 paste("AR(",1:maxP,")",sep=""),
                 paste("MA(",1:maxP,")",sep=""),
                 "P","A","FN","FP")


for(l in 1:S){
  for(k in 1:length(noise)){
    for(j in 1:length(N)){
      
      ####################################################################################################################################
      #Simulate Dataset
      model<-Arima(ts(rnorm(N[j]+1000,0,1),freq=6),order=c(0,0,1),seasonal=c(0,0,1),include.mean=F,fixed=c(0.8,0.7))
      x<-simulate(model,nsim=(N[j]+1000),future=F,innov = rnorm((1000+N[j]),0,noise[k]))[-(1:1000)]
      ####################################################################################################################################
      

      
      ####################################################################################################################################
      #Estimate Using Adaptive Enet Selecting Tuning Parameter Via Minimization of AIC/BIC
      
      ##Use AIC first and AIC second (no update to moving average terms)
      adlasso.aic=cv.adaptlasso.arma.func(x,h=1,cores=1,long.ar.select=T,maxP=maxP,maxQ=maxQ,max.pq=NULL,K=NULL,updateMA=F,
                                        test.per=0.2,BIC1=F,BIC2=F,eta=2,CV.method=c("AIC/BIC"))
      coef.adlasso.aic=round(adlasso.aic$final.mod.coef,4)
      names(coef.adlasso.aic)=c(paste("ar",1:14,sep=""),paste("ma",1:14,sep=""))
      coef.adlasso.aic=c(coef.adlasso.aic,fulleval.func(truecoef=coef.true,estcoef=coef.adlasso.aic))
      
      ##Use AIC first and BIC second (no update to moving average terms)
      adlasso.aicbic=cv.adaptlasso.arma.func(x,h=1,cores=1,long.ar.select=T,maxP=maxP,maxQ=maxQ,max.pq=NULL,K=NULL,updateMA=F,
                                           test.per=0.2,BIC1=F,BIC2=T,eta=2,CV.method=c("AIC/BIC"))
      coef.adlasso.aicbic=round(adlasso.aicbic$final.mod.coef,4)
      names(coef.adlasso.aicbic)=c(paste("ar",1:14,sep=""),paste("ma",1:14,sep=""))
      coef.adlasso.aicbic=c(coef.adlasso.aicbic,fulleval.func(truecoef=coef.true,estcoef=coef.adlasso.aicbic))
      
      ##Use BIC first and BIC second (no update to moving average terms)
      adlasso.bic=cv.adaptlasso.arma.func(x,h=1,cores=1,long.ar.select=T,maxP=maxP,maxQ=maxQ,max.pq=NULL,K=NULL,updateMA=F,
                                        test.per=0.2,BIC1=T,BIC2=T,eta=2,CV.method=c("AIC/BIC"))
      coef.adlasso.bic=round(adlasso.bic$final.mod.coef,4)
      names(coef.adlasso.bic)=c(paste("ar",1:14,sep=""),paste("ma",1:14,sep=""))
      coef.adlasso.bic=c(coef.adlasso.bic,fulleval.func(truecoef=coef.true,estcoef=coef.adlasso.bic))
      ####################################################################################################################################
      
      
      
      save.image("ADLASSOSIM3AICBICB.Rdata")
      
      
      ####################################################################################################################################
      #Save Results into Data Frame
      init.RESULTS=as.data.frame(rbind(
        coef.adlasso.aic,
        coef.adlasso.aicbic,
        coef.adlasso.bic
      ))
      row.names(init.RESULTS)=NULL
      init.RESULTS2=cbind(methods,N[j],noise[k],l,init.RESULTS)
      names(init.RESULTS2)=names(RESULTS)
      RESULTS=rbind(RESULTS,init.RESULTS2)
      ####################################################################################################################################
      
      save.image("ADLASSOSIM3AICBICB.Rdata")
    }
  }
}

save.image("ADLASSOSIM3AICBICB.Rdata")
