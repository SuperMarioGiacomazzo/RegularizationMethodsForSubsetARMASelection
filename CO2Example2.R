options(scipen=9999)

#Necessary Packages
library(bayesreg)
library(doParallel)
library(foreach)
library(forecast)
library(TSA)
library(gplots)

#Set Working Directory
setwd("D:/Mario Documents/Graduate School/Ph.D Dissertation/ChenChanSimulation")

#Function Code for Different Methods
source("subsetARMACode.R")

################################
#Real Application Using CO2 Data
################################

data(co2)

#Obtain Dataset and Apply Seasonal and Single Lag Differencing
original.co2=as.numeric(co2)
original.time=as.numeric(time(co2))
seasdiff.co2=c(rep(NA,12),diff(original.co2,12))
final.co2=c(rep(NA,1),diff(seasdiff.co2,1))

#Split Into Modeling and Validation Periods
MODEL.PERIOD=original.time<2003
VALIDATION.PERIOD=original.time>=2003



#Plot of Modified time Series
png(filename="co2plots2.png",width=500,height=700)
par(mfrow=c(3,1),mar=c(2.5,5,1,1),cex.axis=2,cex.lab=2)
plot(x=original.time,y=original.co2,xlab="",type="n",ylab="Original")
points(x=original.time[MODEL.PERIOD],y=original.co2[MODEL.PERIOD],type="l",col="black")
points(x=original.time[VALIDATION.PERIOD],y=original.co2[VALIDATION.PERIOD],type="l",col="black",lty=3)
plot(x=original.time,y=seasdiff.co2,xlab="",type="n",ylab="Seasonal Difference")
points(x=original.time[MODEL.PERIOD],y=seasdiff.co2[MODEL.PERIOD],type="l",col="black")
points(x=original.time[VALIDATION.PERIOD],y=seasdiff.co2[VALIDATION.PERIOD],type="l",col="black",lty=3)
plot(x=original.time,y=final.co2,xlab="",type="n",ylab="Regular Difference")
points(x=original.time[MODEL.PERIOD],y=final.co2[MODEL.PERIOD],type="l",col="black")
points(x=original.time[VALIDATION.PERIOD],y=final.co2[VALIDATION.PERIOD],type="l",col="black",lty=3)
dev.off()

#Split Data Into Modeling Period
train.co2=as.numeric(na.omit(final.co2[MODEL.PERIOD]))
train.mae=mean(abs(diff(train.co2)))

MODELS=list()



#ADLASSO (Methods: m=1,2,3,4,6,8,9,11) (Total:8 methods)

######################################################################################################################################
#Estimate ARMA Using Adaptive Lasso Using AIC for Both Steps
MODELS$adlasso.mod.AIC=cv.adaptlasso.arma.func(x=train.co2,h=1,long.ar.select=F,maxP=14,maxQ=14,max.pq=14,K=NULL,updateMA=F,
                                               test.per=0.2,BIC1=F,BIC2=F,eta=2,CV.method="AIC/BIC")

#Estimate ARMA Using Adaptive Lasso Using AIC in first step and BIC in second step
MODELS$adlasso.mod.AICBIC=cv.adaptlasso.arma.func(x=train.co2,h=1,long.ar.select=F,maxP=14,maxQ=14,max.pq=14,K=NULL,updateMA=F,
                                                  test.per=0.2,BIC1=F,BIC2=T,eta=2,CV.method="AIC/BIC")

#Estimate ARMA Using Adaptive Lasso Using BIC for Both Steps
MODELS$adlasso.mod.BIC=cv.adaptlasso.arma.func(x=train.co2,h=1,long.ar.select=F,maxP=14,maxQ=14,max.pq=14,K=NULL,updateMA=F,
                                               test.per=0.2,BIC1=T,BIC2=T,eta=2,CV.method="AIC/BIC")
######################################################################################################################################


######################################################################################################################################
#Estimate ARMA Using Adaptive Lasso Using BIC for Both Steps
MODELS$adlasso.mod.oosindep=cv.adaptlasso.arma.func(x=train.co2,h=1,long.ar.select=F,maxP=14,maxQ=14,max.pq=14,K=10,updateMA=F,
                                                    test.per=0.2,BIC1=T,BIC2=T,eta=2,CV.method="OOSINDEP")
######################################################################################################################################


######################################################################################################################################
#Estimate ARMA Using Adaptive Lasso With 5 Fold Regular Cross Validation
MODELS$adlasso.mod.cv5=cv.adaptlasso.arma.func(x=train.co2,h=1,long.ar.select=F,maxP=14,maxQ=14,max.pq=14,K=10,updateMA=F,
                                               test.per=0.2,BIC1=T,BIC2=T,eta=2,CV.method="RegCV")
#Estimate ARMA Using Adaptive Lasso With Leave-one-out Regular Cross Validation
MODELS$adlasso.mod.loocv=cv.adaptlasso.arma.func(x=train.co2,h=1,long.ar.select=F,maxP=14,maxQ=14,max.pq=14,K=NULL,updateMA=F,
                                                 test.per=0.2,BIC1=T,BIC2=T,eta=2,CV.method="RegCV")
######################################################################################################################################


######################################################################################################################################
#Estimate ARMA Using Adaptive Lasso With 10 Fold NonDependent Blocked Cross Validation
MODELS$adlasso.mod.bcv5=cv.adaptlasso.arma.func(x=train.co2,h=1,long.ar.select=F,maxP=14,maxQ=14,max.pq=14,K=5,updateMA=F,
                                                test.per=0.2,BIC1=F,BIC2=T,eta=2,CV.method="Bergmeir")
#Estimate ARMA Using Adaptive Lasso With Leave-one-block-out Cross Validation
MODELS$adlasso.mod.lobocv=cv.adaptlasso.arma.func(x=train.co2,h=1,long.ar.select=F,maxP=14,maxQ=14,max.pq=14,K=NULL,updateMA=F,
                                                  test.per=0.2,BIC1=F,BIC2=F,eta=2,CV.method="LOBOCV")
######################################################################################################################################



#ADENET (Methods: m=1,2,3,4,6,8,9,11) (Total:8 methods)

######################################################################################################################################
#Estimate ARMA Using Adaptive Enet Using AIC for Both Steps
MODELS$adnet.mod.AIC=cv.adaptenet.arma.func(x=train.co2,h=1,long.ar.select=F,maxP=14,maxQ=14,max.pq=14,K=NULL,updateMA=F,
                                            test.per=0.2,BIC1=F,BIC2=F,eta=2,CV.method="AIC/BIC",alpha=seq(0,1,0.1))
#Estimate ARMA Using Adaptive Enet Using AIC in first step and BIC in second step
MODELS$adnet.mod.AICBIC=cv.adaptenet.arma.func(x=train.co2,h=1,long.ar.select=F,maxP=14,maxQ=14,max.pq=14,K=NULL,updateMA=F,
                                               test.per=0.2,BIC1=F,BIC2=T,eta=2,CV.method="AIC/BIC",alpha=seq(0,1,0.1))
#Estimate ARMA Using Adaptive Enet Using BIC for Both Steps
MODELS$adnet.mod.BIC=cv.adaptenet.arma.func(x=train.co2,h=1,long.ar.select=F,maxP=14,maxQ=14,max.pq=14,K=NULL,updateMA=F,
                                            test.per=0.2,BIC1=T,BIC2=T,eta=2,CV.method="AIC/BIC",alpha=seq(0,1,0.1))
######################################################################################################################################


######################################################################################################################################
#Estimate ARMA Using Independent OOS
MODELS$adenet.mod.oosindep=cv.adaptenet.arma.func(x=train.co2,h=1,long.ar.select=F,maxP=14,maxQ=14,max.pq=14,K=10,updateMA=F,
                                                  test.per=0.2,BIC1=T,BIC2=T,eta=2,CV.method="OOSINDEP",alpha=seq(0,1,0.1),CV2=c("min"))
######################################################################################################################################


######################################################################################################################################
#Estimate ARMA Using Adaptive Enet With 5 Fold Regular Cross Validation
MODELS$adenet.mod.cv5=cv.adaptenet.arma.func(x=train.co2,h=1,long.ar.select=F,maxP=14,maxQ=14,max.pq=14,K=5,updateMA=F,
                                             test.per=0.2,BIC1=T,BIC2=T,eta=2,CV.method="RegCV",alpha=seq(0,1,0.1),CV2=c("min"))
#Estimate ARMA Using Adaptive Enet With Leave-one-out Regular Cross Validation
MODELS$adenet.mod.loocv=cv.adaptenet.arma.func(x=train.co2,h=1,long.ar.select=F,maxP=14,maxQ=14,max.pq=14,K=NULL,updateMA=F,
                                               test.per=0.2,BIC1=T,BIC2=T,eta=2,CV.method="RegCV",alpha=seq(0,1,0.1),CV2=c("min"))
######################################################################################################################################


######################################################################################################################################
#Estimate ARMA Using Adaptive Enet With 10 Fold NonDependent Blocked Cross Validation
MODELS$adenet.mod.bcv5=cv.adaptenet.arma.func(x=train.co2,h=1,long.ar.select=F,maxP=14,maxQ=14,max.pq=14,K=5,updateMA=F,
                                              test.per=0.2,BIC1=F,BIC2=T,eta=2,CV.method="Bergmeir",alpha=seq(0,1,0.1),CV2=c("min"))
#Estimate ARMA Using Adaptive Enet With Leave-one-block-out Cross Validation
MODELS$adenet.mod.lobocv=cv.adaptenet.arma.func(x=train.co2,h=1,long.ar.select=F,maxP=14,maxQ=14,max.pq=14,K=NULL,updateMA=F,
                                                test.per=0.2,BIC1=F,BIC2=F,eta=2,CV.method="LOBOCV",alpha=seq(0,1,0.1),CV2=c("min"))
######################################################################################################################################



#BHS (Methods: m=1,2,3,4) (Total:4 methods)

######################################################################################################################################
pms.mod1=pms.arma.func(x=train.co2,h=1,maxP=14,maxQ=14,KL.threshold=c(0.90,0.95,0.98),prior.choice="hs",updateMA=F)

MODELS$pms.mod.hs.90=list(final.mod.coef=pms.mod1$final.mod.coef[[1]],final.mod.int=pms.mod1$final.mod.int[[1]],
                          final.mod.s2=pms.mod1$final.mod.s2[[1]],nonzero.select=as.numeric(na.omit(pms.mod1$KL.select[,1]-1))[-1])

MODELS$pms.mod.hs.95=list(final.mod.coef=pms.mod1$final.mod.coef[[2]],final.mod.int=pms.mod1$final.mod.int[[2]],
                          final.mod.s2=pms.mod1$final.mod.s2[[2]],nonzero.select=as.numeric(na.omit(pms.mod1$KL.select[,2]-1))[-1])

MODELS$pms.mod.hs.98=list(final.mod.coef=pms.mod1$final.mod.coef[[3]],final.mod.int=pms.mod1$final.mod.int[[3]],
                          final.mod.s2=pms.mod1$final.mod.s2[[3]],nonzero.select=as.numeric(na.omit(pms.mod1$KL.select[,3]-1))[-1])
######################################################################################################################################


######################################################################################################################################
indhs<-cv.pms.arma.func(x=train.co2,h=1,maxP=14,maxQ=14,KL.stop=0.98,cv.type="ind",
                        test.per=0.2,prior.choice="hs",updateMA=F)
MODELS$pms.mod.indhs=list(final.mod.coef=indhs$final.mod.coef,final.mod.int=indhs$final.mod.int,
                          final.mod.s2=indhs$final.mod.s2,nonzero.select=which(indhs$final.mod.coef!=0))
######################################################################################################################################



#BHS+ (Methods: m=1,2,3,4) (Total:4 methods)

######################################################################################################################################
pms.mod2=pms.arma.func(x=train.co2,h=1,maxP=14,maxQ=14,KL.threshold=c(0.90,0.95,0.98),prior.choice="hs+",updateMA=F)

MODELS$pms.mod.hsp.90=list(final.mod.coef=pms.mod2$final.mod.coef[[1]],final.mod.int=pms.mod2$final.mod.int[[1]],
                           final.mod.s2=pms.mod2$final.mod.s2[[1]],nonzero.select=as.numeric(na.omit(pms.mod2$KL.select[,1]-1))[-1])

MODELS$pms.mod.hsp.95=list(final.mod.coef=pms.mod2$final.mod.coef[[2]],final.mod.int=pms.mod2$final.mod.int[[2]],
                           final.mod.s2=pms.mod2$final.mod.s2[[2]],nonzero.select=as.numeric(na.omit(pms.mod2$KL.select[,2]-1))[-1])

MODELS$pms.mod.hsp.98=list(final.mod.coef=pms.mod2$final.mod.coef[[3]],final.mod.int=pms.mod2$final.mod.int[[3]],
                           final.mod.s2=pms.mod2$final.mod.s2[[3]],nonzero.select=as.numeric(na.omit(pms.mod2$KL.select[,3]-1))[-1])
######################################################################################################################################


######################################################################################################################################
indhsp<-cv.pms.arma.func(x=train.co2,h=1,maxP=14,maxQ=14,KL.stop=0.98,cv.type="ind",
                         test.per=0.2,prior.choice="hs+",updateMA=F)
MODELS$pms.mod.indhsp=list(final.mod.coef=indhsp$final.mod.coef,final.mod.int=indhsp$final.mod.int,
                           final.mod.s2=indhsp$final.mod.s2,nonzero.select=which(indhsp$final.mod.coef!=0))
######################################################################################################################################












######################################################################################################################################
#Obtain Forecasts for All Models Estimated on Training Data

nMODELS=length(MODELS)
SELECT=matrix(0,nMODELS,ncol=28)
COEF=matrix(0,nMODELS,ncol=28)
RMSFE=rep(NA,nMODELS)
MASE=rep(NA,nMODELS)
MFB=rep(NA,nMODELS)
MDFB=rep(NA,nMODELS)
FC.LOWER=matrix(NA,sum(VALIDATION.PERIOD),nMODELS)
FC.MEAN=matrix(NA,sum(VALIDATION.PERIOD),nMODELS)
FC.UPPER=matrix(NA,sum(VALIDATION.PERIOD),nMODELS)

for(k in 1:nMODELS){
  temp.mod=MODELS[[k]]
  SELECT[k,temp.mod$nonzero.select]=1
  COEF[k,]=temp.mod$final.mod.coef
  
  fc.lower=rep(NA,length(final.co2))
  fc.mean=rep(NA,length(final.co2))
  fc.upper=rep(NA,length(final.co2))
  
  error=rep(0,length(final.co2))
  
  for(j in 30:length(final.co2)){
    fc=as.numeric(final.co2[(j-1):(j-14)]%*%temp.mod$final.mod.coef[1:14]+
                    error[(j-1):(j-14)]%*%temp.mod$final.mod.coef[-(1:14)]) +
      rnorm(100000,mean=temp.mod$final.mod.int,sd=sqrt(temp.mod$final.mod.s2))
    fc.lower[j]=quantile(fc,0.05)
    fc.mean[j]=mean(fc)
    fc.upper[j]=quantile(fc,0.95)
    error[j]=final.co2[j]-fc.mean[j]
  }
  
  RMSFE[k]=sqrt(mean((final.co2[VALIDATION.PERIOD]-fc.mean[VALIDATION.PERIOD])^2))
  MASE[k]=mean(abs((final.co2[VALIDATION.PERIOD]-fc.mean[VALIDATION.PERIOD])/train.mae))
  MFB[k]=mean(final.co2[VALIDATION.PERIOD]-fc.mean[VALIDATION.PERIOD])
  MDFB[k]=(sum((final.co2[VALIDATION.PERIOD]-fc.mean[VALIDATION.PERIOD])>0)-
             sum((final.co2[VALIDATION.PERIOD]-fc.mean[VALIDATION.PERIOD])<0))/sum(VALIDATION.PERIOD)
  
  FC.LOWER[,k]=fc.lower[VALIDATION.PERIOD]
  FC.MEAN[,k]=fc.mean[VALIDATION.PERIOD]
  FC.UPPER[,k]=fc.upper[VALIDATION.PERIOD]
  
}
######################################################################################################################################




######################################################################################################################################
#Create Table of Results
OUT.RESULTS=cbind(c(1,2,3,4,6,8,9,11,1,2,3,4),
                  RMSFE[c(1:8,17:20)],RMSFE[-c(1:8,17:20)],
                  MASE[c(1:8,17:20)],MASE[-c(1:8,17:20)],
                  MFB[c(1:8,17:20)],MFB[-c(1:8,17:20)],
                  MDFB[c(1:8,17:20)],MDFB[-c(1:8,17:20)])

OUT.RESULTS2=round(OUT.RESULTS,2)

#Identify Bad Forecasts Based on Naive Model
BAD.PREDICT=lag.func(final.co2,1)
BAD.RMSFE=sqrt(mean((final.co2[VALIDATION.PERIOD]-BAD.PREDICT[VALIDATION.PERIOD])^2))
BAD.MASE=mean(abs((final.co2[VALIDATION.PERIOD]-BAD.PREDICT[VALIDATION.PERIOD])/train.mae))
BAD.MFB=mean(final.co2[VALIDATION.PERIOD]-BAD.PREDICT[VALIDATION.PERIOD])
BAD.MDFB=(sum((final.co2[VALIDATION.PERIOD]-BAD.PREDICT[VALIDATION.PERIOD])>0)-
            sum((final.co2[VALIDATION.PERIOD]-BAD.PREDICT[VALIDATION.PERIOD])<0))/sum(VALIDATION.PERIOD)
OUT.RESULTS3=rbind(OUT.RESULTS2,round(c(100,BAD.RMSFE,BAD.RMSFE,BAD.MASE,BAD.MASE,BAD.MFB,BAD.MFB,BAD.MDFB,BAD.MDFB),2))

#Check Stationarity and Invertibility of Estimates
stationarity=rep(NA,nMODELS)
invertibility=rep(NA,nMODELS)
for(k in 1:nMODELS){
  ar.coef=COEF[k,1:14]
  st.poly=c(1,-ar.coef)
  st.root=polyroot(st.poly)
  st.check=sum(abs(st.root)>1)==length(st.root)
  stationarity[k]=st.check
  
  ma.coef=COEF[k,-(1:14)]
  ma.poly=c(1,ma.coef)
  ma.root=polyroot(ma.poly)
  ma.check=sum(abs(ma.root)>1)==length(ma.root)
  invertibility[k]=ma.check
}

#Print Out Table
OUT.RESULTS.TAB=xtable(OUT.RESULTS3)
align(OUT.RESULTS.TAB)=rep("c",10)
print(OUT.RESULTS.TAB,include.rownames=F,include.colnames=F)
######################################################################################################################################



######################################################################################################################################
#Heatmap
colors <- c("maroon", "gold")
SELECT2=melt(SELECT)
SELECT2$value=factor(SELECT2$value)
levels(SELECT2$value)=c("Out ","In")

png(filename="alertselect.png",width=800,height=800)
ggplot(SELECT2, aes(x = factor(Var2), y = factor(Var1,levels=24:1), fill = factor(value))) + 
  geom_tile() + xlab("")+ ylab("") +
  scale_fill_manual("Final Model Selection  ",values=colors) +
  geom_vline(xintercept=seq(0.5,28.5,1),color="white",size=0.1)+
  geom_hline(yintercept=seq(0.5,24.5,1),color="white",size=0.1)+
  scale_x_discrete(breaks=1:28,labels=c(expression(phi[1]),expression(phi[2]),
                                        expression(phi[3]),expression(phi[4]),expression(phi[5]),expression(phi[6]),
                                        expression(phi[7]),expression(phi[8]),expression(phi[9]),expression(phi[10]),
                                        expression(phi[11]),expression(phi[12]),expression(phi[13]),expression(phi[14]),
                                        expression(theta[1]),expression(theta[2]),
                                        expression(theta[3]),expression(theta[4]),expression(theta[5]),expression(theta[6]),
                                        expression(theta[7]),expression(theta[8]),expression(theta[9]),expression(theta[10]),
                                        expression(theta[11]),expression(theta[12]),expression(theta[13]),expression(theta[14])))+
  scale_y_discrete(breaks=1:24,labels=c(expression("AL"[1]),expression("AL"[2]),
                                        expression("AL"[3]),expression("AL"[4]),expression("AL"[6]),
                                        expression("AL"[8]),expression("AL"[9]),expression("AL"[11]),
                                        expression("AE"[1]),expression("AE"[2]),
                                        expression("AE"[3]),expression("AE"[4]),expression("AE"[6]),
                                        expression("AE"[8]),expression("AE"[9]),expression("AE"[11]),
                                        expression("BHS"[1]),expression("BHS"[2]),
                                        expression("BHS"[3]),expression("BHS"[4]),
                                        expression("BHS"[1]^"+"),expression("BHS"[2]^"+"),
                                        expression("BHS"[3]^"+"),expression("BHS"[4]^"+")))+
  theme(legend.position="bottom",line=element_blank(),legend.title=element_text(size=16),
        plot.margin = unit(c(0.5, 0.5,1, 0), "cm"),legend.justification="center",legend.text=element_text(size=14))+
  theme(axis.text.x=element_text(size=14),axis.text.y=element_text(size=14))
dev.off()
######################################################################################################################################



######################################################################################################################################
#Forecasts from auto.arima
REF.ARMA.MOD=auto.arima(train.co2,d=0,D=0,max.p=14,max.q=14,seasonal=F,ic="aic")
REF.ARMA.FIT=Arima(final.co2,model=REF.ARMA.MOD)
REF.ARMA.FC=REF.ARMA.FIT$fitted[VALIDATION.PERIOD]
REF.ARMA.RMSFE=sqrt(mean((final.co2[VALIDATION.PERIOD]-REF.ARMA.FC)^2))
REF.ARMA.MASE=mean(abs((final.co2[VALIDATION.PERIOD]-REF.ARMA.FC)/train.mae))
REF.ARMA.MFB=mean(final.co2[VALIDATION.PERIOD]-REF.ARMA.FC)
REF.ARMA.MDFB=(sum((final.co2[VALIDATION.PERIOD]-REF.ARMA.FC)>0)-
                 sum((final.co2[VALIDATION.PERIOD]-REF.ARMA.FC)<0))/sum(VALIDATION.PERIOD)
c(REF.ARMA.RMSFE,REF.ARMA.MASE,REF.ARMA.MFB,REF.ARMA.MDFB)

REF.SARMA.MOD=auto.arima(ts(train.co2,freq=12),d=0,D=0,max.P=14,max.Q=14,max.p=14,max.q=14,seasonal=T,ic="aic")
REF.SARMA.FIT=Arima(final.co2,model=REF.SARMA.MOD)
REF.SARMA.FC=REF.SARMA.FIT$fitted[VALIDATION.PERIOD]
REF.SARMA.RMSFE=sqrt(mean((final.co2[VALIDATION.PERIOD]-REF.SARMA.FC)^2))
REF.SARMA.MASE=mean(abs((final.co2[VALIDATION.PERIOD]-REF.SARMA.FC)/train.mae))
REF.SARMA.MFB=mean(final.co2[VALIDATION.PERIOD]-REF.SARMA.FC)
REF.SARMA.MDFB=(sum((final.co2[VALIDATION.PERIOD]-REF.SARMA.FC)>0)-
                  sum((final.co2[VALIDATION.PERIOD]-REF.SARMA.FC)<0))/sum(VALIDATION.PERIOD)
c(REF.SARMA.RMSFE,REF.SARMA.MASE,REF.SARMA.MFB,REF.SARMA.MDFB)

######################################################################################################################################






######################################################################################################################################
MODELS2=list()
adshrink.check.func<-function(x,h=1,long.ar.select=F,maxP,maxQ,updateMA=F,
                           BIC1=T,BIC2=T,eta=2,alpha=seq(0,1,0.1),
                           Method=c("ADLASSO","ADENET")){
  
  #Package Required
  require(glmnet) #Performs Ridge, Lasso, Elastic Net Estimation
  
  Method=match.arg(Method)
  
  Nt=length(x) #Length of Input Time Series
  
  #Fit Long AR Model to Estimate Innovations
  max.ar.order=ceiling(10*log10(Nt)) #Maximum Autoregressive Order
  init.mod.est=ar(x,aic=long.ar.select, #Allows Stepwise Selection
                  order.max=max.ar.order,demean=T)
  init.mod.error=residuals(init.mod.est)
  init.mod.order=length(which(is.na(init.mod.error)))
  
  #Create Model Matrix of AR and MA terms
  dataP=foreach(p=1:maxP,.combine=cbind)%do%{
    lag.func(x,k=(p+h-1))
  }
  dataQ=foreach(q=1:maxQ,.combine=cbind)%do%{
    lag.func(init.mod.error,k=(q+h-1))
  }
  first.modX=as.matrix(cbind(dataP,
                             dataQ))[-(1:(init.mod.order+max(maxP,maxQ)+h-1)),]
  first.y=x[-(1:(init.mod.order+max(maxP,maxQ)+h-1))]
  
  #Number of Alphas For Elastic Net
  n.alpha=length(alpha)
  
  #Estimation Via ADLASSO
  if(Method=="ADLASSO"){
    #Initial LASSO Weights
    first.mod.est=glmnet(y=first.y,x=first.modX,standardize=T,alpha=1)
    first.mod.RSS=colSums((first.y-predict(first.mod.est,
                                           newx=first.modX))^2)
    if(BIC1){
      first.bic.out=log(length(first.y))*first.mod.est$df+
        length(first.y)*log(as.vector(first.mod.RSS)/length(first.y))
    }else{
      first.aic.out=2*first.mod.est$df+
        length(first.y)*log(as.vector(first.mod.RSS)/length(first.y))
    }
    
    if(BIC1){
      stationarity=rep(NA,length(first.bic.out))
      invertibility=rep(NA,length(first.bic.out))
      for(k in 1:length(first.bic.out)){
        ar.coef=as.numeric(coef(first.mod.est,s=first.mod.est$lambda[k],method="lambda"))[2:15]
        st.poly=c(1,-ar.coef)
        st.root=polyroot(st.poly)
        st.check=sum(abs(st.root)>1)==length(st.root)
        stationarity[k]=st.check
        
        ma.coef=as.numeric(coef(first.mod.est,s=first.mod.est$lambda[k],method="lambda"))[16:29]
        ma.poly=c(1,ma.coef)
        ma.root=polyroot(ma.poly)
        ma.check=sum(abs(ma.root)>1)==length(ma.root)
        invertibility[k]=ma.check
      }    
    }else{
      stationarity=rep(NA,length(first.aic.out))
      invertibility=rep(NA,length(first.aic.out))
      for(k in 1:length(first.aic.out)){
        ar.coef=as.numeric(coef(first.mod.est,s=first.mod.est$lambda[k],method="lambda"))[2:15]
        st.poly=c(1,-ar.coef)
        st.root=polyroot(st.poly)
        st.check=sum(abs(st.root)>1)==length(st.root)
        stationarity[k]=st.check
        
        ma.coef=as.numeric(coef(first.mod.est,s=first.mod.est$lambda[k],method="lambda"))[16:29]
        ma.poly=c(1,ma.coef)
        ma.root=polyroot(ma.poly)
        ma.check=sum(abs(ma.root)>1)==length(ma.root)
        invertibility[k]=ma.check
      }
    }
      
    if(BIC1){
      which.good=which(stationarity & invertibility)
      good.bic=first.bic.out[which.good]
      good.lambdas=first.mod.est$lambda[which.good]
      first.mod.lambda=good.lambdas[which.min(good.bic)]
    }else{
      which.good=which(stationarity & invertibility)
      good.aic=first.aic.out[which.good]
      good.lambdas=first.mod.est$lambda[which.good]
      first.mod.lambda=good.lambdas[which.min(good.aic)]
    } 
      
    first.mod.coef=as.numeric(coef(first.mod.est,
                                   s=first.mod.lambda,method="lambda"))[-1]
    first.mod.mu=as.numeric(coef(first.mod.est,
                                 s=first.mod.lambda,method="lambda"))[1]
    weights=abs(first.mod.coef+1/length(first.y))^(-eta)
    
    #Update Model Matrix of AR and MA terms Based off Initial Estimation
    if(updateMA){
      update.mod.predict=rep(NA,length(x))
      update.mod.error=rep(0,length(x))
      for(v in (h+max(maxP,maxQ)):Nt){
        update.mod.predict[v]=first.mod.mu+
          x[(v-h):(v-maxP-h+1)]%*%first.mod.coef[1:maxP]+
          update.mod.error[(v-h):(v-maxQ-h+1)]%*%first.mod.coef[-(1:maxP)]
        update.mod.error[v]=x[v]-update.mod.predict[v]
      }
      update.dataQ=foreach(q=1:maxQ,.combine=cbind)%do%{
        lag.func(update.mod.error,k=(q+h-1))
      }
      second.modX=as.matrix(cbind(dataP,
                                  update.dataQ))[-(1:(max(maxP,maxQ)+h-1)),]
      second.y=x[-(1:(max(maxP,maxQ)+h-1))]
    }else{
      second.modX=first.modX
      second.y=first.y
    }
    
    second.mod.est=glmnet(y=second.y,x=second.modX,standardize=T,alpha=1,
                          thresh=1e-16, penalty.factor=weights)
    second.mod.RSS=colSums((second.y-predict(second.mod.est,
                                             newx=second.modX))^2)
    if(BIC2){
      second.bic.out=log(length(second.y))*second.mod.est$df+
        length(second.y)*log(as.vector(second.mod.RSS)/length(second.y))
    }else{
      second.aic.out=2*second.mod.est$df+
        length(second.y)*log(as.vector(second.mod.RSS)/length(second.y))
    }
    
    if(BIC2){
      stationarity=rep(NA,length(second.bic.out))
      invertibility=rep(NA,length(second.bic.out))
      for(k in 1:length(second.bic.out)){
        ar.coef=as.numeric(coef(second.mod.est,s=second.mod.est$lambda[k],method="lambda"))[2:15]
        st.poly=c(1,-ar.coef)
        st.root=polyroot(st.poly)
        st.check=sum(abs(st.root)>1)==length(st.root)
        stationarity[k]=st.check
        
        ma.coef=as.numeric(coef(second.mod.est,s=second.mod.est$lambda[k],method="lambda"))[16:29]
        ma.poly=c(1,ma.coef)
        ma.root=polyroot(ma.poly)
        ma.check=sum(abs(ma.root)>1)==length(ma.root)
        invertibility[k]=ma.check
      }    
    }else{
      stationarity=rep(NA,length(second.aic.out))
      invertibility=rep(NA,length(second.aic.out))
      for(k in 1:length(second.aic.out)){
        ar.coef=as.numeric(coef(second.mod.est,s=second.mod.est$lambda[k],method="lambda"))[2:15]
        st.poly=c(1,-ar.coef)
        st.root=polyroot(st.poly)
        st.check=sum(abs(st.root)>1)==length(st.root)
        stationarity[k]=st.check
        
        ma.coef=as.numeric(coef(second.mod.est,s=second.mod.est$lambda[k],method="lambda"))[16:29]
        ma.poly=c(1,ma.coef)
        ma.root=polyroot(ma.poly)
        ma.check=sum(abs(ma.root)>1)==length(ma.root)
        invertibility[k]=ma.check
      }
    }
    
    if(BIC2){
      which.good=which(stationarity & invertibility)
      good.bic=second.bic.out[which.good]
      good.lambdas=second.mod.est$lambda[which.good]
      second.mod.lambda=good.lambdas[which.min(good.bic)]
    }else{
      which.good=which(stationarity & invertibility)
      good.aic=second.aic.out[which.good]
      good.lambdas=second.mod.est$lambda[which.good]
      second.mod.lambda=good.lambdas[which.min(good.aic)]
    } 
    
    final.mod.coef=as.numeric(coef(second.mod.est,s=second.mod.lambda,
                                   method="lambda"))[-1]
    nonzero.select=which(final.mod.coef!=0)
    final.mod.int=as.numeric(coef(second.mod.est,s=second.mod.lambda,
                                  method="lambda"))[1]
    final.mod.s2=sum((second.y-predict(second.mod.est,newx=second.modX,
                                       s=second.mod.lambda,method="lambda"))^2)/(length(second.y)-
                                                                                   sum(final.mod.coef[nonzero.select]!=0)-1)
    
    out=list(final.mod.coef=final.mod.coef, #Final Selection of Coefficients
             final.mod.int=final.mod.int,   #Final Estimated Intercept
             final.mod.s2=final.mod.s2,     #Final Estimated Noise Variance
             nonzero.select=nonzero.select) #Identifies the Nonzero Parameters
  }
  
  #Estimation Via ADENET
  if(Method=="ADENET"){
    #Initial LASSO Weights
    first.mod.est=glmnet(y=first.y,x=first.modX,standardize=T,alpha=1)
    first.mod.RSS=colSums((first.y-predict(first.mod.est,newx=first.modX))^2)
    
    if(BIC1){
      first.bic.out=log(length(first.y))*first.mod.est$df+
        length(first.y)*log(as.vector(first.mod.RSS)/length(first.y))
    }else{
      first.aic.out=2*first.mod.est$df+
        length(first.y)*log(as.vector(first.mod.RSS)/length(first.y))
    }
    
    if(BIC1){
      stationarity=rep(NA,length(first.bic.out))
      invertibility=rep(NA,length(first.bic.out))
      for(k in 1:length(first.bic.out)){
        ar.coef=as.numeric(coef(first.mod.est,s=first.mod.est$lambda[k],method="lambda"))[2:15]
        st.poly=c(1,-ar.coef)
        st.root=polyroot(st.poly)
        st.check=sum(abs(st.root)>1)==length(st.root)
        stationarity[k]=st.check
        
        ma.coef=as.numeric(coef(first.mod.est,s=first.mod.est$lambda[k],method="lambda"))[16:29]
        ma.poly=c(1,ma.coef)
        ma.root=polyroot(ma.poly)
        ma.check=sum(abs(ma.root)>1)==length(ma.root)
        invertibility[k]=ma.check
      }    
    }else{
      stationarity=rep(NA,length(first.aic.out))
      invertibility=rep(NA,length(first.aic.out))
      for(k in 1:length(first.aic.out)){
        ar.coef=as.numeric(coef(first.mod.est,s=first.mod.est$lambda[k],method="lambda"))[2:15]
        st.poly=c(1,-ar.coef)
        st.root=polyroot(st.poly)
        st.check=sum(abs(st.root)>1)==length(st.root)
        stationarity[k]=st.check
        
        ma.coef=as.numeric(coef(first.mod.est,s=first.mod.est$lambda[k],method="lambda"))[16:29]
        ma.poly=c(1,ma.coef)
        ma.root=polyroot(ma.poly)
        ma.check=sum(abs(ma.root)>1)==length(ma.root)
        invertibility[k]=ma.check
      }
    }
    
    if(BIC1){
      which.good=which(stationarity & invertibility)
      good.bic=first.bic.out[which.good]
      good.lambdas=first.mod.est$lambda[which.good]
      first.mod.lambda=good.lambdas[which.min(good.bic)]
    }else{
      which.good=which(stationarity & invertibility)
      good.aic=first.aic.out[which.good]
      good.lambdas=first.mod.est$lambda[which.good]
      first.mod.lambda=good.lambdas[which.min(good.aic)]
    } 
    
    first.mod.alpha=1
    first.mod.lambda=first.mod.lambda
    first.mod.est=glmnet(y=first.y,x=first.modX,standardize=T,
                         alpha=first.mod.alpha,lambda=first.mod.lambda)
    first.mod.coef=as.numeric(coef(first.mod.est))[-1]
    first.mod.mu=as.numeric(coef(first.mod.est))[1]
    
    weights=abs(first.mod.coef+1/length(first.y))^(-eta)
    
    #Update Model Matrix of AR and MA terms Based off Initial Estimation
    if(updateMA){
      update.mod.predict=rep(NA,length(x))
      update.mod.error=rep(0,length(x))
      for(v in (h+max(maxP,maxQ)):Nt){
        update.mod.predict[v]=first.mod.mu+
          x[(v-h):(v-maxP-h+1)]%*%first.mod.coef[1:maxP]+
          update.mod.error[(v-h):(v-maxQ-h+1)]%*%first.mod.coef[-(1:maxP)]
        update.mod.error[v]=x[v]-update.mod.predict[v]
      }
      update.dataQ=foreach(q=1:maxQ,.combine=cbind)%do%{
        lag.func(update.mod.error,k=(q+h-1))
      }
      second.modX=as.matrix(cbind(dataP,
                                  update.dataQ))[-(1:(max(maxP,maxQ)+h-1)),]
      second.y=x[-(1:(max(maxP,maxQ)+h-1))]
    }else{
      second.modX=first.modX
      second.y=first.y
    }
    
    #Final Elastic Net Estimates (Search Through All Lambdas)
    second.cv.out=foreach(a=1:n.alpha,.combine=rbind)%do%{
      second.mod.est=glmnet(y=second.y,x=second.modX,standardize=T,
                            alpha=alpha[a],penalty.factor=weights)
      second.mod.RSS=colSums((second.y-predict(second.mod.est,
                                               newx=second.modX))^2)
      
      if(BIC2){
        second.bic.out=log(length(second.y))*second.mod.est$df+
          length(second.y)*log(as.vector(second.mod.RSS)/length(second.y))
      }else{
        second.aic.out=2*second.mod.est$df+
          length(second.y)*log(as.vector(second.mod.RSS)/length(second.y))
      }
      
      if(BIC2){
        stationarity=rep(NA,length(second.bic.out))
        invertibility=rep(NA,length(second.bic.out))
        for(k in 1:length(second.bic.out)){
          ar.coef=as.numeric(coef(second.mod.est,s=second.mod.est$lambda[k],method="lambda"))[2:15]
          st.poly=c(1,-ar.coef)
          st.root=polyroot(st.poly)
          st.check=sum(abs(st.root)>1)==length(st.root)
          stationarity[k]=st.check
          
          ma.coef=as.numeric(coef(second.mod.est,s=second.mod.est$lambda[k],method="lambda"))[16:29]
          ma.poly=c(1,ma.coef)
          ma.root=polyroot(ma.poly)
          ma.check=sum(abs(ma.root)>1)==length(ma.root)
          invertibility[k]=ma.check
        }    
      }else{
        stationarity=rep(NA,length(second.aic.out))
        invertibility=rep(NA,length(second.aic.out))
        for(k in 1:length(second.aic.out)){
          ar.coef=as.numeric(coef(second.mod.est,s=second.mod.est$lambda[k],method="lambda"))[2:15]
          st.poly=c(1,-ar.coef)
          st.root=polyroot(st.poly)
          st.check=sum(abs(st.root)>1)==length(st.root)
          stationarity[k]=st.check
          
          ma.coef=as.numeric(coef(second.mod.est,s=second.mod.est$lambda[k],method="lambda"))[16:29]
          ma.poly=c(1,ma.coef)
          ma.root=polyroot(ma.poly)
          ma.check=sum(abs(ma.root)>1)==length(ma.root)
          invertibility[k]=ma.check
        }
      }
      
      if(BIC2){
        which.good=which(stationarity & invertibility)
        good.bic=second.bic.out[which.good]
        good.lambdas=second.mod.est$lambda[which.good]
        second.mod.lambda=good.lambdas[which.min(good.bic)]
        result=c(alpha[a],second.mod.lambda,min(second.bic.out))
      }else{
        which.good=which(stationarity & invertibility)
        good.aic=second.aic.out[which.good]
        good.lambdas=second.mod.est$lambda[which.good]
        second.mod.lambda=good.lambdas[which.min(good.aic)]
        result=c(alpha[a],second.mod.lambda,min(second.aic.out))
      } 
      
      
      
      
      
      
      result
    }  
    
    second.mod.alpha=alpha[which.min(second.cv.out[,3])]
    second.mod.lambda=second.cv.out[which.min(second.cv.out[,3]),2]
    second.mod.est=glmnet(y=second.y,x=second.modX,
                          standardize=T,alpha=second.mod.alpha,
                          lambda=second.mod.lambda,penalty.factor=weights)
    second.mod.coef=as.numeric(coef(second.mod.est))[-1]
    second.mod.mu=as.numeric(coef(second.mod.est))[1]
    
    final.mod.coef=second.mod.coef
    nonzero.select=which(final.mod.coef!=0)
    final.mod.int=second.mod.mu
    final.mod.s2=sum((second.y-predict(second.mod.est,
                                       newx=second.modX))^2)/(length(second.y)-
                                                                sum(final.mod.coef[nonzero.select]!=0)-1)
    
    out=list(final.mod.coef=final.mod.coef, #Final Selection of Coefficients
             final.mod.int=final.mod.int,   #Final Estimated Intercept
             final.mod.s2=final.mod.s2,     #Final Estimated Noise Variance
             nonzero.select=nonzero.select) #Identifies the Nonzero Parameters
  }
  
  return(out)
}


MODELS2$adlasso.mod.AIC=adshrink.check.func(x=train.co2,h=1,long.ar.select=F,maxP=14,maxQ=14,updateMA=F,
                                            BIC1=F,BIC2=F,eta=2,Method="ADLASSO")
MODELS2$adenet.mod.AIC=adshrink.check.func(x=train.co2,h=1,long.ar.select=F,maxP=14,maxQ=14,updateMA=F,
                                            BIC1=F,BIC2=F,eta=2,Method="ADENET",alpha=seq(0,1,0.1))


nMODELS2=length(MODELS2)
SELECT=matrix(0,nMODELS2,ncol=28)
COEF=matrix(0,nMODELS2,ncol=28)
RMSFE=rep(NA,nMODELS2)
MASE=rep(NA,nMODELS2)
MFB=rep(NA,nMODELS2)
MDFB=rep(NA,nMODELS2)
FC.LOWER=matrix(NA,sum(VALIDATION.PERIOD),nMODELS2)
FC.MEAN=matrix(NA,sum(VALIDATION.PERIOD),nMODELS2)
FC.UPPER=matrix(NA,sum(VALIDATION.PERIOD),nMODELS2)

for(k in 1:nMODELS2){
  temp.mod=MODELS2[[k]]
  SELECT[k,temp.mod$nonzero.select]=1
  COEF[k,]=temp.mod$final.mod.coef
  
  fc.lower=rep(NA,length(final.co2))
  fc.mean=rep(NA,length(final.co2))
  fc.upper=rep(NA,length(final.co2))
  
  error=rep(0,length(final.co2))
  
  for(j in 30:length(final.co2)){
    fc=as.numeric(final.co2[(j-1):(j-14)]%*%temp.mod$final.mod.coef[1:14]+
                    error[(j-1):(j-14)]%*%temp.mod$final.mod.coef[-(1:14)]) +
      rnorm(100000,mean=temp.mod$final.mod.int,sd=sqrt(temp.mod$final.mod.s2))
    fc.lower[j]=quantile(fc,0.05)
    fc.mean[j]=mean(fc)
    fc.upper[j]=quantile(fc,0.95)
    error[j]=final.co2[j]-fc.mean[j]
  }
  
  RMSFE[k]=sqrt(mean((final.co2[VALIDATION.PERIOD]-fc.mean[VALIDATION.PERIOD])^2))
  MASE[k]=mean(abs((final.co2[VALIDATION.PERIOD]-fc.mean[VALIDATION.PERIOD])/train.mae))
  MFB[k]=mean(final.co2[VALIDATION.PERIOD]-fc.mean[VALIDATION.PERIOD])
  MDFB[k]=(sum((final.co2[VALIDATION.PERIOD]-fc.mean[VALIDATION.PERIOD])>0)-
             sum((final.co2[VALIDATION.PERIOD]-fc.mean[VALIDATION.PERIOD])<0))/sum(VALIDATION.PERIOD)
  
  FC.LOWER[,k]=fc.lower[VALIDATION.PERIOD]
  FC.MEAN[,k]=fc.mean[VALIDATION.PERIOD]
  FC.UPPER[,k]=fc.upper[VALIDATION.PERIOD]
  
}
######################################################################################################################################

