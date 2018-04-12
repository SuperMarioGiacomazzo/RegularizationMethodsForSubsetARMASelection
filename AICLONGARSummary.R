options(scipen=999)

#Packages Required
library(xtable)

#Set Working Directory
setwd("D:/Mario Documents/Graduate School/Github Repositories/RegularizationMethodsForSubsetARMASelection")

source("subsetARMACode.R")
load("ADLASSOSIM1AICBIC.Rdata")
coef.true1=coef.true
RESULTS1.LONG=RESULTS[-1,]
load("ADLASSOSIM2AICBIC.Rdata")
coef.true2=coef.true
RESULTS2.LONG=RESULTS[-1,]
load("ADLASSOSIM3AICBIC.Rdata")
coef.true3=coef.true
RESULTS3.LONG=RESULTS[-1,]
load("ADLASSOSIM1AICBICB.Rdata")
RESULTS1.SHORT=RESULTS[-1,]
load("ADLASSOSIM2AICBICB.Rdata")
RESULTS2.SHORT=RESULTS[-1,]
load("ADLASSOSIM3AICBICB.Rdata")
RESULTS3.SHORT=RESULTS[-1,]

#Vector of Column Names
colheads=names(RESULTS1.LONG)

#Function to Create Tables
create.table.func<-function(dataset1,dataset2){
  AIC1=subset(dataset1,Method=="ADLASSO AIC")
  AICBIC1=subset(dataset1,Method=="ADLASSO AIC/BIC")
  BIC1=subset(dataset1,Method=="ADLASSO BIC")
  
  AIC2=subset(dataset2,Method=="ADLASSO AIC")
  AICBIC2=subset(dataset2,Method=="ADLASSO AIC/BIC")
  BIC2=subset(dataset2,Method=="ADLASSO BIC")
  
  time.spec=unique(dataset1$Length)
  n.time=length(time.spec)
  noise.spec=unique(dataset1$Noise)
  n.noise=length(noise.spec)
  
  FINAL1=NULL
  for(j in 1:n.time){
      sumAIC1=colMeans(subset(AIC1,Length==time.spec[j] & Noise==noise.spec[2])[,colheads[33:36]])
      sumAIC2=colMeans(subset(AIC2,Length==time.spec[j] & Noise==noise.spec[2])[,colheads[33:36]])
      FINAL1=rbind(FINAL1,c(time.spec[j],sumAIC1,sumAIC2))
  }
  
  FINAL2=NULL
  for(j in 1:n.time){
    sumAICBIC1=colMeans(subset(AICBIC1,Length==time.spec[j] & Noise==noise.spec[2])[,colheads[33:36]])
    sumAICBIC2=colMeans(subset(AICBIC2,Length==time.spec[j] & Noise==noise.spec[2])[,colheads[33:36]])
    FINAL2=rbind(FINAL2,c(time.spec[j],sumAICBIC1,sumAICBIC2))
  }
  
  FINAL3=NULL
  for(j in 1:n.time){
    sumBIC1=colMeans(subset(BIC1,Length==time.spec[j] & Noise==noise.spec[2])[,colheads[33:36]])
    sumBIC2=colMeans(subset(BIC2,Length==time.spec[j] & Noise==noise.spec[2])[,colheads[33:36]])
    FINAL3=rbind(FINAL3,c(time.spec[j],sumBIC1,sumBIC2))
  }
      
  FINAL=round(rbind(FINAL1,FINAL2,FINAL3),2)
  
  rownames(FINAL)=NULL
  colnames(FINAL)=NULL
  
  return(FINAL)
}

#Creation of Tables for Simulation 1
sim1.tab=create.table.func(RESULTS1.LONG,RESULTS1.SHORT)
sim1.tab.latex=xtable(sim1.tab,include.rownames=F)
align(sim1.tab.latex)=rep("c",10)
print(sim1.tab.latex,include.rownames=F,include.colnames=F)

#Creation of Tables for Simulation 2
sim2.tab=create.table.func(RESULTS2.LONG,RESULTS2.SHORT)
sim2.tab.latex=xtable(sim2.tab,include.rownames=F)
align(sim2.tab.latex)=rep("c",10)
print(sim2.tab.latex,include.rownames=F,include.colnames=F)

#Creation of Tables for Simulation 3
sim3.tab=create.table.func(RESULTS3.LONG,RESULTS3.SHORT)
sim3.tab.latex=xtable(sim3.tab,include.rownames=F)
align(sim3.tab.latex)=rep("c",10)
print(sim3.tab.latex,include.rownames=F,include.colnames=F)





















