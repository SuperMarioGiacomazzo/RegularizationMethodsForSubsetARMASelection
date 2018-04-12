options(scipen=999)

#Necessary Packages
library(ggplot2)
library(reshape2)
library(gridExtra)
library(cowplot)

#Set Working Directory
setwd("D:/Mario Documents/Graduate School/Github Repositories/RegularizationMethodsForSubsetARMASelection")
source("subsetARMACode.R")

#Length of Usable Observations => Number of Rows in Model Matrix
N=100

#Simulated Data
x=runif(N)

#Assumed Strength of Serial Correlation 
p=4

#Plot for Out of Sample Cross Validation Under Independence
percent.out=0.20
y=rep("Train",N)
y[((1-percent.out)*N +1):N]="Test"
data.oos=data.frame(t=as.factor(1:N),vert=factor("OOS"),select=factor(y,levels=c("Train","Test")))

oos1<-ggplot(data.oos,aes(x=t,y=vert,fill=select))+ 
  geom_tile()+
  theme(panel.grid.minor = element_line(colour="black", size=1)) +
  scale_fill_manual(values=c("maroon","gold"))+
  ylab("")+xlab("")+
  scale_x_discrete(breaks=paste(seq(0,100,by=10))) +
  geom_vline(xintercept=seq(0.5,100.5,1),color="white",size=0.5)+
  labs(title="Full Time Series (T=100)")+
  theme(legend.position=c(50,1),line=element_blank(),legend.title=element_blank(),
        plot.margin = unit(c(0, 1, -0.9, 0), "cm"))
  

#Plot for Out of Sample Cross Validation for Under Dependence
percent.out=0.20
y=rep("Skip",N)
y[1:((1-percent.out)*N-p)]="Train     "
y[((1-percent.out)*N +1):N]="Test     "

data.oos2=data.frame(t=as.factor(1:N),vert=factor("depOOS"),select=factor(y,levels=c("Train     ","Test     ","Skip")))

oos2<-ggplot(data.oos2,aes(x=t,y=vert,fill=select))+ 
  geom_tile()+
  scale_fill_manual(values=c("maroon","gold","black"))+
  ylab("")+ xlab("Time")+
  labs(title="")+
  scale_x_discrete(breaks=paste(seq(0,100,by=10)),position="top") +
  geom_vline(xintercept=seq(0.5,100.5,1),color="white",size=0.5) +
  theme(legend.position="bottom",line=element_blank(),legend.title=element_blank(),
        plot.margin = unit(c(0, 1,1, 0), "cm"),legend.justification="center")

#Combine OOS Plots for Paper
png(filename="oosplots.png",width=600,height=400)
plot_grid(oos1, oos2, align = "v", nrow = 2, rel_heights = c(0.65, 1))
dev.off()



#Plot for 5-Fold CV for Independent Data
K=5
set.seed(5)
kcv5=sample(1:K,N,replace=T)
X=matrix("Train",N,K)
for(k in 1:K){
  X[which(kcv5==k),k]="Test"
}

X=melt(X)
X$Var1=as.factor(X$Var1)
X$Var2=as.factor(X$Var2)
X$value=factor(X$value,levels=c("Train","Test"))

data.5kcv=data.frame(X)
names(data.5kcv)=c("Time","Fold","Select")

kcv5<-ggplot(data.5kcv,aes(x=Time,y=Fold,fill=Select))+
  geom_tile()+
  scale_fill_manual(values=c("maroon","gold"))+
  ylab("Fold")+ xlab("")+
  scale_x_discrete(breaks=paste(seq(0,100,by=10))) +
  labs(title="Full Time Series (T=100)")+
  geom_vline(xintercept=seq(0.5,100.5,1),color="white",size=0.5) +
  geom_hline(yintercept=seq(-0.5,5.5,1),color="white",size=1.5) +
  theme(legend.position=c(50,1),line=element_blank(),legend.title=element_blank(),
        plot.margin = unit(c(0, 1, -0.9, 0), "cm"))
  
#Plot for 10-Fold CV for Independent Data
K=10
set.seed(10)
kcv10=sample(1:K,N,replace=T)
X=matrix("Train     ",N,K)
for(k in 1:K){
  X[which(kcv10==k),k]="Test"
}

X=melt(X)
X$Var1=as.factor(X$Var1)
X$Var2=as.factor(X$Var2)
X$value=factor(X$value,levels=c("Train     ","Test"))

data.10kcv=data.frame(X)
names(data.10kcv)=c("Time","Fold","Select")

kcv10<-ggplot(data.10kcv,aes(x=Time,y=Fold,fill=Select))+
  geom_tile()+
  scale_fill_manual(values=c("maroon","gold"))+
  ylab("Fold")+ xlab("Time")+
  scale_x_discrete(breaks=paste(seq(0,100,by=10)),position="top") +
  labs(title="")+
  geom_vline(xintercept=seq(0.5,100.5,1),color="white",size=0.5) +
  geom_hline(yintercept=seq(-0.5,10.5,1),color="white",size=1.5) +
  theme(legend.position="bottom",line=element_blank(),legend.title=element_blank(),
        plot.margin = unit(c(0, 1,1, 0), "cm"),legend.justification="center")   
  

#Combine kcv plots for Paper
png(filename="kcvplots.png",width=600,height=400)
plot_grid(kcv5, kcv10, align = "v", nrow = 2, rel_heights = c(0.4, 1))
dev.off()




#Plot for 5-Fold CV for Dependent Data
K=5
set.seed(5)
kcv5=sample(1:K,N,replace=T)
X=matrix("Train     ",N,K)
for(k in 1:K){
  loc=which(kcv5==k)
  X[loc,k]="Test     "
  lb=loc-p
  ub=loc+p
  loc2=which(kcv5!=k)
  
  for(v in 1:length(loc)){
    loc3=which(loc2>=lb[v] & loc2<=ub[v])
    X[loc2[loc3],k]="Skip"
  }
}

X=melt(X)
X$Var1=as.factor(X$Var1)
X$Var2=as.factor(X$Var2)
X$value=factor(X$value,levels=c("Train     ","Test     ","Skip"))

data.5kcv=data.frame(X)
names(data.5kcv)=c("Time","Fold","Select")

depkcv5<-ggplot(data.5kcv,aes(x=Time,y=Fold,fill=Select))+
  geom_tile()+
  scale_fill_manual(values=c("maroon","gold","black"))+
  ylab("Fold")+ xlab("")+
  scale_x_discrete(breaks=paste(seq(0,100,by=10))) +
  labs(title="Full Time Series (T=100)")+
  geom_vline(xintercept=seq(0.5,100.5,1),color="white",size=0.5) +
  geom_hline(yintercept=seq(-0.5,5.5,1),color="white",size=1.5) +
  theme(legend.position=c(50,1),line=element_blank(),legend.title=element_blank(),
        plot.margin = unit(c(0, 1, -0.9, 0), "cm"))




#Plot for 10-Fold CV for Dependent Data
K=10
set.seed(10)
kcv10=sample(1:K,N,replace=T)
X=matrix("Train     ",N,K)
for(k in 1:K){
  loc=which(kcv10==k)
  X[loc,k]="Test     "
  lb=loc-p
  ub=loc+p
  loc2=which(kcv10!=k)
  
  for(v in 1:length(loc)){
    loc3=which(loc2>=lb[v] & loc2<=ub[v])
    X[loc2[loc3],k]="Skip"
  }
}

X=melt(X)
X$Var1=as.factor(X$Var1)
X$Var2=as.factor(X$Var2)
X$value=factor(X$value,levels=c("Train     ","Test     ","Skip"))

data.10kcv=data.frame(X)
names(data.10kcv)=c("Time","Fold","Select")

depkcv10<-ggplot(data.10kcv,aes(x=Time,y=Fold,fill=Select))+
  geom_tile()+
  scale_fill_manual(values=c("maroon","gold","black"))+
  ylab("Fold")+ xlab("Time")+
  scale_x_discrete(breaks=paste(seq(0,100,by=10)),position="top") +
  labs(title="")+
  geom_vline(xintercept=seq(0.5,100.5,1),color="white",size=0.5) +
  geom_hline(yintercept=seq(-0.5,10.5,1),color="white",size=1.5) +
  theme(legend.position="bottom",line=element_blank(),legend.title=element_blank(),
        plot.margin = unit(c(0, 1,1, 0), "cm"),legend.justification="center")   

#Combine depkcv plots for Paper
png(filename="depkcvplots.png",width=600,height=400)
plot_grid(depkcv5, depkcv10, align = "v", nrow = 2, rel_heights = c(0.4, 1))
dev.off()




#My Method
set.seed(5)
X.LOOBCV=NonDepCV1.func(x=x,max.pq=p)
X2.LOOBCV=melt(X.LOOBCV)
X2.LOOBCV$Var1=as.factor(X2.LOOBCV$Var1)
X2.LOOBCV$Var2=as.factor(X2.LOOBCV$Var2)
X2.LOOBCV$value[X2.LOOBCV$value==1]="Test     "
X2.LOOBCV$value[X2.LOOBCV$value==0]="Train     "
X2.LOOBCV$value[is.na(X2.LOOBCV$value)]="Skip"
X2.LOOBCV$value=factor(X2.LOOBCV$value,levels=c("Train     ","Test     ","Skip"))
names(X2.LOOBCV)=c("Time","Fold","Select")

loobcv=ggplot(X2.LOOBCV,aes(x=Time,y=Fold,fill=Select))+
  geom_tile()+
  scale_fill_manual(values=c("maroon","gold","black"))+
  ylab("Fold")+ xlab("Time")+
  scale_x_discrete(breaks=paste(seq(0,100,by=4))) +
  scale_y_discrete(breaks=paste(seq(0,25,by=5))) +
  labs(title="Full Time Series (T=100)")+
  geom_vline(xintercept=seq(0.5,100.5,4),color="white",size=1.5) +
  geom_hline(yintercept=seq(-0.5,25.5,1),color="white",size=1.5) +
  theme(legend.position="bottom",line=element_blank(),legend.title=element_blank(),
        plot.margin = unit(c(0, 1,0, 0), "cm"),legend.justification="center")  

png(filename="lobocvplots.png",width=600,height=400)
loobcv
dev.off()




#Bergmeir Method for 5-Blocks and 10-Blocks
set.seed(5)
X.BCV5=NonDepCV2.func(x=x,max.pq=p,K=5)
X2.BCV5=melt(X.BCV5)
X2.BCV5$Var1=as.factor(X2.BCV5$Var1)
X2.BCV5$Var2=as.factor(X2.BCV5$Var2)
X2.BCV5$value[X2.BCV5$value==1]="Test     "
X2.BCV5$value[X2.BCV5$value==0]="Train     "
X2.BCV5$value[is.na(X2.BCV5$value)]="Skip"
X2.BCV5$value=factor(X2.BCV5$value,levels=c("Train     ","Test     ","Skip"))
names(X2.BCV5)=c("Time","Fold","Select")

BCV5=ggplot(X2.BCV5,aes(x=Time,y=Fold,fill=Select))+
  geom_tile()+
  scale_fill_manual(values=c("maroon","gold","black"))+
  ylab("Fold")+ xlab("")+
  scale_x_discrete(breaks=paste(seq(0,100,by=20))) +
  scale_y_discrete(breaks=paste(seq(1,5,by=1))) +
  labs(title="Full Time Series (T=100)")+
  geom_vline(xintercept=seq(0.5,100.5,20),color="white",size=1.5) +
  geom_hline(yintercept=seq(-0.5,25.5,1),color="white",size=1.5) +
  theme(legend.position="none",line=element_blank(),legend.title=element_blank(),
        plot.margin = unit(c(0, 1, -1, 0), "cm"))  

set.seed(5)
X.BCV10=NonDepCV2.func(x=x,max.pq=p,K=10)
X2.BCV10=melt(X.BCV10)
X2.BCV10$Var1=as.factor(X2.BCV10$Var1)
X2.BCV10$Var2=as.factor(X2.BCV10$Var2)
X2.BCV10$value[X2.BCV10$value==1]="Test     "
X2.BCV10$value[X2.BCV10$value==0]="Train     "
X2.BCV10$value[is.na(X2.BCV10$value)]="Skip"
X2.BCV10$value=factor(X2.BCV10$value,levels=c("Train     ","Test     ","Skip"))
names(X2.BCV10)=c("Time","Fold","Select")

BCV10=ggplot(X2.BCV10,aes(x=Time,y=Fold,fill=Select))+
  geom_tile()+
  scale_fill_manual(values=c("maroon","gold","black"))+
  ylab("Fold")+ xlab("Time")+
  scale_x_discrete(breaks=paste(seq(0,100,by=10)),position="top") +
  scale_y_discrete(breaks=paste(seq(1,10,by=1))) +
  labs(title="")+
  geom_vline(xintercept=seq(0.5,100.5,10),color="white",size=1.5) +
  geom_hline(yintercept=seq(-0.5,25.5,1),color="white",size=1.5) +
  theme(legend.position="bottom",line=element_blank(),legend.title=element_blank(),
        plot.margin = unit(c(0, 1,1, 0), "cm"),legend.justification="center")  

png(filename="bcvplots.png",width=600,height=400)
plot_grid(BCV5, BCV10, align = "v", nrow = 2, rel_heights = c(0.4, 1))
dev.off()









