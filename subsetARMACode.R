##############################
#Function to lag a time series
#Arguments: x=time series
#           k=lag (Defaults to 1)
##############################

lag.func<-function(x,k=1){
  t=length(x)
  y=c(rep(NA,t))
  for(i in (k+1):t){
    y[i]=x[i-k]
  }
  return(y)
}
##############################




########################################################################
#Function to Obtain Projected Posterior Distribution from Full Posterior
#Arguments: fullpost=Posterior Distribution from BHS estimated model
#           X=Model Matrix of Full Model
#           indproj=Vector Indicating which Columns are Included
#Source: Piironen and Vehtari (2015)
########################################################################

pms.proj.func<-function(fullpost,X,indproj){
  lres.coef=t(fullpost[,-dim(fullpost)[2]])   #Transpose Matrix of Posterior Samples of ARMA Coefficients
  lres.s2=(fullpost[,dim(fullpost)[2]])  #Vector of Posterior Samples of Variance Parameter
  S=length(lres.s2)   #Number of Posterior Samples Obtained
  N=dim(X)[1]   #Number of Points in Data Set
  
  P=dim(X)[2]
  COEF.PROJ=matrix(0,S,P)
  X.proj=X[,indproj]
  pred.proj=X%*%lres.coef
  
  coef.proj=solve(t(X.proj)%*%X.proj)%*%t(X.proj)%*%pred.proj
  var.proj=c(lres.s2)+colMeans((pred.proj-X.proj%*%coef.proj)^2)
  KL.PROJ=0.5*log(var.proj/lres.s2)
  COEF.PROJ[,indproj]=t(coef.proj)
  KL.MEAN=mean(KL.PROJ)
  COEF.MEAN=colMeans(COEF.PROJ)
  VAR.PROJ=var.proj
  VAR.MEAN=mean(var.proj)
  return(list(KL.MEAN=KL.MEAN,COEF.MEAN=COEF.MEAN,COEF.PROJ=COEF.PROJ,VAR.MEAN=VAR.MEAN,VAR.PROJ=VAR.PROJ))
} 
########################################################################




#############################################################################################
#Function to Conduct ARMA Selection via Bayesian Projection Posterior Predictive Distribution
#   Implementing Relative Efficiency for Final Model Selection
#Arguments: x=Time Series to Be Modeled Using ARMA Process
#           h=Horizon Specific Model (Defaults to 1)
#           cores=Number of Cores for Parallel Processing (Defaults to 3)
#           maxP=Maximum Autoregressive Order 
#           maxQ=Maximum Moving Average Order
#           KL.threshold=Single value or vector of chosen thresholds for stopping rule
#                        based on Relative Efficiency Comparing Submodel to Full Model
#                        based on Kullback Leibler Divergence
#                        (Defaults to c(0.9,0.95,0.99))
#Source: Piironen and Vehtari (2015)
#############################################################################################

pms.arma.func<-function(x,h=1,cores=1,maxP,maxQ,KL.threshold=c(0.90,0.95,0.99),prior.choice=c("hs","hs+"),updateMA=F){
  require(MCMCpack)
  require(bayesreg)
  require(doParallel)
  require(foreach)
  
  prior.choice=match.arg(prior.choice)
  
  N=length(x)
  max.ar.order=ceiling(10*log10(N))
  
  cl<-makeCluster(min(cores,max.ar.order))
  registerDoParallel(cl)
  init.modX=foreach(init.ar=1:max.ar.order,.combine=cbind,.export="lag.func")%dopar%{
    lag.func(x,k=(init.ar+h-1))
  }
  stopCluster(cl)
  
  init.data=data.frame(y=x,init.modX)
  init.data=init.data[-(1:(max.ar.order+h-1)),]
  
  init.mod.est=MCMCregress(y~.,data=init.data,mcmc=2000,thin=10,burnin=10000)
  muBeta0=mean(init.mod.est[,1])
  muBeta=colMeans(init.mod.est[,-c(1,dim(init.mod.est)[2])])
  init.mod.error=init.data$y-(as.numeric(muBeta0) + as.matrix(init.data[,-1])%*%as.vector(muBeta))
  
  cl<-makeCluster(min(cores,maxP))
  registerDoParallel(cl)
  dataP=foreach(p=1:maxP,.combine=cbind,.export="lag.func")%dopar%{
    lag.func(init.data$y,k=(p+h-1))
  }
  stopCluster(cl)
  
  cl<-makeCluster(min(cores,maxQ))
  registerDoParallel(cl)
  dataQ=foreach(q=1:maxQ,.combine=cbind,.export="lag.func")%dopar%{
    lag.func(init.mod.error,k=(q+h-1))
  }
  stopCluster(cl)
  
  full.data=data.frame(y=init.data$y,dataP=dataP,dataQ=dataQ)
  full.data=full.data[-(1:(max(maxP,maxQ)+h-1)),]
  
  full.mod.est=bayesreg(y~.,data=full.data,prior=prior.choice,nsamples=2000,thin=10,burnin=10000)  
  full.mod.posterior=cbind(as.vector(full.mod.est$beta0),t(as.matrix(full.mod.est$beta)),as.vector(full.mod.est$sigma2))
  
  full.mod.int=c(full.mod.est$muBeta0)
  full.mod.coef=c(full.mod.est$muBeta)
  full.mod.s2=full.mod.est$muSigma2
  
  if(updateMA){
    full.mod.predict=rep(NA,length(x))
    full.mod.error=rep(0,length(x))
    for(k in (h+max(maxP,maxQ)):N){
      full.mod.predict[k]=full.mod.int+x[(k-h):(k-maxP-h+1)]%*%full.mod.coef[1:maxP]+full.mod.error[(k-h):(k-maxQ-h+1)]%*%full.mod.coef[-(1:maxP)]
      full.mod.error[k]=x[k]-full.mod.predict[k]
    }
    
    cl<-makeCluster(min(cores,maxP))
    registerDoParallel(cl)
    full.dataP=foreach(p=1:maxP,.combine=cbind,.export="lag.func")%dopar%{
      lag.func(x,k=(p+h-1))
    }
    stopCluster(cl)
    
    cl<-makeCluster(min(cores,maxQ))
    registerDoParallel(cl)
    update.dataQ=foreach(p=1:maxQ,.combine=cbind,.export="lag.func")%dopar%{
      lag.func(full.mod.error,k=(p+h-1))
    }
    stopCluster(cl)
  
    full.mod.X=cbind(1,full.dataP,update.dataQ)
    full.mod.X=full.mod.X[-(1:(max(maxP,maxQ)+h-1)),]
    nMod=dim(full.mod.X)[2]
  }else{
    full.mod.X=as.matrix(cbind(1,full.data[,-1]))
    nMod=dim(full.mod.X)[2]
  }
    
    
  KL=rep(NA,nMod)
  chosen=1
  notchosen=setdiff(1:nMod,chosen)
  FIRST=pms.proj.func(fullpost=full.mod.posterior,X=full.mod.X,indproj=chosen)
  KL[1]=FIRST$KL.MEAN
  
  for(modnum in 2:nMod){
    nleft<-length(notchosen)
    cl<-makeCluster(cores,outfile="debug.txt")
    registerDoParallel(cl)
    val=foreach(j=1:nleft,.combine=c,.export="pms.proj.func")%dopar%{
      ind<-sort(c(chosen,notchosen[j]))
      NEXT<-tryCatch({pms.proj.func(fullpost=full.mod.posterior,X=full.mod.X,indproj=ind)$KL.MEAN},
                     error=function(e){return(NA)})
      NEXT
    }
    stopCluster(cl)
    minval<-which.min(val)
    chosen<-c(chosen,notchosen[minval])
    notchosen<-setdiff(1:nMod,chosen)
    KL[modnum]<-val[minval]
  }
  
  nKL.threshold=length(KL.threshold)
  KL.select=matrix(rep(chosen,nKL.threshold),ncol=nKL.threshold)
  final.mod.posterior=list()
  final.mod.int=list()
  final.mod.coef=list()
  final.mod.s2=list()
  
  for(k in 1:nKL.threshold){
    TEMP1=min(which((1-KL/KL[1])>KL.threshold[k]))
    TEMP2=pms.proj.func(fullpost=full.mod.posterior,X=full.mod.X,indproj=chosen[1:TEMP1])
    KL.select[-(1:TEMP1),k]=NA
    final.mod.posterior[[k]]=cbind(TEMP2$COEF.PROJ,TEMP2$VAR.PROJ)
    final.mod.int[[k]]=TEMP2$COEF.MEAN[1]
    final.mod.coef[[k]]=TEMP2$COEF.MEAN[-1]
    final.mod.s2[[k]]=TEMP2$VAR.MEAN
  }
  
  return(list(CHOICE=chosen,KL=KL,KL.threshold=KL.threshold,KL.select=KL.select,
              full.mod.int=full.mod.int,full.mod.coef=full.mod.coef,full.mod.s2=full.mod.s2,
              final.mod.int=final.mod.int,final.mod.coef=final.mod.coef,final.mod.s2=final.mod.s2))
}    
#############################################################################################






#############################################################################################
#Functions to Evaluate Subset Model Selection of Various Methods
#Evaluation Measures: id.sig.coef.func=Identifies if Non-zero Parameters Selected
#                     id.true.coef.func=Identifies if True Model Chosen
#                     false.neg.func=Proportion of zero parameters identified as Non-zero
#                     false.pos.func=Proportion of Non-zero parameters selected
#                     fulleval.func=Measures all of the Above
#Arguments: truecoef=true known coefficients
#           estcoef=Estimated Coefficients
#############################################################################################

id.sig.coef.func<-function(truecoef,estcoef){
  all(which(truecoef!=0) %in% which(estcoef!=0))
}
id.true.coef.func<-function(truecoef,estcoef){
  if(length(which(truecoef!=0))==length(which(estcoef!=0))){
    all(sort(which(truecoef!=0))==as.numeric(sort(which(estcoef!=0))))
  }else{
    FALSE
  }
}
false.neg.func<-function(truecoef,estcoef){
  mean(estcoef[which(truecoef!=0)]==0)
}
false.pos.func<-function(truecoef,estcoef){
  mean(estcoef[which(truecoef==0)]!=0)
}
fulleval.func<-function(truecoef,estcoef){
  return(c(id.sig.coef.func(truecoef,estcoef),id.true.coef.func(truecoef,estcoef),
           false.neg.func(truecoef,estcoef),false.pos.func(truecoef,estcoef)))
}
#############################################################################################





#############################################################################################
#Functions Used  to Perform K-fold Non-Dependent Cross Validation. All Three methods are
#   based on dividing time series into blocks and performing cross validation with the blocks
#   removing data that has mutual dependence with test and training sets
#Arguments: x = time series
#           max.pq = maximum ar/ma order to consider (covers temporal dependence)
#           K = Number of Folds to Consider
#############################################################################################

#Method 1: My Method 1 (Partition Data Into Blocks of Length max.pq)
NonDepCV1.func<-function(x,max.pq=NULL,K=NULL){
  N=length(x)
  if(is.null(max.pq)) max.pq=ceiling(10*log10(N))
  if(max.pq>ceiling(10*log10(N))) warning("Maximum ARMA Order Considered Too Large Based on Time Series Length")
  nblocks=floor(N/max.pq)
  if(is.null(K)) K=nblocks
  blocks=rep(nblocks,N)
  for(b in 1:(nblocks-1)){
    blocks[(b*max.pq-max.pq+1):(b*max.pq)]=b
  }
  if(K>nblocks) warning("Maximum Number of Folds Surpassed Based on Choice of max.pq and Time Series Length")
  test.blocks=sort(sample(1:nblocks,K))
  Block.Matrix=matrix(0,N,K) 
  for(v in 1:K){
    in.test=which(blocks==test.blocks[v]) 
    Block.Matrix[in.test,v]=1 
    if(test.blocks[v]==1){
      out.ignore=which(blocks==2)
      Block.Matrix[out.ignore,v]=NA
    }else if(test.blocks[v]==nblocks){
      out.ignore=which(blocks==(nblocks-1))
      Block.Matrix[out.ignore,v]=NA
    }else{
      out.ignore=which(blocks %in% c(test.blocks[v]-1,test.blocks[v]+1))
      Block.Matrix[out.ignore,v]=NA
    }
  }
  return(Block.Matrix)
} 

#Method 2: Begmeir Method (Divide Data Into K Blocks and Remove Boundary Dependent Data)
NonDepCV2.func<-function(x,max.pq=NULL,K=NULL){
  N=length(x)
  obs.per.block=floor(N/K)
  if(is.null(max.pq)) max.pq=ceiling(10*log10(obs.per.block))
  if((max.pq/2)>obs.per.block) warning("Choice of max.pq and K Lead to No Observations")
  
  blocks=rep(NA,N)
  for(b in 1:K){
    blocks[(b*obs.per.block-obs.per.block+ceiling(max.pq/2)+1):(b*obs.per.block-ceiling(max.pq)/2)]=b
  }
  
  Block.Matrix=matrix(NA,N,K) 
  for(v in 1:K){
    in.test=which(blocks==(1:K)[v]) 
    in.train=which(blocks %in% (1:K)[-v])
    Block.Matrix[in.test,v]=1 
    Block.Matrix[in.train,v]=0
  }
  
  return(Block.Matrix)
}    

#Method 3: Always Forecast Forward Method
NonDepCV3.func<-function(x,max.pq=NULL,K=NULL){
  N=length(x)
  obs.per.block=floor(N/K)
  if(is.null(max.pq)) max.pq=ceiling(10*log10(obs.per.block))
  if(max.pq>(obs.per.block)) warning("Choice of max.pq and K Lead to No Observations")
  
  Block.Matrix=matrix(NA,N,K) 
  for(v in 1:K){
    Block.Matrix[(1:(obs.per.block*k-max.pq)),v]=0
    Block.Matrix[-(1:(obs.per.block*k)),v]=1
  }
  
  return(Block.Matrix)
}
#############################################################################################




#############################################################################################
#Functions Used  to Partition Data into a Train Set containing the first (1-test.per)x100%
#   of data and a Test Set containing the last (test.per)x100% of data. The first function 
#   institutes a gap of length max.pq to remove dependencies between train and test set
#Arguments: x = time series
#           max.pq = maximum ar/ma order to consider (covers temporal dependence)
#           test.per = Percent of Data to Consider in Test Set
#############################################################################################

OOS.DepCV.func<-function(x,max.pq=NULL,test.per=0.20){
  N=length(x)
  test.N=ceiling(test.per*N)
  if(is.null(max.pq)) max.pq=ceiling(10*log10(test.N))
  if(test.per>0.5) warning("Less than Half the Dataset Is Being Used for Training ")
  if(max.pq>(N-test.N)) warning("max.pq > Number of Obs in Training")
  Block.Vector=rep(0,N)
  Block.Vector[(N-test.N+1):N]=1
  Block.Vector[(N-test.N-max.pq+1):(N-test.N)]=NA
  return(as.matrix(Block.Vector))
}

OOS.IndepCV.func<-function(x,test.per=0.20){
  N=length(x)
  test.N=ceiling(test.per*N)
  if(test.per>0.5) warning("Less than Half the Dataset Is Being Used for Training ")
  Block.Vector=rep(0,N)
  Block.Vector[(N-test.N+1):N]=1
  return(as.matrix(Block.Vector))
}
#############################################################################################





#############################################################################################
#Function to Conduct ARMA Selection via Adaptive Shrinkage Methods With Cross Validation
#Arguments: x=Time Series to Be Modeled Using ARMA Process
#           h=Horizon Specific Model (Defaults to 1)
#           cores=Number of Cores for Parallel Processing (Defaults to 3)
#           long.ar.select=Indicator Determining if Model Selection Should Be Performed in 
#               the Initial Modeling of Long AR Process
#           maxP=Maximum Autoregressive Order 
#           maxQ=Maximum Moving Average Order
#           update=Should Moving Average Terms Be Updated After Initial Coefficients Selected
#                   (Default to True)
#           BIC1=Should BIC be used to choose shrinkage parameter in initial lasso step
#           BIC2=Should BIC be used to choose shrinkage parameter in adaptive lasso step
#           alpha=Elastic Net Mixing Parameter (0=Ridge, 1=Lasso, Other=Elastic Net)
#           eta=Exponent Applied to Weights for Adaptive LASSO
#           test.per = percent of data left out for OOS Test Methods (Default 20%)
#           max.pq = Assumed Maximum ARMA order to Remove Dependence in CV methods
#           K = Number of Folds to Consider in CV methods
#           CV.method= Choose 1 Method of Cross Validation   (5 Methods)  
#                      (AIC/BIC) -> No Cross Validation AIC and/or BIC is used for selection
#                                   of tuning parameters for both stages of adaptive lasso
#Source: Wang and Leng (2007) and Efron et al. (2004) and Zou and Hastie (2005)
#############################################################################################

#Performs Adaptive Lasso
cv.adaptlasso.arma.func<-function(x,h=1,cores=1,long.ar.select=F,maxP,maxQ,max.pq=NULL,K=NULL,updateMA=F,
    test.per=0.2,BIC1=F,BIC2=F,eta=2,CV.method=c("AIC/BIC","LOBOCV","Bergmeir","RegCV","OOSINDEP","OOSDEP")){
  
  #########################
  #Require Certain Packages
  #########################
  require(glmnet) #Performs Ridge, Lasso, Elastic Net
  require(lars)   #Use to Get Lambda Sequence
  require(forecast)
  require(doParallel)
  require(foreach)
  
  CV.method=match.arg(CV.method)
  
  Nt=length(x)
  max.ar.order=ceiling(10*log10(Nt))
  
  init.mod.est=ar(x,aic=long.ar.select,order.max=max.ar.order,demean=T)
  init.mod.error=residuals(init.mod.est)
  init.mod.order=length(which(is.na(init.mod.error)))
  
  cl<-makeCluster(min(cores,maxP))
  registerDoParallel(cl)
  dataP=foreach(p=1:maxP,.combine=cbind,.export="lag.func")%dopar%{
    lag.func(x,k=(p+h-1))
  }
  stopCluster(cl)
  
  cl<-makeCluster(min(cores,maxQ))
  registerDoParallel(cl)
  dataQ=foreach(q=1:maxQ,.combine=cbind,.export="lag.func")%dopar%{
    lag.func(init.mod.error,k=(q+h-1))
  }
  stopCluster(cl)
  
  first.modX=as.matrix(cbind(dataP,dataQ))[-(1:(init.mod.order+max(maxP,maxQ)+h-1)),]
  first.y=x[-(1:(init.mod.order+max(maxP,maxQ)+h-1))]
  
  if(CV.method=="AIC/BIC"){
    first.mod.est=glmnet(y=first.y,x=first.modX,standardize=T,alpha=1)
    first.mod.RSS=colSums((first.y-predict(first.mod.est,newx=first.modX))^2)
    if(BIC1){
      first.bic.out=log(length(first.y))*first.mod.est$df+length(first.y)*log(as.vector(first.mod.RSS)/length(first.y))
      first.mod.lambda=first.mod.est$lambda[which.min(first.bic.out)]
    }else{
      first.aic.out=2*first.mod.est$df+length(first.y)*log(as.vector(first.mod.RSS)/length(first.y))
      first.mod.lambda=first.mod.est$lambda[which.min(first.aic.out)]
    }
    first.mod.coef=as.numeric(coef(first.mod.est,s=first.mod.lambda,method="lambda"))[-1]
    first.mod.mu=as.numeric(coef(first.mod.est,s=first.mod.lambda,method="lambda"))[1]
    weights=abs(first.mod.coef+1/length(first.y))^(-eta)

    if(updateMA){
      update.mod.predict=rep(NA,length(x))
      update.mod.error=rep(0,length(x))
      for(v in (h+max(maxP,maxQ)):Nt){
        update.mod.predict[v]=first.mod.mu+x[(v-h):(v-maxP-h+1)]%*%first.mod.coef[1:maxP]+
          update.mod.error[(v-h):(v-maxQ-h+1)]%*%first.mod.coef[-(1:maxP)]
        update.mod.error[v]=x[v]-update.mod.predict[v]
      }
      cl<-makeCluster(min(cores,maxQ))
      registerDoParallel(cl)
      update.dataQ=foreach(q=1:maxQ,.combine=cbind,.export="lag.func")%dopar%{
        lag.func(update.mod.error,k=(q+h-1))
      }
      stopCluster(cl)
      second.modX=as.matrix(cbind(dataP,update.dataQ))[-(1:(max(maxP,maxQ)+h-1)),]
      second.y=x[-(1:(max(maxP,maxQ)+h-1))]
    }else{
      second.modX=first.modX
      second.y=first.y
    }
    
    second.mod.est=glmnet(y=second.y,x=second.modX,standardize=T,alpha=1,thresh=1e-16,penalty.factor=weights)
    second.mod.RSS=colSums((second.y-predict(second.mod.est,newx=second.modX))^2)
    if(BIC2){
      second.bic.out=log(length(second.y))*second.mod.est$df+length(second.y)*log(as.vector(second.mod.RSS)/length(second.y))
      second.mod.lambda=second.mod.est$lambda[which.min(second.bic.out)]
    }else{
      second.aic.out=2*second.mod.est$df+length(second.y)*log(as.vector(second.mod.RSS)/length(second.y))
      second.mod.lambda=second.mod.est$lambda[which.min(second.aic.out)]
    }
    
    final.mod.coef=as.numeric(coef(second.mod.est,s=second.mod.lambda,method="lambda"))[-1]
    nonzero.select=which(final.mod.coef!=0)
    final.mod.int=as.numeric(coef(second.mod.est,s=second.mod.lambda,method="lambda"))[1]
    final.mod.s2=sum((second.y-predict(second.mod.est,newx=second.modX,s=second.mod.lambda,method="lambda"))^2)/
      (length(second.y)-sum(final.mod.coef[nonzero.select]!=0)-1)
    
    out=list(final.mod.coef=final.mod.coef,final.mod.int=final.mod.int,final.mod.s2=final.mod.s2,nonzero.select=nonzero.select)
  }
  
  if(CV.method=="LOBOCV"){
    Fold.Matrix=NonDepCV1.func(first.y,max.pq=max.pq,K=K)
    nfolds=dim(Fold.Matrix)[2]
    
    lambda1.seq=cv.glmnet(y=first.y,x=first.modX,standardize=T,alpha=1)$lambda
    
    cl<-makeCluster(cores,outfile="debug.txt")
    registerDoParallel(cl)
    SQDEV1=foreach(f=1:nfolds,.packages=c("glmnet"),.combine=rbind)%dopar%{
          in.train=which(Fold.Matrix[,f]==0)
          in.test=which(Fold.Matrix[,f]==1)
          first.mod.est=glmnet(y=first.y[in.train],x=first.modX[in.train,],standardize=T,alpha=1,lambda=lambda1.seq)
          first.mod.res=(first.y[in.test]-predict(first.mod.est,newx=first.modX[in.test,]))^2
          first.mod.res
    }
    stopCluster(cl)
    
    CVM1=apply(SQDEV1,2,mean)
    lambda1.min=lambda1.seq[which.min(CVM1)]
    lambda1.1se=lambda1.seq[min(which(CVM1<(min(CVM1)+sd(CVM1)/sqrt(length(CVM1)))))]
   
    first.mod.est=glmnet(y=first.y,x=first.modX,standardize=T,alpha=1,lambda=lambda1.min)
    first.mod.coef=as.numeric(coef(first.mod.est))[-1]
    first.mod.mu=as.numeric(coef(first.mod.est))[1]
    weights=abs(first.mod.coef+1/length(first.y))^(-eta)
    
    update.mod.predict=rep(NA,length(x))
    update.mod.error=rep(0,length(x))
    for(v in (h+max(maxP,maxQ)):Nt){
      update.mod.predict[v]=first.mod.mu+x[(v-h):(v-maxP-h+1)]%*%first.mod.coef[1:maxP]+
        update.mod.error[(v-h):(v-maxQ-h+1)]%*%first.mod.coef[-(1:maxP)]
      update.mod.error[v]=x[v]-update.mod.predict[v]
    }
    cl<-makeCluster(min(cores,maxQ))
    registerDoParallel(cl)
    update.dataQ=foreach(q=1:maxQ,.combine=cbind,.export="lag.func")%dopar%{
      lag.func(update.mod.error,k=(q+h-1))
    }
    stopCluster(cl)
    
    if(updateMA){
      update.mod.predict=rep(NA,length(x))
      update.mod.error=rep(0,length(x))
      for(v in (h+max(maxP,maxQ)):Nt){
        update.mod.predict[v]=first.mod.mu+x[(v-h):(v-maxP-h+1)]%*%first.mod.coef[1:maxP]+
          update.mod.error[(v-h):(v-maxQ-h+1)]%*%first.mod.coef[-(1:maxP)]
        update.mod.error[v]=x[v]-update.mod.predict[v]
      }
      cl<-makeCluster(min(cores,maxQ))
      registerDoParallel(cl)
      update.dataQ=foreach(q=1:maxQ,.combine=cbind,.export="lag.func")%dopar%{
        lag.func(update.mod.error,k=(q+h-1))
      }
      stopCluster(cl)
      second.modX=as.matrix(cbind(dataP,update.dataQ))[-(1:(max(maxP,maxQ)+h-1)),]
      second.y=x[-(1:(max(maxP,maxQ)+h-1))]
    }else{
      second.modX=first.modX
      second.y=first.y
    }
    
    lambda2.seq=cv.glmnet(y=second.y,x=second.modX,standardize=T,alpha=1,penalty.factor=weights)$lambda
    
    cl<-makeCluster(cores,outfile="debug.txt")
    registerDoParallel(cl)
    SQDEV2=foreach(f=1:nfolds,.packages=c("glmnet"),.combine=rbind)%dopar%{
      in.train=which(Fold.Matrix[,f]==0)
      in.test=which(Fold.Matrix[,f]==1)
      first.mod.est=glmnet(y=second.y[in.train],x=second.modX[in.train,],standardize=T,alpha=1,lambda=lambda2.seq,penalty.factor=weights)
      first.mod.res=(second.y[in.test]-predict(first.mod.est,newx=second.modX[in.test,]))^2
      first.mod.res
    }
    stopCluster(cl)
    
    CVM2=apply(SQDEV2,2,mean)
    lambda2.min=lambda2.seq[which.min(CVM2)]
    lambda2.1se=lambda2.seq[min(which(CVM2<(min(CVM2)+sd(CVM2)/sqrt(length(CVM2)))))]
    
    final.mod.est=glmnet(y=second.y,x=second.modX,standardize=T,alpha=1,lambda=lambda2.1se,penalty.factor=weights)
    final.mod.coef=as.numeric(coef(final.mod.est))[-1]
    nonzero.select=which(final.mod.coef!=0)
    final.mod.int=as.numeric(coef(final.mod.est))[1]
    final.mod.s2=sum((second.y-predict(final.mod.est,newx=second.modX,s=lambda2.1se,method="lambda"))^2)/
      (length(second.y)-sum(final.mod.coef[nonzero.select]!=0)-1)
    
    out=list(final.mod.coef=final.mod.coef,final.mod.int=final.mod.int,final.mod.s2=final.mod.s2,nonzero.select=nonzero.select)
  }
        
  if(CV.method=="Bergmeir"){
    Fold.Matrix=NonDepCV2.func(x=first.y,max.pq=max.pq,K=K)
    nfolds=dim(Fold.Matrix)[2]
    
    lambda1.seq=cv.glmnet(y=first.y,x=first.modX,standardize=T,alpha=1)$lambda
    
    cl<-makeCluster(cores,outfile="debug.txt")
    registerDoParallel(cl)
    SQDEV1=foreach(f=1:nfolds,.packages=c("glmnet"),.combine=rbind)%dopar%{
      in.train=which(Fold.Matrix[,f]==0)
      in.test=which(Fold.Matrix[,f]==1)
      first.mod.est=glmnet(y=first.y[in.train],x=first.modX[in.train,],standardize=T,alpha=1,lambda=lambda1.seq)
      first.mod.res=(first.y[in.test]-predict(first.mod.est,newx=first.modX[in.test,]))^2
      first.mod.res
    }
    stopCluster(cl)
    
    CVM1=apply(SQDEV1,2,mean)
    lambda1.min=lambda1.seq[which.min(CVM1)]
    lambda1.1se=lambda1.seq[min(which(CVM1<(min(CVM1)+sd(CVM1)/sqrt(length(CVM1)))))]
    
    first.mod.est=glmnet(y=first.y,x=first.modX,standardize=T,alpha=1,lambda=lambda1.min)
    first.mod.coef=as.numeric(coef(first.mod.est))[-1]
    first.mod.mu=as.numeric(coef(first.mod.est))[1]
    weights=abs(first.mod.coef+1/length(first.y))^(-eta)
    
    if(updateMA){
      update.mod.predict=rep(NA,length(x))
      update.mod.error=rep(0,length(x))
      for(v in (h+max(maxP,maxQ)):Nt){
        update.mod.predict[v]=first.mod.mu+x[(v-h):(v-maxP-h+1)]%*%first.mod.coef[1:maxP]+
          update.mod.error[(v-h):(v-maxQ-h+1)]%*%first.mod.coef[-(1:maxP)]
        update.mod.error[v]=x[v]-update.mod.predict[v]
      }
      cl<-makeCluster(min(cores,maxQ))
      registerDoParallel(cl)
      update.dataQ=foreach(q=1:maxQ,.combine=cbind,.export="lag.func")%dopar%{
        lag.func(update.mod.error,k=(q+h-1))
      }
      stopCluster(cl)
      second.modX=as.matrix(cbind(dataP,update.dataQ))[-(1:(max(maxP,maxQ)+h-1)),]
      second.y=x[-(1:(max(maxP,maxQ)+h-1))]
    }else{
      second.modX=first.modX
      second.y=first.y
    }
    
    lambda2.seq=cv.glmnet(y=second.y,x=second.modX,standardize=T,alpha=1,penalty.factor=weights)$lambda
    
    cl<-makeCluster(cores,outfile="debug.txt")
    registerDoParallel(cl)
    SQDEV2=foreach(f=1:nfolds,.packages=c("glmnet"),.combine=rbind)%dopar%{
      in.train=which(Fold.Matrix[,f]==0)
      in.test=which(Fold.Matrix[,f]==1)
      first.mod.est=glmnet(y=second.y[in.train],x=second.modX[in.train,],standardize=T,alpha=1,lambda=lambda2.seq,penalty.factor=weights)
      first.mod.res=(second.y[in.test]-predict(first.mod.est,newx=second.modX[in.test,]))^2
      first.mod.res
    }
    stopCluster(cl)
    
    CVM2=apply(SQDEV2,2,mean)
    lambda2.min=lambda2.seq[which.min(CVM2)]
    lambda2.1se=lambda2.seq[min(which(CVM2<(min(CVM2)+sd(CVM2)/sqrt(length(CVM2)))))]
    
    final.mod.est=glmnet(y=second.y,x=second.modX,standardize=T,alpha=1,lambda=lambda2.1se,penalty.factor=weights)
    final.mod.coef=as.numeric(coef(final.mod.est))[-1]
    nonzero.select=which(final.mod.coef!=0)
    final.mod.int=as.numeric(coef(final.mod.est))[1]
    final.mod.s2=sum((second.y-predict(final.mod.est,newx=second.modX,s=lambda2.1se,method="lambda"))^2)/
      (length(second.y)-sum(final.mod.coef[nonzero.select]!=0)-1)
    
    out=list(final.mod.coef=final.mod.coef,final.mod.int=final.mod.int,final.mod.s2=final.mod.s2,nonzero.select=nonzero.select)
  }
  
  if(CV.method=="RegCV"){
    
    if(is.null(K)){
      cl<-makeCluster(cores,outfile="debug.txt")
      registerDoParallel(cl)
      first.mod.est=cv.glmnet(y=first.y,x=first.modX,standardize=T,alpha=1,parallel=T)
      stopCluster(cl)
    }else{
      cl<-makeCluster(cores,outfile="debug.txt")
      registerDoParallel(cl)
      first.mod.est=cv.glmnet(y=first.y,x=first.modX,standardize=T,alpha=1,nfolds=K,parallel=T)
      stopCluster(cl)
    }
    
    first.mod.coef=coef(first.mod.est,s=first.mod.est$lambda.min)[-1]
    first.mod.mu=coef(first.mod.est,s=first.mod.est$lambda.min)[1]
    weights=abs(first.mod.coef+1/length(first.y))^(-2)
    
    if(updateMA){
      update.mod.predict=rep(NA,length(x))
      update.mod.error=rep(0,length(x))
      for(v in (h+max(maxP,maxQ)):Nt){
        update.mod.predict[v]=first.mod.mu+x[(v-h):(v-maxP-h+1)]%*%first.mod.coef[1:maxP]+
          update.mod.error[(v-h):(v-maxQ-h+1)]%*%first.mod.coef[-(1:maxP)]
        update.mod.error[v]=x[v]-update.mod.predict[v]
      }
      cl<-makeCluster(min(cores,maxQ))
      registerDoParallel(cl)
      update.dataQ=foreach(q=1:maxQ,.combine=cbind,.export="lag.func")%dopar%{
        lag.func(update.mod.error,k=(q+h-1))
      }
      stopCluster(cl)
      second.modX=as.matrix(cbind(dataP,update.dataQ))[-(1:(max(maxP,maxQ)+h-1)),]
      second.y=x[-(1:(max(maxP,maxQ)+h-1))]
    }else{
      second.modX=first.modX
      second.y=first.y
    }
    
    if(is.null(K)){
      cl<-makeCluster(cores,outfile="debug.txt")
      registerDoParallel(cl)
      second.mod.est=cv.glmnet(y=second.y,x=second.modX,standardize=T,alpha=1,penalty.factor=weights)
      stopCluster(cl)
    }else{
      cl<-makeCluster(cores,outfile="debug.txt")
      registerDoParallel(cl)
      second.mod.est=cv.glmnet(y=second.y,x=second.modX,standardize=T,alpha=1,nfolds=K,penalty.factor=weights)
      stopCluster(cl)
      
    }
    
    final.mod.coef=coef(second.mod.est,s=second.mod.est$lambda.1se)[-1]
    nonzero.select=which(final.mod.coef!=0)
    final.mod.int=coef(second.mod.est,s=second.mod.est$lambda.1se)[1]
    final.mod.s2=sum((second.y-predict(second.mod.est,newx=second.modX,s=second.mod.est$lambda.1se,method="lambda"))^2)/
      (length(second.y)-sum(final.mod.coef[nonzero.select]!=0)-1)
    
    out=list(final.mod.coef=final.mod.coef,final.mod.int=final.mod.int,final.mod.s2=final.mod.s2,nonzero.select=nonzero.select)
  }
  
  if(CV.method=="OOSINDEP"){
    Fold.Vector=OOS.IndepCV.func(x=first.y,test.per=test.per)
    
    in.train=which(Fold.Vector==0)
    in.test=which(Fold.Vector==1)
    
    first.mod.est=glmnet(y=first.y[in.train],x=first.modX[in.train,],standardize=T,alpha=1)
    first.mod.res=(first.y[in.test]-predict(first.mod.est,newx=first.modX[in.test,]))^2
    
    OOS.MSE1=apply(first.mod.res,2,mean)
    lambda1.min=first.mod.est$lambda[which.min(OOS.MSE1)]
    lambda1.1se=first.mod.est$lambda[min(which(OOS.MSE1<(min(OOS.MSE1)+sd(OOS.MSE1)/sqrt(length(OOS.MSE1)))))]
    
    first.mod.coef=as.numeric(coef(first.mod.est,s=lambda1.min))[-1]
    first.mod.mu=as.numeric(coef(first.mod.est,s=lambda1.min))[1]
    weights=abs(first.mod.coef+1/length(first.y))^(-eta)
    
    if(updateMA){
      update.mod.predict=rep(NA,length(x))
      update.mod.error=rep(0,length(x))
      for(v in (h+max(maxP,maxQ)):Nt){
        update.mod.predict[v]=first.mod.mu+x[(v-h):(v-maxP-h+1)]%*%first.mod.coef[1:maxP]+
          update.mod.error[(v-h):(v-maxQ-h+1)]%*%first.mod.coef[-(1:maxP)]
        update.mod.error[v]=x[v]-update.mod.predict[v]
      }
      cl<-makeCluster(min(cores,maxQ))
      registerDoParallel(cl)
      update.dataQ=foreach(q=1:maxQ,.combine=cbind,.export="lag.func")%dopar%{
        lag.func(update.mod.error,k=(q+h-1))
      }
      stopCluster(cl)
      second.modX=as.matrix(cbind(dataP,update.dataQ))[-(1:(max(maxP,maxQ)+h-1)),]
      second.y=x[-(1:(max(maxP,maxQ)+h-1))]
    }else{
      second.modX=first.modX
      second.y=first.y
    }
    
    second.mod.est=glmnet(y=second.y[in.train],x=second.modX[in.train,],standardize=T,alpha=1,penalty.factor=weights)
    second.mod.res=(second.y[in.test]-predict(second.mod.est,newx=second.modX[in.test,]))^2
    
    OOS.MSE2=apply(second.mod.res,2,mean)
    lambda2.min=second.mod.est$lambda[which.min(OOS.MSE2)]
    lambda2.1se=second.mod.est$lambda[min(which(OOS.MSE2<(min(OOS.MSE2)+sd(OOS.MSE2)/sqrt(length(OOS.MSE2)))))]
    
    second.mod.coef=as.numeric(coef(second.mod.est,s=lambda2.1se))[-1]
    second.mod.mu=as.numeric(coef(second.mod.est,s=lambda2.1se))[1]

    final.mod.coef=second.mod.coef
    nonzero.select=which(final.mod.coef!=0)
    final.mod.int=second.mod.mu
    final.mod.s2=sum((second.y-predict(second.mod.est,newx=second.modX,s=lambda2.1se,method="lambda"))^2)/
      (length(second.y)-sum(final.mod.coef[nonzero.select]!=0)-1)
    
    out=list(final.mod.coef=final.mod.coef,final.mod.int=final.mod.int,final.mod.s2=final.mod.s2,nonzero.select=nonzero.select)
  }
  
  if(CV.method=="OOSDEP"){
    Fold.Vector=OOS.DepCV.func(x=first.y,test.per=test.per)
    
    in.train=which(Fold.Vector==0)
    in.test=which(Fold.Vector==1)
    
    first.mod.est=glmnet(y=first.y[in.train],x=first.modX[in.train,],standardize=T,alpha=1)
    first.mod.res=(first.y[in.test]-predict(first.mod.est,newx=first.modX[in.test,]))^2
    
    OOS.MSE1=apply(first.mod.res,2,mean)
    lambda1.min=first.mod.est$lambda[which.min(OOS.MSE1)]
    lambda1.1se=first.mod.est$lambda[min(which(OOS.MSE1<(min(OOS.MSE1)+sd(OOS.MSE1)/sqrt(length(OOS.MSE1)))))]
    
    first.mod.coef=as.numeric(coef(first.mod.est,s=lambda1.min))[-1]
    first.mod.mu=as.numeric(coef(first.mod.est,s=lambda1.min))[1]
    weights=abs(first.mod.coef+1/length(first.y))^(-eta)
    
    if(updateMA){
      update.mod.predict=rep(NA,length(x))
      update.mod.error=rep(0,length(x))
      for(v in (h+max(maxP,maxQ)):Nt){
        update.mod.predict[v]=first.mod.mu+x[(v-h):(v-maxP-h+1)]%*%first.mod.coef[1:maxP]+
          update.mod.error[(v-h):(v-maxQ-h+1)]%*%first.mod.coef[-(1:maxP)]
        update.mod.error[v]=x[v]-update.mod.predict[v]
      }
      cl<-makeCluster(min(cores,maxQ))
      registerDoParallel(cl)
      update.dataQ=foreach(q=1:maxQ,.combine=cbind,.export="lag.func")%dopar%{
        lag.func(update.mod.error,k=(q+h-1))
      }
      stopCluster(cl)
      second.modX=as.matrix(cbind(dataP,update.dataQ))[-(1:(max(maxP,maxQ)+h-1)),]
      second.y=x[-(1:(max(maxP,maxQ)+h-1))]
    }else{
      second.modX=first.modX
      second.y=first.y
    }
    
    second.mod.est=glmnet(y=second.y[in.train],x=second.modX[in.train,],standardize=T,alpha=1,penalty.factor=weights)
    second.mod.res=(second.y[in.test]-predict(second.mod.est,newx=second.modX[in.test,]))^2
    
    OOS.MSE2=apply(second.mod.res,2,mean)
    lambda2.min=second.mod.est$lambda[which.min(OOS.MSE2)]
    lambda2.1se=second.mod.est$lambda[min(which(OOS.MSE2<(min(OOS.MSE2)+sd(OOS.MSE2)/sqrt(length(OOS.MSE2)))))]
    
    second.mod.coef=as.numeric(coef(second.mod.est,s=lambda2.1se))[-1]
    second.mod.mu=as.numeric(coef(second.mod.est,s=lambda2.1se))[1]
    
    final.mod.coef=second.mod.coef
    nonzero.select=which(final.mod.coef!=0)
    final.mod.int=second.mod.mu
    final.mod.s2=sum((second.y-predict(second.mod.est,newx=second.modX,s=lambda2.1se,method="lambda"))^2)/
      (length(second.y)-sum(final.mod.coef[nonzero.select]!=0)-1)
    
    out=list(final.mod.coef=final.mod.coef,final.mod.int=final.mod.int,final.mod.s2=final.mod.s2,nonzero.select=nonzero.select)
  }
  
  return(out)
}

#Performs Adaptive Elastic Net
cv.adaptenet.arma.func<-function(x,h=1,long.ar.select=F,maxP,maxQ,max.pq=NULL,K=NULL,updateMA=F,alpha=seq(0,1,0.01),
    test.per=0.2,BIC1=F,BIC2=F,eta=2,CV.method=c("AIC/BIC","LOBOCV","Bergmeir","RegCV","OOSINDEP","OOSDEP"),CV2=c("min","1se")){
  
  #########################
  #Require Certain Packages
  #########################
  require(glmnet) #Performs Ridge, Lasso, Elastic Net
  require(lars)   #Use to Get Lambda Sequence
  require(forecast)
  
  CV.method=match.arg(CV.method)
  CV2=match.arg(CV2)
  
  Nt=length(x)
  max.ar.order=ceiling(10*log10(Nt))
  
  init.mod.est=ar(x,aic=long.ar.select,order.max=max.ar.order,demean=T)
  init.mod.error=residuals(init.mod.est)
  init.mod.order=length(which(is.na(init.mod.error)))
  
  dataP=foreach(p=1:maxP,.combine=cbind)%do%{
    lag.func(x,k=(p+h-1))
  }
  dataQ=foreach(q=1:maxQ,.combine=cbind)%do%{
    lag.func(init.mod.error,k=(q+h-1))
  }
  
  first.modX=as.matrix(cbind(dataP,dataQ))[-(1:(init.mod.order+max(maxP,maxQ)+h-1)),]
  first.y=x[-(1:(init.mod.order+max(maxP,maxQ)+h-1))]
  
  n.alpha=length(alpha)
  
  if(CV.method=="AIC/BIC"){
    
    first.mod.est=glmnet(y=first.y,x=first.modX,standardize=T,alpha=1)
    first.mod.RSS=colSums((first.y-predict(first.mod.est,newx=first.modX))^2)
    if(BIC1){
      first.bic.out=log(length(first.y))*first.mod.est$df+length(first.y)*log(as.vector(first.mod.RSS)/length(first.y))
      first.mod.lambda=first.mod.est$lambda[which.min(first.bic.out)]
      first.cv.out=c(1,first.mod.lambda,min(first.bic.out))
    }else{
      first.aic.out=2*first.mod.est$df+length(first.y)*log(as.vector(first.mod.RSS)/length(first.y))
      first.mod.lambda=first.mod.est$lambda[which.min(first.aic.out)]
      first.cv.out=c(1,first.mod.lambda,min(first.aic.out))
    }
    
    first.mod.alpha=1
    first.mod.lambda=first.cv.out[2]
    first.mod.est=glmnet(y=first.y,x=first.modX,standardize=T,alpha=first.mod.alpha,lambda=first.mod.lambda)
    first.mod.coef=as.numeric(coef(first.mod.est))[-1]
    first.mod.mu=as.numeric(coef(first.mod.est))[1]
    
    weights=abs(first.mod.coef+1/length(first.y))^(-eta)
    
    if(updateMA){
      update.mod.predict=rep(NA,length(x))
      update.mod.error=rep(0,length(x))
      for(v in (h+max(maxP,maxQ)):Nt){
        update.mod.predict[v]=first.mod.mu+x[(v-h):(v-maxP-h+1)]%*%first.mod.coef[1:maxP]+
          update.mod.error[(v-h):(v-maxQ-h+1)]%*%first.mod.coef[-(1:maxP)]
        update.mod.error[v]=x[v]-update.mod.predict[v]
      }
      update.dataQ=foreach(q=1:maxQ,.combine=cbind)%do%{
        lag.func(update.mod.error,k=(q+h-1))
      }
      second.modX=as.matrix(cbind(dataP,update.dataQ))[-(1:(max(maxP,maxQ)+h-1)),]
      second.y=x[-(1:(max(maxP,maxQ)+h-1))]
    }else{
      second.modX=first.modX
      second.y=first.y
    }
    
    second.cv.out=foreach(a=1:n.alpha,.combine=rbind)%do%{
      second.mod.est=glmnet(y=second.y,x=second.modX,standardize=T,alpha=alpha[a],penalty.factor=weights)
      second.mod.RSS=colSums((second.y-predict(second.mod.est,newx=second.modX))^2)
      if(BIC2){
        second.bic.out=log(length(second.y))*second.mod.est$df+length(second.y)*log(as.vector(second.mod.RSS)/length(second.y))
        second.mod.lambda=second.mod.est$lambda[which.min(second.bic.out)]
        result=c(alpha[a],second.mod.lambda,min(second.bic.out))
      }else{
        second.aic.out=2*second.mod.est$df+length(second.y)*log(as.vector(second.mod.RSS)/length(second.y))
        second.mod.lambda=second.mod.est$lambda[which.min(second.aic.out)]
        result=c(alpha[a],second.mod.lambda,min(second.aic.out))
      }
      result
    }  
    
    second.mod.alpha=alpha[which.min(second.cv.out[,3])]
    second.mod.lambda=second.cv.out[which.min(second.cv.out[,3]),2]
    second.mod.est=glmnet(y=second.y,x=second.modX,standardize=T,alpha=second.mod.alpha,lambda=second.mod.lambda,penalty.factor=weights)
    second.mod.coef=as.numeric(coef(second.mod.est))[-1]
    second.mod.mu=as.numeric(coef(second.mod.est))[1]
    
    final.mod.coef=second.mod.coef
    nonzero.select=which(final.mod.coef!=0)
    final.mod.int=second.mod.mu
    final.mod.s2=sum((second.y-predict(second.mod.est,newx=second.modX))^2)/(length(second.y)-sum(final.mod.coef[nonzero.select]!=0)-1)
    
    out=list(final.mod.coef=final.mod.coef,final.mod.int=final.mod.int,final.mod.s2=final.mod.s2,nonzero.select=nonzero.select)
  }
  
  if(CV.method=="LOBOCV"){
    Fold.Matrix=NonDepCV1.func(first.y,max.pq=max.pq,K=K)
    nfolds=dim(Fold.Matrix)[2]
    
    first.cv.out=NULL
    lambda1.seq=cv.glmnet(parallel=F,y=first.y,x=first.modX,standardize=T,alpha=1)$lambda
    SQDEV1=foreach(f=1:nfolds,.combine=rbind)%do%{
        in.train=which(Fold.Matrix[,f]==0)
        in.test=which(Fold.Matrix[,f]==1)
        first.mod.est=glmnet(y=first.y[in.train],x=first.modX[in.train,],standardize=T,alpha=1,lambda=lambda1.seq)
        first.mod.res=(first.y[in.test]-predict(first.mod.est,newx=first.modX[in.test,]))^2
        first.mod.res
    }
    CVM1=apply(SQDEV1,2,mean)
    lambda1.min=lambda1.seq[which.min(CVM1)]
    lambda1.1se=lambda1.seq[min(which(CVM1<(min(CVM1)+sd(CVM1)/sqrt(length(CVM1)))))]
    first.cv.out=rbind(first.cv.out,c(1,lambda1.min,min(CVM1)))
    
    
    first.mod.alpha=1
    first.mod.lambda=first.cv.out[which.min(first.cv.out[,3]),2]
    first.mod.est=glmnet(y=first.y,x=first.modX,standardize=T,alpha=first.mod.alpha,lambda=first.mod.lambda)
    first.mod.coef=as.numeric(coef(first.mod.est))[-1]
    first.mod.mu=as.numeric(coef(first.mod.est))[1]
    
    weights=abs(first.mod.coef+1/length(first.y))^(-eta)
  
    if(updateMA){
      update.mod.predict=rep(NA,length(x))
      update.mod.error=rep(0,length(x))
      for(v in (h+max(maxP,maxQ)):Nt){
        update.mod.predict[v]=first.mod.mu+x[(v-h):(v-maxP-h+1)]%*%first.mod.coef[1:maxP]+
          update.mod.error[(v-h):(v-maxQ-h+1)]%*%first.mod.coef[-(1:maxP)]
        update.mod.error[v]=x[v]-update.mod.predict[v]
      }
      update.dataQ=foreach(q=1:maxQ,.combine=cbind)%do%{
        lag.func(update.mod.error,k=(q+h-1))
      }
      second.modX=as.matrix(cbind(dataP,update.dataQ))[-(1:(max(maxP,maxQ)+h-1)),]
      second.y=x[-(1:(max(maxP,maxQ)+h-1))]
    }else{
      second.modX=first.modX
      second.y=first.y
    }
    
    second.cv.out=NULL
    for(a in 1:n.alpha){
      lambda2.seq=cv.glmnet(parallel=F,y=second.y,x=second.modX,standardize=T,alpha=alpha[a],penalty.factor=weights)$lambda
      SQDEV2=foreach(f=1:nfolds,.combine=rbind)%do%{
        in.train=which(Fold.Matrix[,f]==0)
        in.test=which(Fold.Matrix[,f]==1)
        second.mod.est=glmnet(y=second.y[in.train],x=second.modX[in.train,],standardize=T,
                              alpha=alpha[a],penalty.factor=weights,lambda=lambda2.seq)
        second.mod.res=(second.y[in.test]-predict(second.mod.est,newx=second.modX[in.test,]))^2
        second.mod.res
      }
      CVM2=apply(SQDEV2,2,mean)
      lambda2.min=lambda2.seq[which.min(CVM2)]
      lambda2.1se=lambda2.seq[min(which(CVM2<(min(CVM2)+sd(CVM2)/sqrt(length(CVM2)))))]
      if(CV2=="min") second.cv.out=rbind(second.cv.out,c(alpha[a],lambda2.min,min(CVM2)))
      if(CV2=="1se") second.cv.out=rbind(second.cv.out,c(alpha[a],lambda2.1se,CVM2[min(which(CVM2<(min(CVM2)+sd(CVM2)/sqrt(length(CVM2)))))]))
    }
    
    second.mod.alpha=alpha[which.min(second.cv.out[,3])]
    second.mod.lambda=second.cv.out[which.min(second.cv.out[,3]),2]
    second.mod.est=glmnet(y=second.y,x=second.modX,standardize=T,alpha=second.mod.alpha,lambda=second.mod.lambda,penalty.factor=weights)
    second.mod.coef=as.numeric(coef(second.mod.est))[-1]
    second.mod.mu=as.numeric(coef(second.mod.est))[1]
    
    final.mod.coef=second.mod.coef
    nonzero.select=which(final.mod.coef!=0)
    final.mod.int=second.mod.mu
    final.mod.s2=sum((second.y-predict(second.mod.est,newx=second.modX))^2)/(length(second.y)-sum(final.mod.coef[nonzero.select]!=0)-1)
    
    out=list(final.mod.coef=final.mod.coef,final.mod.int=final.mod.int,final.mod.s2=final.mod.s2,nonzero.select=nonzero.select)
  }
  
  if(CV.method=="Bergmeir"){
    Fold.Matrix=NonDepCV2.func(first.y,max.pq=max.pq,K=K)
    nfolds=dim(Fold.Matrix)[2]
    
    first.cv.out=NULL
    lambda1.seq=cv.glmnet(parallel=F,y=first.y,x=first.modX,standardize=T,alpha=1)$lambda
    SQDEV1=foreach(f=1:nfolds,.combine=rbind)%do%{
      in.train=which(Fold.Matrix[,f]==0)
      in.test=which(Fold.Matrix[,f]==1)
      first.mod.est=glmnet(y=first.y[in.train],x=first.modX[in.train,],standardize=T,alpha=1,lambda=lambda1.seq)
      first.mod.res=(first.y[in.test]-predict(first.mod.est,newx=first.modX[in.test,]))^2
      first.mod.res
    }
    CVM1=apply(SQDEV1,2,mean)
    lambda1.min=lambda1.seq[which.min(CVM1)]
    lambda1.1se=lambda1.seq[min(which(CVM1<(min(CVM1)+sd(CVM1)/sqrt(length(CVM1)))))]
    first.cv.out=rbind(first.cv.out,c(1,lambda1.min,min(CVM1)))
    
    
    first.mod.alpha=1
    first.mod.lambda=first.cv.out[which.min(first.cv.out[,3]),2]
    first.mod.est=glmnet(y=first.y,x=first.modX,standardize=T,alpha=first.mod.alpha,lambda=first.mod.lambda)
    first.mod.coef=as.numeric(coef(first.mod.est))[-1]
    first.mod.mu=as.numeric(coef(first.mod.est))[1]
    
    weights=abs(first.mod.coef+1/length(first.y))^(-eta)
    
    if(updateMA){
      update.mod.predict=rep(NA,length(x))
      update.mod.error=rep(0,length(x))
      for(v in (h+max(maxP,maxQ)):Nt){
        update.mod.predict[v]=first.mod.mu+x[(v-h):(v-maxP-h+1)]%*%first.mod.coef[1:maxP]+
          update.mod.error[(v-h):(v-maxQ-h+1)]%*%first.mod.coef[-(1:maxP)]
        update.mod.error[v]=x[v]-update.mod.predict[v]
      }
      update.dataQ=foreach(q=1:maxQ,.combine=cbind)%do%{
        lag.func(update.mod.error,k=(q+h-1))
      }
      second.modX=as.matrix(cbind(dataP,update.dataQ))[-(1:(max(maxP,maxQ)+h-1)),]
      second.y=x[-(1:(max(maxP,maxQ)+h-1))]
    }else{
      second.modX=first.modX
      second.y=first.y
    }
    
    second.cv.out=NULL
    for(a in 1:n.alpha){
      lambda2.seq=cv.glmnet(parallel=F,y=second.y,x=second.modX,standardize=T,alpha=alpha[a],penalty.factor=weights)$lambda
      SQDEV2=foreach(f=1:nfolds,.combine=rbind)%do%{
        in.train=which(Fold.Matrix[,f]==0)
        in.test=which(Fold.Matrix[,f]==1)
        second.mod.est=glmnet(y=second.y[in.train],x=second.modX[in.train,],standardize=T,
                              alpha=alpha[a],penalty.factor=weights,lambda=lambda2.seq)
        second.mod.res=(second.y[in.test]-predict(second.mod.est,newx=second.modX[in.test,]))^2
        second.mod.res
      }
      CVM2=apply(SQDEV2,2,mean)
      lambda2.min=lambda2.seq[which.min(CVM2)]
      lambda2.1se=lambda2.seq[min(which(CVM2<(min(CVM2)+sd(CVM2)/sqrt(length(CVM2)))))]
      if(CV2=="min") second.cv.out=rbind(second.cv.out,c(alpha[a],lambda2.min,min(CVM2)))
      if(CV2=="1se") second.cv.out=rbind(second.cv.out,c(alpha[a],lambda2.1se,CVM2[min(which(CVM2<(min(CVM2)+sd(CVM2)/sqrt(length(CVM2)))))]))
    }
    
    second.mod.alpha=alpha[which.min(second.cv.out[,3])]
    second.mod.lambda=second.cv.out[which.min(second.cv.out[,3]),2]
    second.mod.est=glmnet(y=second.y,x=second.modX,standardize=T,alpha=second.mod.alpha,lambda=second.mod.lambda,penalty.factor=weights)
    second.mod.coef=as.numeric(coef(second.mod.est))[-1]
    second.mod.mu=as.numeric(coef(second.mod.est))[1]
    
    final.mod.coef=second.mod.coef
    nonzero.select=which(final.mod.coef!=0)
    final.mod.int=second.mod.mu
    final.mod.s2=sum((second.y-predict(second.mod.est,newx=second.modX))^2)/(length(second.y)-sum(final.mod.coef[nonzero.select]!=0)-1)
    
    out=list(final.mod.coef=final.mod.coef,final.mod.int=final.mod.int,final.mod.s2=final.mod.s2,nonzero.select=nonzero.select)
  }
  
  if(CV.method=="RegCV"){
   
    if(is.null(K)){
        first.mod.est=cv.glmnet(parallel=F,y=first.y,x=first.modX,standardize=T,alpha=1)
        first.cv.out=c(1,first.mod.est$lambda.min,first.mod.est$cvm[which(first.mod.est$lambda==first.mod.est$lambda.min)])
    }else{
        first.mod.est=cv.glmnet(parallel=F,y=first.y,x=first.modX,standardize=T,alpha=1,nfolds=K)
        first.cv.out=c(1,first.mod.est$lambda.min,first.mod.est$cvm[which(first.mod.est$lambda==first.mod.est$lambda.min)])
    }
    
    first.mod.alpha=1
    first.mod.lambda=first.cv.out[2]
    first.mod.est=glmnet(y=first.y,x=first.modX,standardize=T,alpha=first.mod.alpha,lambda=first.mod.lambda)
    first.mod.coef=as.numeric(coef(first.mod.est))[-1]
    first.mod.mu=as.numeric(coef(first.mod.est))[1]
    
    weights=abs(first.mod.coef+1/length(first.y))^(-eta)
    
    
    if(updateMA){
      update.mod.predict=rep(NA,length(x))
      update.mod.error=rep(0,length(x))
      for(v in (h+max(maxP,maxQ)):Nt){
        update.mod.predict[v]=first.mod.mu+x[(v-h):(v-maxP-h+1)]%*%first.mod.coef[1:maxP]+
          update.mod.error[(v-h):(v-maxQ-h+1)]%*%first.mod.coef[-(1:maxP)]
        update.mod.error[v]=x[v]-update.mod.predict[v]
      }
      update.dataQ=foreach(q=1:maxQ,.combine=cbind)%do%{
        lag.func(update.mod.error,k=(q+h-1))
      }
      second.modX=as.matrix(cbind(dataP,update.dataQ))[-(1:(max(maxP,maxQ)+h-1)),]
      second.y=x[-(1:(max(maxP,maxQ)+h-1))]
    }else{
      second.modX=first.modX
      second.y=first.y
    }
    
    n.alpha=length(alpha)
    
    if(is.null(K)){
      second.cv.out=NULL
      for(a in 1:n.alpha){  
        second.mod.est=cv.glmnet(parallel=F,y=second.y,x=second.modX,standardize=T,alpha=alpha[a],penalty.factor=weights)
        if(CV2=="min") second.cv.out=rbind(second.cv.out,c(alpha[a],second.mod.est$lambda.min,second.mod.est$cvm[which(second.mod.est$lambda==second.mod.est$lambda.min)]))
        if(CV2=="1se") second.cv.out=rbind(second.cv.out,c(alpha[a],second.mod.est$lambda.1se,second.mod.est$cvm[which(second.mod.est$lambda==second.mod.est$lambda.1se)]))
      }
    }else{
      second.cv.out=NULL
      for(a in 1:n.alpha){  
        second.mod.est=cv.glmnet(parallel=F,y=second.y,x=second.modX,standardize=T,alpha=alpha[a],nfolds=K,penalty.factor=weights)
        if(CV2=="min") second.cv.out=rbind(second.cv.out,c(alpha[a],second.mod.est$lambda.min,second.mod.est$cvm[which(second.mod.est$lambda==second.mod.est$lambda.min)]))
        if(CV2=="1se") second.cv.out=rbind(second.cv.out,c(alpha[a],second.mod.est$lambda.1se,second.mod.est$cvm[which(second.mod.est$lambda==second.mod.est$lambda.1se)]))
      }
    }
    
    second.mod.alpha=alpha[which.min(second.cv.out[,3])]
    second.mod.lambda=second.cv.out[which.min(second.cv.out[,3]),2]
    second.mod.est=glmnet(y=second.y,x=second.modX,standardize=T,alpha=second.mod.alpha,lambda=second.mod.lambda,penalty.factor=weights)
    second.mod.coef=as.numeric(coef(second.mod.est))[-1]
    second.mod.mu=as.numeric(coef(second.mod.est))[1]
    
    final.mod.coef=second.mod.coef
    nonzero.select=which(final.mod.coef!=0)
    final.mod.int=second.mod.mu
    final.mod.s2=sum((second.y-predict(second.mod.est,newx=second.modX))^2)/(length(second.y)-sum(final.mod.coef[nonzero.select]!=0)-1)
    
    out=list(final.mod.coef=final.mod.coef,final.mod.int=final.mod.int,final.mod.s2=final.mod.s2,nonzero.select=nonzero.select)
  }
  
  if(CV.method=="OOSINDEP"){
    Fold.Vector=OOS.IndepCV.func(x=first.y,test.per=test.per)
    
    in.train=which(Fold.Vector==0)
    in.test=which(Fold.Vector==1)
   
    first.mod.est=glmnet(y=first.y[in.train],x=first.modX[in.train,],standardize=T,alpha=1)
    first.mod.res=(first.y[in.test]-predict(first.mod.est,newx=first.modX[in.test,]))^2
    OOS.MSE1=apply(first.mod.res,2,mean)
    lambda1.min=first.mod.est$lambda[which.min(OOS.MSE1)]
    lambda1.1se=first.mod.est$lambda[min(which(OOS.MSE1<(min(OOS.MSE1)+sd(OOS.MSE1)/sqrt(length(OOS.MSE1)))))]
    first.cv.out=c(1,lambda1.min,OOS.MSE1[which.min(OOS.MSE1)])
    
    first.mod.alpha=1
    first.mod.lambda=first.cv.out[2]
    first.mod.est=glmnet(y=first.y,x=first.modX,standardize=T,alpha=first.mod.alpha,lambda=first.mod.lambda)
    first.mod.coef=as.numeric(coef(first.mod.est))[-1]
    first.mod.mu=as.numeric(coef(first.mod.est))[1]
    
    weights=abs(first.mod.coef+1/length(first.y))^(-eta)
    
    if(updateMA){
      update.mod.predict=rep(NA,length(x))
      update.mod.error=rep(0,length(x))
      for(v in (h+max(maxP,maxQ)):Nt){
        update.mod.predict[v]=first.mod.mu+x[(v-h):(v-maxP-h+1)]%*%first.mod.coef[1:maxP]+
          update.mod.error[(v-h):(v-maxQ-h+1)]%*%first.mod.coef[-(1:maxP)]
        update.mod.error[v]=x[v]-update.mod.predict[v]
      }
      update.dataQ=foreach(q=1:maxQ,.combine=cbind)%do%{
        lag.func(update.mod.error,k=(q+h-1))
      }
      second.modX=as.matrix(cbind(dataP,update.dataQ))[-(1:(max(maxP,maxQ)+h-1)),]
      second.y=x[-(1:(max(maxP,maxQ)+h-1))]
    }else{
      second.modX=first.modX
      second.y=first.y
    }
    
    second.cv.out=foreach(a=1:n.alpha,.combine=rbind)%do%{
      second.mod.est=glmnet(y=second.y[in.train],x=second.modX[in.train,],standardize=T,alpha=alpha[a])
      second.mod.res=(second.y[in.test]-predict(second.mod.est,newx=second.modX[in.test,]))^2
      OOS.MSE2=apply(second.mod.res,2,mean)
      lambda2.min=second.mod.est$lambda[which.min(OOS.MSE2)]
      lambda2.1se=second.mod.est$lambda[min(which(OOS.MSE2<(min(OOS.MSE2)+sd(OOS.MSE2)/sqrt(length(OOS.MSE2)))))]
      if(CV2=="min") out=c(alpha[a],lambda2.min,min(OOS.MSE2))
      if(CV2=="1se") out=c(alpha[a],lambda2.1se,OOS.MSE2[min(which(OOS.MSE2<(min(OOS.MSE2)+sd(OOS.MSE2)/sqrt(length(OOS.MSE2)))))])
      out
    }
    
    second.mod.alpha=alpha[which.min(second.cv.out[,3])]
    second.mod.lambda=second.cv.out[which.min(second.cv.out[,3]),2]
    second.mod.est=glmnet(y=second.y,x=second.modX,standardize=T,alpha=second.mod.alpha,lambda=second.mod.lambda,penalty.factor=weights)
    second.mod.coef=as.numeric(coef(second.mod.est))[-1]
    second.mod.mu=as.numeric(coef(second.mod.est))[1]
    
    final.mod.coef=second.mod.coef
    nonzero.select=which(final.mod.coef!=0)
    final.mod.int=second.mod.mu
    final.mod.s2=sum((second.y-predict(second.mod.est,newx=second.modX))^2)/(length(second.y)-sum(final.mod.coef[nonzero.select]!=0)-1)
    
    out=list(final.mod.coef=final.mod.coef,final.mod.int=final.mod.int,final.mod.s2=final.mod.s2,nonzero.select=nonzero.select)
  }
  
  if(CV.method=="OOSDEP"){
    Fold.Vector=OOS.DepCV.func(x=first.y,test.per=test.per)
    
    in.train=which(Fold.Vector==0)
    in.test=which(Fold.Vector==1)
    
    first.mod.est=glmnet(y=first.y[in.train],x=first.modX[in.train,],standardize=T,alpha=1)
    first.mod.res=(first.y[in.test]-predict(first.mod.est,newx=first.modX[in.test,]))^2
    OOS.MSE1=apply(first.mod.res,2,mean)
    lambda1.min=first.mod.est$lambda[which.min(OOS.MSE1)]
    lambda1.1se=first.mod.est$lambda[min(which(OOS.MSE1<(min(OOS.MSE1)+sd(OOS.MSE1)/sqrt(length(OOS.MSE1)))))]
    first.cv.out=c(1,lambda1.min,OOS.MSE1[which.min(OOS.MSE1)])
    
    first.mod.alpha=1
    first.mod.lambda=first.cv.out[2]
    first.mod.est=glmnet(y=first.y,x=first.modX,standardize=T,alpha=first.mod.alpha,lambda=first.mod.lambda)
    first.mod.coef=as.numeric(coef(first.mod.est))[-1]
    first.mod.mu=as.numeric(coef(first.mod.est))[1]
    
    weights=abs(first.mod.coef+1/length(first.y))^(-eta)
    
    if(updateMA){
      update.mod.predict=rep(NA,length(x))
      update.mod.error=rep(0,length(x))
      for(v in (h+max(maxP,maxQ)):Nt){
        update.mod.predict[v]=first.mod.mu+x[(v-h):(v-maxP-h+1)]%*%first.mod.coef[1:maxP]+
          update.mod.error[(v-h):(v-maxQ-h+1)]%*%first.mod.coef[-(1:maxP)]
        update.mod.error[v]=x[v]-update.mod.predict[v]
      }
      update.dataQ=foreach(q=1:maxQ,.combine=cbind)%do%{
        lag.func(update.mod.error,k=(q+h-1))
      }
      second.modX=as.matrix(cbind(dataP,update.dataQ))[-(1:(max(maxP,maxQ)+h-1)),]
      second.y=x[-(1:(max(maxP,maxQ)+h-1))]
    }else{
      second.modX=first.modX
      second.y=first.y
    }
    
    second.cv.out=foreach(a=1:n.alpha,.combine=rbind)%do%{
      second.mod.est=glmnet(y=second.y[in.train],x=second.modX[in.train,],standardize=T,alpha=alpha[a])
      second.mod.res=(second.y[in.test]-predict(second.mod.est,newx=second.modX[in.test,]))^2
      OOS.MSE2=apply(second.mod.res,2,mean)
      lambda2.min=second.mod.est$lambda[which.min(OOS.MSE2)]
      lambda2.1se=second.mod.est$lambda[min(which(OOS.MSE2<(min(OOS.MSE2)+sd(OOS.MSE2)/sqrt(length(OOS.MSE2)))))]
      if(CV2=="min") out=c(alpha[a],lambda2.min,min(OOS.MSE2))
      if(CV2=="1se") out=c(alpha[a],lambda2.1se,OOS.MSE2[min(which(OOS.MSE2<(min(OOS.MSE2)+sd(OOS.MSE2)/sqrt(length(OOS.MSE2)))))])
      out
    }
    
    second.mod.alpha=alpha[which.min(second.cv.out[,3])]
    second.mod.lambda=second.cv.out[which.min(second.cv.out[,3]),2]
    second.mod.est=glmnet(y=second.y,x=second.modX,standardize=T,alpha=second.mod.alpha,lambda=second.mod.lambda,penalty.factor=weights)
    second.mod.coef=as.numeric(coef(second.mod.est))[-1]
    second.mod.mu=as.numeric(coef(second.mod.est))[1]
    
    final.mod.coef=second.mod.coef
    nonzero.select=which(final.mod.coef!=0)
    final.mod.int=second.mod.mu
    final.mod.s2=sum((second.y-predict(second.mod.est,newx=second.modX))^2)/(length(second.y)-sum(final.mod.coef[nonzero.select]!=0)-1)
    
    out=list(final.mod.coef=final.mod.coef,final.mod.int=final.mod.int,final.mod.s2=final.mod.s2,nonzero.select=nonzero.select)
  }
  
  return(out)
}
#############################################################################################





#############################################################################################
#Function to Conduct ARMA Selection via Bayesian Projection Posterior Predictive Distribution
#   Implementing Cross Validation for Final Model Selection
#Arguments: x=Time Series to Be Modeled Using ARMA Process
#           h=Horizon Specific Model (Defaults to 1)
#           cores=Number of Cores for Parallel Processing (Defaults to 3)
#           maxP=Maximum Autoregressive Order 
#           maxQ=Maximum Moving Average Order
#           KL.threshold=Single value or vector of chosen thresholds for stopping rule
#                        based on Relative Efficiency Comparing Submodel to Full Model
#                        based on Kullback Leibler Divergence
#                        (Defaults to c(0.9,0.95,0.99))
#Source: Piironen and Vehtari (2015)
#############################################################################################

cv.pms.arma.func<-function(x,h=1,cores=1,maxP,maxQ,KL.stop=0.98,cv.type=c("ind","dep"),
                           test.per=0.2,prior.choice=c("hs","hs+"),updateMA=F){
  require(MCMCpack)
  require(bayesreg)
  require(doParallel)
  require(foreach)
  
  prior.choice=match.arg(prior.choice)
  cv.type=match.arg(cv.type)
  
  if(cv.type=="ind"){
    cv.vector=OOS.IndepCV.func(x,test.per=test.per)
  }else if(cv.type=="dep"){
    cv.vector=OOS.DepCV.func(x,test.per=test.per,max.pq=max(maxP,maxQ))
  }
  
  x.train=x[cv.vector==0]
  x.test=x[cv.vector==1]
  
  Nt=length(x.train)
  max.ar.order=ceiling(10*log10(Nt))
  
  cl<-makeCluster(min(cores,max.ar.order))
  registerDoParallel(cl)
  init.modX=foreach(init.ar=1:max.ar.order,.combine=cbind,.export="lag.func")%dopar%{
    lag.func(x.train,k=(init.ar+h-1))
  }
  stopCluster(cl)
  
  init.data=data.frame(y=x.train,init.modX)
  init.data=init.data[-(1:(max.ar.order+h-1)),]
  
  init.mod.est=MCMCregress(y~.,data=init.data,mcmc=2000,thin=10,burnin=10000)
  muBeta0=mean(init.mod.est[,1])
  muBeta=colMeans(init.mod.est[,-c(1,dim(init.mod.est)[2])])
  init.mod.error=init.data$y-(as.numeric(muBeta0) + as.matrix(init.data[,-1])%*%as.vector(muBeta))
  
  cl<-makeCluster(min(cores,maxP))
  registerDoParallel(cl)
  dataP=foreach(p=1:maxP,.combine=cbind,.export="lag.func")%dopar%{
    lag.func(init.data$y,k=(p+h-1))
  }
  stopCluster(cl)
  
  cl<-makeCluster(min(cores,maxQ))
  registerDoParallel(cl)
  dataQ=foreach(q=1:maxQ,.combine=cbind,.export="lag.func")%dopar%{
    lag.func(init.mod.error,k=(q+h-1))
  }
  stopCluster(cl)
  
  full.data=data.frame(y=init.data$y,dataP=dataP,dataQ=dataQ)
  full.data=full.data[-(1:(max(maxP,maxQ)+h-1)),]
  
  #############################################
  #xc.mean=as.numeric(colMeans(full.data))
  #xc.sd=as.numeric(apply(full.data,2,sd))
  #full.data=as.data.frame(scale(full.data))
  ##############################################
  
  full.mod.est=bayesreg(y~.,data=full.data,prior=prior.choice,nsamples=2000,thin=10,burnin=10000)  
  full.mod.posterior=cbind(as.vector(full.mod.est$beta0),t(as.matrix(full.mod.est$beta)),as.vector(full.mod.est$sigma2))
  
  full.mod.int=c(full.mod.est$muBeta0)
  full.mod.coef=c(full.mod.est$muBeta)
  full.mod.s2=full.mod.est$muSigma2
  
  if(updateMA){
    full.mod.predict=rep(NA,length(x.train))
    full.mod.error=rep(0,length(x.train))
    for(k in (h+max(maxP,maxQ)):Nt){
      full.mod.predict[k]=full.mod.int+x.train[(k-h):(k-maxP-h+1)]%*%full.mod.coef[1:maxP]+full.mod.error[(k-h):(k-maxQ-h+1)]%*%full.mod.coef[-(1:maxP)]
      full.mod.error[k]=x.train[k]-full.mod.predict[k]
    }
    
    cl<-makeCluster(min(cores,maxP))
    registerDoParallel(cl)
    full.dataP=foreach(p=1:maxP,.combine=cbind,.export="lag.func")%dopar%{
      lag.func(x.train,k=(p+h-1))
    }
    stopCluster(cl)
    
    cl<-makeCluster(min(cores,maxQ))
    registerDoParallel(cl)
    update.dataQ=foreach(p=1:maxQ,.combine=cbind,.export="lag.func")%dopar%{
      lag.func(full.mod.error,k=(p+h-1))
    }
    stopCluster(cl)
    
    full.mod.X=cbind(1,full.dataP,update.dataQ)
    full.mod.X=full.mod.X[-(1:(max(maxP,maxQ)+h-1)),]
    nMod=dim(full.mod.X)[2]
  }else{
    full.mod.X=as.matrix(cbind(1,full.data[,-1]))
    nMod=dim(full.mod.X)[2]
  }
  
  KL=rep(NA,nMod)
  chosen=1
  notchosen=setdiff(1:nMod,chosen)
  FIRST=pms.proj.func(fullpost=full.mod.posterior,X=full.mod.X,indproj=chosen)
  KL[1]=FIRST$KL.MEAN
  modnum=2
  check=0
  
  while(check<KL.stop | modnum==nMod){
    nleft<-length(notchosen)
    cl<-makeCluster(cores)
    registerDoParallel(cl)
    val=foreach(j=1:nleft,.combine=c,.export="pms.proj.func")%dopar%{
      ind<-sort(c(chosen,notchosen[j]))
      NEXT<-tryCatch({pms.proj.func(fullpost=full.mod.posterior,X=full.mod.X,indproj=ind)$KL.MEAN},
                     error=function(e){return(NA)})
      NEXT
    }
    stopCluster(cl)
    minval<-which.min(val)
    chosen<-c(chosen,notchosen[minval])
    notchosen<-setdiff(1:nMod,chosen)
    KL[modnum]<-val[minval]
    check=1-KL[modnum]/KL[1]
    modnum=modnum+1
  }
  
  KL=KL[!is.na(KL)]
  chosen=chosen[!is.na(chosen)]
  
  nMod2=length(KL)
  MSE=rep(NA,nMod2)
  
  for(modnum in 1:nMod2){
    out=pms.proj.func(fullpost=full.mod.posterior,X=full.mod.X,indproj=chosen[1:modnum])
    coef=out$COEF.MEAN
    int=coef[1]
    coef.ar=coef[2:(1+maxP)]
    coef.ma=coef[-(1:(maxP+1))]
    
    predictx.test=rep(NA,length(x.test))
    errorx.test=rep(0,length(x.test))
    for(k in (h+max(maxP,maxQ)):length(x.test)){
      predictx.test[k]=int+x.test[(k-h):(k-maxP-h+1)]%*%coef.ar+errorx.test[(k-h):(k-maxQ-h+1)]%*%coef.ma
      errorx.test[k]=x.test[k]-predictx.test[k]
    }
    
    MSE[modnum]=mean((x.test-predictx.test)^2,na.rm=T)
  }
  best.mod=chosen[1:which.min(MSE)]
  out.mod=pms.proj.func(fullpost=full.mod.posterior,X=full.mod.X,indproj=best.mod)
  
  final.mod.int=out.mod$COEF.MEAN[1]
  final.mod.coef=out.mod$COEF.MEAN[-1]
  final.mod.s2=out.mod$VAR.MEAN
  
  return(list(CHOICE=chosen,KL=KL,full.mod.int=full.mod.int,full.mod.coef=full.mod.coef,full.mod.s2=full.mod.s2,
              final.mod=best.mod,final.mod.int=final.mod.int,final.mod.coef=final.mod.coef,final.mod.s2=final.mod.s2))
}    
#############################################################################################



