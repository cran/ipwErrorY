#' Doubly Robust Estimation of ATE with Known Error
#'
#' Doubly robust estimation of average treatment effect with known outcome misclassification probabilities, i.e., known sensitivity and specificity
#'@param data The dataset to be analyzed in the form of R data frame
#'@param indA A column name indicating the treatment variable
#'@param indYerror A column name indicating the misclassified outcome variable
#'@param indXtrt A vector of column names indicating the covariates included in the treatment model
#'@param indXout A vector of column names indicating the covariates included in the outcome model
#'@param sensitivity The specified sensitivity between 0 and 1
#'@param specificity The specified specificity between 0 and 1
#'@param numBoot The specified number of bootstrap replicates for variance estimation
#'@param sharePara if the treated and untreated groups share parameters for covariates in the logistic outcome model (i.e., assuming Y~ T+X), then set \code{sharePara=TRUE}; if not (i.e., modeling Y~ X for the treated and untreated groups separately), then set \code{sharePara=FALSE}. By default,  \code{sharePara=FALSE}
#'@param confidence The confidence level between 0 and 1; the default is 0.95 corresponding to a 95 per cent confidence interval
#'@return A list of the estimate of average treatment effect, bootstrap standard error and confidence interval
#'@import nleqslv
#'@import stats
#'
#'@examples
#'#create a dataset with sensitivity=0.95 and specificity=0.85
#'set.seed(100)
#'X=rnorm(2000)   
#'xx=X^2
#'A=rbinom(2000,1,1/(1+exp(-0.1-X-0.2*xx)))
#'Y=rbinom(2000,1,1/(1+exp(1-A-0.5*X-xx)))
#'y1=which(Y==1)
#'y0=which(Y==0) 
#'Y[y1]=rbinom(length(y1),1,0.95)
#'Y[y0]=rbinom(length(y0),1,0.15)
#'Yast=Y
#'da=data.frame(A=A,X=X,xx=xx,Yast=Yast)
#'head(da)
#'#apply the doubly robust correction method with sensitivity=0.95 and specificity=0.85
#'set.seed(100)
#'KnownErrorDR(da,"A","Yast",c("X","xx"),c("X","xx"),0.95,0.85,50,FALSE,0.95)
#'
#'@export
KnownErrorDR<-function(data,indA,indYerror,indXtrt,indXout,sensitivity,specificity,numBoot,sharePara=FALSE,confidence=0.95){
  p11=sensitivity
  p10=1-specificity 
  n=nrow(data)
  
  once<-function(da){
  n=nrow(da)
  A=da[,indA]
  iAi=which(A==1)
  daTrt=da[iAi,]
  daOut=da[-iAi,]
  Y=da[,indYerror]
  nXtrt=length(indXtrt)+1
  nXout=length(indXout)+1
  trtX0=da[,indXtrt]
  trtX=cbind(rep(1,n),trtX0)
  outX0=da[,indXout]
  outX=as.matrix(cbind(rep(1,n),outX0))
  nT=sum(A)
  Ytrt=daTrt[,indYerror]
  outX0Trt=daTrt[,indXout]
  outXTrt=as.matrix(cbind(rep(1,nT),outX0Trt))
  Yout=daOut[,indYerror]
  outX0Out=daOut[,indXout]
  outXOut=as.matrix(cbind(rep(1,n-nT),outX0Out))
  
  gg=glm(A~.,family="binomial",data=as.data.frame(cbind(A,trtX0)))
  psfit=predict(gg, type = "response")
  
  if (sharePara==FALSE){
    #for treament group
  ggout=glm(Ytrt~.,family="binomial",data=as.data.frame(cbind(Ytrt,outX0Trt)))
  intc=ggout$coef
    
    fc=function(gma){
      gma=as.matrix(gma)
      -sum(log(1/(1+exp(-outXTrt%*%gma))*(p11*Ytrt+(1-p11)*(1-Ytrt))
          +(1-1/(1+exp(-outXTrt%*%gma)))*(p10*Ytrt+(1-p10)*(1-Ytrt))))
    } 
    estgma=optim(intc,fc)$par
  
  m1=1/(1+exp(-outX%*%estgma))
  
  ggout=glm(Yout~.,family="binomial",data=as.data.frame(cbind(Yout,outX0Out)))
  intc=ggout$coef
  
  fc=function(gma){
    
    gma=as.matrix(gma)
    -sum(log(1/(1+exp(-outXOut%*%gma))*(p11*Yout+(1-p11)*(1-Yout))
             +(1-1/(1+exp(-outXOut%*%gma)))*(p10*Yout+(1-p10)*(1-Yout))))
  } 
  estgma=optim(intc,fc)$par
  m0=1/(1+exp(-outX%*%estgma))
  
    mu1=mean(A*Y/(psfit*(p11-p10))-(A-psfit)/(psfit)*m1-
               A/psfit*(p10/(p11-p10)))
    mu0=mean((1-A)*Y/((1-psfit)*(p11-p10))+(A-psfit)/(1-psfit)*m0-
               (1-A)/(1-psfit)*(p10/(p11-p10)))
   mu1-mu0   
  }
  
  else if (sharePara==TRUE){
    ggout=glm(Y~.,family="binomial",data=as.data.frame(cbind(Y,outX0,A)))
    intc=ggout$coef
    
    fc=function(gma){
      Agma=gma[nXout+1]
      gmaX=gma[1:nXout]
      gmaX=as.matrix(gmaX)
      -sum(log(1/(1+exp(-outX%*%gmaX-as.matrix(Agma*A)))*(p11*Y+(1-p11)*(1-Y))
               +(1-1/(1+exp(-outX%*%gmaX-as.matrix(Agma*A))))*(p10*Y+(1-p10)*(1-Y))))
    } 
    estgma=optim(intc,fc)$par
    Aestgma=estgma[nXout+1]
    estgmaX=estgma[1:nXout]
    m1=1/(1+exp(-outX%*%estgmaX-as.matrix(Aestgma*rep(1,n))))
    m0=1/(1+exp(-outX%*%estgmaX))
    
    mu1=mean(A*Y/(psfit*(p11-p10))-(A-psfit)/(psfit)*m1-
               A/psfit*(p10/(p11-p10)))
    mu0=mean((1-A)*Y/((1-psfit)*(p11-p10))+(A-psfit)/(1-psfit)*m0-
               (1-A)/(1-psfit)*(p10/(p11-p10)))
    mu1-mu0    
  }
}
cor=once(da=data)
STD=sd(replicate(numBoot,once(da=data[sample(1:n,replace=TRUE),])))
low=cor-qnorm(1-(1-confidence)/2)*STD
up=cor+qnorm(1-(1-confidence)/2)*STD

lsthere=list("Estimate"=cor,"Std.Error"=STD)
lsthere[[paste(confidence*100,"% Confidence Interval", sep ="")]] <- c(low,up)
lsthere   
}  
  
  