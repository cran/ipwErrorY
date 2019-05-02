#' Doubly Robust Estimation of ATE with Known Error
#'
#' Doubly robust estimation of average treatment effect with known outcome misclassification probabilities, i.e., known sensitivity and specificity
#'@param data The dataset to be analyzed in the form of R data frame without missing data
#'@param indA A column name indicating the binary treatment variable
#'@param indYerror A column name indicating the misclassified binary outcome variable
#'@param indXtrt A vector of column names indicating the covariates included in the treatment model
#'@param indXout A vector of column names indicating the covariates included in the outcome model
#'@param sensitivity The specified sensitivity between 0 and 1
#'@param specificity The specified specificity between 0 and 1
#'@param sharePara if the treated and untreated groups share parameters for covariates in the logistic outcome model (i.e., assuming Y~ T+X), then set \code{sharePara=TRUE}; if not (i.e., modeling Y~ X for the treated and untreated groups separately), then set \code{sharePara=FALSE}. By default,  \code{sharePara=FALSE}
#'@param confidence The confidence level between 0 and 1; the default is 0.95 corresponding to a 95 per cent confidence interval
#'@return A list of the estimate of average treatment effect, sandwich-variance-based standard error and confidence interval
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
#'KnownErrorDR(da,"A","Yast",c("X","xx"),c("X","xx"),0.95,0.85,FALSE,0.95)
#'
#'@export
KnownErrorDR<-function(data,indA,indYerror,indXtrt,indXout,sensitivity,specificity,sharePara=FALSE,confidence=0.95){
  
  if(sum(apply(data, 2, function(x) as.integer(any(is.na(x)))))>0) stop('invalid dataset with NAs (missing data detected)')
  
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
  outXA=as.matrix(cbind(rep(1,n),outX0,A))
  outXX1=as.matrix(cbind(rep(1,n),outX0,rep(1,n)))
  outXX0=as.matrix(cbind(rep(1,n),outX0,rep(0,n)))
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
  
  ggout=glm(Ytrt~.,family="binomial",data=as.data.frame(cbind(Ytrt,outX0Trt)))
  intc=ggout$coef
    
    fc=function(gma){
      gma=as.matrix(gma)
      -sum(log(1/(1+exp(-outXTrt%*%gma))*(p11*Ytrt+(1-p11)*(1-Ytrt))
          +(1-1/(1+exp(-outXTrt%*%gma)))*(p10*Ytrt+(1-p10)*(1-Ytrt))))
    } 
    estgma=optim(intc,fc)$par
  
  m1=1/(1+exp(-outX%*%estgma))
  gstar1=as.vector(1/(1+exp(-outXTrt%*%estgma)))
  
  ggout=glm(Yout~.,family="binomial",data=as.data.frame(cbind(Yout,outX0Out)))
  intc=ggout$coef
  
  fc=function(gma){
    
    gma=as.matrix(gma)
    -sum(log(1/(1+exp(-outXOut%*%gma))*(p11*Yout+(1-p11)*(1-Yout))
             +(1-1/(1+exp(-outXOut%*%gma)))*(p10*Yout+(1-p10)*(1-Yout))))
  } 
  estgma=optim(intc,fc)$par
  m0=1/(1+exp(-outX%*%estgma))
  gstar0=as.vector(1/(1+exp(-outXOut%*%estgma)))
  
    mu1=mean(A*Y/(psfit*(p11-p10))-(A-psfit)/(psfit)*m1-
               A/psfit*(p10/(p11-p10)))
    mu0=mean((1-A)*Y/((1-psfit)*(p11-p10))+(A-psfit)/(1-psfit)*m0-
               (1-A)/(1-psfit)*(p10/(p11-p10)))
    
    cor=mu1-mu0  
  
  
  A11=matrix(0,nXtrt,nXtrt)
  
  
  for(i in 1:nXtrt){
    for(j in 1:nXtrt){
      A11[i,j]=mean(psfit*(1-psfit)*trtX[,i]*trtX[,j])
    }
  }
  
  A22=matrix(0,nXout,nXout)
  A33=A22
  
  for(i in 1:nXout){
    for(j in 1:nXout){
      A22[i,j]=mean((gstar1*(1-gstar1)*(p11-p10))^2*(Ytrt/((p11*gstar1+p10*(1-gstar1))^2)+(1-Ytrt)/(((1-p11)*gstar1+(1-p10)*(1-gstar1))^2))*outXTrt[,i]*outXTrt[,j]-
                      (p11-p10)*gstar1*(1-gstar1)*(1-2*gstar1)*(Ytrt/(p11*gstar1+p10*(1-gstar1))-(1-Ytrt)/((1-p11)*gstar1+(1-p10)*(1-gstar1)))*outXTrt[,i]*outXTrt[,j])
    
      A33[i,j]=mean((gstar0*(1-gstar0)*(p11-p10))^2*(Yout/((p11*gstar0+p10*(1-gstar0))^2)+(1-Yout)/(((1-p11)*gstar0+(1-p10)*(1-gstar0))^2))*outXOut[,i]*outXOut[,j]-
                      (p11-p10)*gstar0*(1-gstar0)*(1-2*gstar0)*(Yout/(p11*gstar0+p10*(1-gstar0))-(1-Yout)/((1-p11)*gstar0+(1-p10)*(1-gstar0)))*outXOut[,i]*outXOut[,j])
      
    }
  }
  
  A41=rep(0,nXtrt)
  A51=A41
  
  for(i in 1:nXtrt){
    A41[i]=mean((1-psfit)/psfit*(A*Y/(p11-p10)-A*m1-A*p10/(p11-p10))*trtX[,i])
    A51[i]=mean(psfit/(1-psfit)*(-(1-A)*Y/(p11-p10)+(1-A)*m0+(1-A)*p10/(p11-p10))*trtX[,i])
  }
  
  A42=rep(0,nXout)
  A52=A42
  
  for(i in 1:nXout){
    A42[i]=mean((A-psfit)/psfit*m1*(1-m1)*outX[,i])
    A52[i]=mean(-(A-psfit)/(1-psfit)*m0*(1-m0)*outX[,i])
  }
  
  Arow1=cbind(A11,matrix(0,nXtrt,nXout+nXout+3))
  Arow2=cbind(matrix(0,nXout,nXtrt),A22,matrix(0,nXout,nXout+3))
  Arow3=cbind(matrix(0,nXout,nXtrt+nXout),A33,matrix(0,nXout,3))
  Arow4=c(A41,A42,rep(0,nXout),1,0,0)
  Arow5=c(A51,rep(0,nXout),A52,0,1,0)
  Arow6=c(rep(0,nXout+nXout+nXtrt),-1,1,1)
  
  Amat=rbind(Arow1,Arow2,Arow3,Arow4,Arow5,Arow6)
  Ainv=solve(Amat)
  
  BB2=matrix(0,n,nXout)
  BB3=BB2
  
  BB1=diag(A-psfit,n,n)%*%as.matrix(trtX)
  BB2in=diag(gstar1*(1-gstar1)*(p11-p10)*(Ytrt/(p11*gstar1+p10*(1-gstar1))-(1-Ytrt)/((1-p11)*gstar1+(1-p10)*(1-gstar1))),nT,nT)%*%as.matrix(outXTrt)
  BB2[iAi,]=n/nT*BB2in
  BB3in=diag(gstar0*(1-gstar0)*(p11-p10)*(Yout/(p11*gstar0+p10*(1-gstar0))-(1-Yout)/((1-p11)*gstar0+(1-p10)*(1-gstar0))),n-nT,n-nT)%*%as.matrix(outXOut)
  BB3[-iAi,]=n/(n-nT)*BB3in
  BB4=A*Y/(psfit*(p11-p10))-(A-psfit)/(psfit)*m1-A/psfit*(p10/(p11-p10))-mu1
  BB5=(1-A)*Y/((1-psfit)*(p11-p10))+(A-psfit)/(1-psfit)*m0-(1-A)/(1-psfit)*(p10/(p11-p10))-mu0
  BB6=mu1-mu0-cor
  
  BBbind=cbind(BB1,BB2,BB3,BB4,BB5,BB6)
  BB=t(BBbind)%*%BBbind/n
  
  varMat=Ainv%*%BB%*%(t(Ainv))/n
  vars=varMat[nXout+nXout+nXtrt+3,nXout+nXout+nXtrt+3]
  STD=sqrt(vars)
  
  c(cor,STD)
  
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
    cor=mu1-mu0   
    
    gstar=as.vector(1/(1+exp(-outX%*%estgmaX-as.matrix(Aestgma*A))))
    
    A11=matrix(0,nXtrt,nXtrt)
    
  
    for(i in 1:nXtrt){
      for(j in 1:nXtrt){
        A11[i,j]=mean(psfit*(1-psfit)*trtX[,i]*trtX[,j])
      }
    }
    
    A22=matrix(0,nXout+1,nXout+1)
  
    
    for(i in 1:(nXout+1)){
      for(j in 1:(nXout+1)){
        A22[i,j]=mean((gstar*(1-gstar)*(p11-p10))^2*(Y/((p11*gstar+p10*(1-gstar))^2)+(1-Y)/(((1-p11)*gstar+(1-p10)*(1-gstar))^2))*outXA[,i]*outXA[,j]-
          (p11-p10)*gstar*(1-gstar)*(1-2*gstar)*(Y/(p11*gstar+p10*(1-gstar))-(1-Y)/((1-p11)*gstar+(1-p10)*(1-gstar)))*outXA[,i]*outXA[,j])
      }
    }
    
    A31=rep(0,nXtrt)
    A41=A31
    
    for(i in 1:nXtrt){
      A31[i]=mean((1-psfit)/psfit*(A*Y/(p11-p10)-A*m1-A*p10/(p11-p10))*trtX[,i])
      A41[i]=mean(psfit/(1-psfit)*(-(1-A)*Y/(p11-p10)+(1-A)*m0+(1-A)*p10/(p11-p10))*trtX[,i])
    }
    
    A32=rep(0,nXout+1)
    A42=A32
    
    for(i in 1:(nXout+1)){
      A32[i]=mean((A-psfit)/psfit*m1*(1-m1)*outXX1[,i])
      A42[i]=mean(-(A-psfit)/(1-psfit)*m0*(1-m0)*outXX0[,i])
    }
    
    Arow1=cbind(A11,matrix(0,nXtrt,nXout+4))
    Arow2=cbind(matrix(0,nXout+1,nXtrt),A22,matrix(0,nXout+1,3))
    Arow3=c(A31,A32,1,0,0)
    Arow4=c(A41,A42,0,1,0)
    Arow5=c(rep(0,nXout+1+nXtrt),-1,1,1)
    
    Amat=rbind(Arow1,Arow2,Arow3,Arow4,Arow5)
    Ainv=solve(Amat)
    
    BB1=diag(A-psfit,n,n)%*%as.matrix(trtX)
    BB2=diag(gstar*(1-gstar)*(p11-p10)*(Y/(p11*gstar+p10*(1-gstar))-(1-Y)/((1-p11)*gstar+(1-p10)*(1-gstar))),n,n)%*%as.matrix(outXA)
    BB3=A*Y/(psfit*(p11-p10))-(A-psfit)/(psfit)*m1-A/psfit*(p10/(p11-p10))-mu1
    BB4=(1-A)*Y/((1-psfit)*(p11-p10))+(A-psfit)/(1-psfit)*m0-(1-A)/(1-psfit)*(p10/(p11-p10))-mu0
    BB5=mu1-mu0-cor
    
    BBbind=cbind(BB1,BB2,BB3,BB4,BB5)
    BB=t(BBbind)%*%BBbind/n
    
    varMat=Ainv%*%BB%*%(t(Ainv))/n
    vars=varMat[nXout+1+nXtrt+3,nXout+1+nXtrt+3]
    STD=sqrt(vars)
    
    c(cor,STD)
    
  }
}

resu=once(da=data)
cor=resu[1]
STD=resu[2]

low=cor-qnorm(1-(1-confidence)/2)*STD
up=cor+qnorm(1-(1-confidence)/2)*STD

lsthere=list("Estimate"=cor,"Std.Error"=STD)
lsthere[[paste(confidence*100,"% Confidence Interval", sep ="")]] <- c(low,up)
lsthere   
}  
  
  