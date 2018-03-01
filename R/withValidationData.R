#' Estimation of ATE with Validation Data
#'
#' Estimation of average treatment effect using the optimal linear combination method when misclassification probabilities are unknown but validation data are available
#'@param maindata The non-validation main data in the form of R data frame
#'@param validationdata The validation data in the form of R data frame
#'@param indA A column name indicating the treatment variable
#'@param indYerror A column name indicating the misclassified outcome variable
#'@param indX A vector of column names indicating the covariates included in the treatment model
#'@param indY A column name indicating the true outcome variable 
#'@param confidence The confidence level between 0 and 1; the default is 0.95 corresponding to a 95 per cent confidence interval
#'@return A list of the estimate of average treatment effect, sandwich standard error, confidence interval, and the estimated sensitivity and specificity
#'@import stats
#'
#'@examples
#'#create main data and validation data with sensitivity=0.95 and specificity=0.85
#'set.seed(100)
#'X1=rnorm(1200)   
#'A=rbinom(1200,1,1/(1+exp(-0.2-X1)))
#'Y=rbinom(1200,1,1/(1+exp(-0.2-A-X1)))
#'y1=which(Y==1)
#'y0=which(Y==0) 
#'Yast=Y
#'Yast[y1]=rbinom(length(y1),1,0.95)
#'Yast[y0]=rbinom(length(y0),1,0.15)
#'mainda=data.frame(A=A,X1=X1,Yast=Yast)
#'X1=rnorm(800)   
#'A=rbinom(800,1,1/(1+exp(-0.2-X1)))
#'Y=rbinom(800,1,1/(1+exp(-0.2-A-X1)))
#'y1=which(Y==1)
#'y0=which(Y==0) 
#'Yast=Y
#'Yast[y1]=rbinom(length(y1),1,0.95)
#'Yast[y0]=rbinom(length(y0),1,0.15)
#'validationda=data.frame(A=A,X1=X1,Y=Y,Yast=Yast)
#'head(mainda)
#'head(validationda)
#'#apply the optimal linear combination correction method
#'EstValidation(mainda,validationda,"A","Yast","X1","Y",0.95)
#'
#'@export
EstValidation<-function(maindata,validationdata,indA,indYerror,indX,indY,confidence=0.95){
  nv=nrow(validationdata)
  n=nrow(maindata)+nv
  AM=maindata[,indA]
  YastM=maindata[,indYerror]
  AV=validationdata[,indA]
  YastV=validationdata[,indYerror]
  nX=length(indX)+1
  trtX0main=maindata[,indX]
  trtX0vali=validationdata[,indX]
  trtXmain=cbind(rep(1,n-nv),trtX0main)
  trtXvali=cbind(rep(1,nv),trtX0vali)
  
  A=c(AV,AM)
   
  if (nX>2) {
    trtX0=rbind(trtX0vali,trtX0main)
  } else {
    trtX0=c(trtX0vali,trtX0main)
  }
  
  trtX=cbind(rep(1,n),trtX0)
  gg=glm(A~.,family="binomial",data=as.data.frame(cbind(A,trtX0)))
  ps=predict(gg, type = "response")
  vps=ps[1:nv]
  mps=ps[(nv+1):n]
  YV=validationdata[,indY]
  Ynai=c(YV,YastM)
  psfit=ps
  nai=mean(A*Ynai/psfit)-mean((1-A)*Ynai/(1-psfit))
  p11=sum(YV*YastV)/sum(YV)
  p10=sum((1-YV)*YastV)/sum(1-YV)
  sens=p11
  speci=1-p10
  tauV=mean(AV*YV/vps)-mean((1-AV)*YV/(1-vps))
  tauN=(mean(AM*YastM/mps)-mean((1-AM)*YastM/(1-mps)))/(p11-p10)
  AA=matrix(0,nX+4,nX+4)
  AA[nX+1,1]=p11-p10
  AA[nX+1,nX+3]=tauN
  AA[nX+1,nX+4]=-tauN
  AA[nX+2,nX+3]=mean(YV)
  AA[nX+3,nX+4]=mean(1-YV)
  AA[nX+4,2]=1
  
  for(i in 1:nX){
    for(j in 3:(nX+2)){
      AA[i,j]=mean(ps*(1-ps)*trtX[,i]*trtX[,j-2])
    }
  }
  
  for(j in 3:(nX+2)){
    AA[nX+1,j]=mean((AM*YastM*(1-mps)/mps+(1-AM)*YastM*mps/(1-mps))*trtXmain[,j-2])
    AA[nX+4,j]=mean((AV*YV*(1-vps)/vps+(1-AV)*YV*vps/(1-vps))*trtXvali[,j-2])
  }
  
  BB=matrix(0,nX+4,nX+4)
  B3=c(rep(0,nv),AM*YastM/mps-(1-AM)*YastM/(1-mps)-(p11-p10)*tauN)*n/(n-nv)
  B4=c(YV*YastV-p11*YV,rep(0,n-nv))*n/nv
  B5=c((1-YV)*YastV-p10*(1-YV), rep(0,n-nv))*n/nv
  B6=c(AV*YV/vps-(1-AV)*YV/(1-vps)-tauV,rep(0,n-nv))*n/nv
  BB[nX+1,nX+1]=mean(B3*B3)
  BB[nX+1,nX+2]=mean(B3*B4)
  BB[nX+1,nX+3]=mean(B3*B5)
  BB[nX+1,nX+4]=mean(B3*B6)
  BB[nX+2,nX+1]=mean(B4*B3)
  BB[nX+2,nX+2]=mean(B4*B4)
  BB[nX+2,nX+3]=mean(B4*B5)
  BB[nX+2,nX+4]=mean(B4*B6)
  BB[nX+3,nX+1]=mean(B5*B3)
  BB[nX+3,nX+2]=mean(B5*B4)
  BB[nX+3,nX+3]=mean(B5*B5)
  BB[nX+3,nX+4]=mean(B5*B6)
  BB[nX+4,nX+1]=mean(B6*B3)
  BB[nX+4,nX+2]=mean(B6*B4)
  BB[nX+4,nX+3]=mean(B6*B5)
  BB[nX+4,nX+4]=mean(B6*B6)
  
  for(i in 1:nX){
    BB[i,nX+1]=mean((A-ps)*trtX[,i]*B3)
    BB[i,nX+2]=mean((A-ps)*trtX[,i]*B4)
    BB[i,nX+3]=mean((A-ps)*trtX[,i]*B5)
    BB[i,nX+4]=mean((A-ps)*trtX[,i]*B6)
  }
  
  for(j in 1:nX){
    BB[nX+1,j]=mean((A-ps)*trtX[,j]*B3)
    BB[nX+2,j]=mean((A-ps)*trtX[,j]*B4)
    BB[nX+3,j]=mean((A-ps)*trtX[,j]*B5)
    BB[nX+4,j]=mean((A-ps)*trtX[,j]*B6)
  }
  
  for(i in 1:nX){
    for(j in 1:nX){
      BB[i,j]=mean((A-ps)*trtX[,i]*(A-ps)*trtX[,j])
    }
  }
  
  mmmw=((solve(AA)%*%BB)%*%t(solve(AA)))/n
  VtauN=mmmw[1,1]
  Cov2=mmmw[1,2]
  VtauV2=mmmw[2,2]
  nSTD=sqrt(VtauN)
  nlow=tauN-1.96*nSTD  
  nup=tauN+1.96*nSTD
  vSTD=sqrt(VtauV2)
  vlow=tauV-1.96*vSTD
  vup=tauV+1.96*vSTD 
  ORal2=ifelse(((VtauN-Cov2)/(VtauN+VtauV2-2*Cov2))*(1-(VtauN-Cov2)/(VtauN+VtauV2-2*Cov2))>=0
               &VtauN+VtauV2-2*Cov2>0,
               (VtauN-Cov2)/(VtauN+VtauV2-2*Cov2),as.integer(VtauV2<VtauN))
  tauW=ORal2*tauV+(1-ORal2)*tauN
  VtauW=ORal2^2*VtauV2+(1-ORal2)^2*VtauN+2*ORal2*(1-ORal2)*Cov2
  wSTD=sqrt(VtauW)
  
  wlow=tauW-qnorm(1-(1-confidence)/2)*wSTD
  wup=tauW+qnorm(1-(1-confidence)/2)*wSTD
  
  ORal3=nv/n
  tau12=ORal3*tauV+(1-ORal3)*tauN
  Vtau12=ORal3^2*VtauV2+(1-ORal3)^2*VtauN+2*ORal3*(1-ORal3)*Cov2
  STDs=sqrt(Vtau12)
  
  EST=c(tauV,tauN,tau12,tauW)
  VARS=c(VtauV2,VtauN,Vtau12,VtauW)
  STD=VARS^0.5


  lsthere=list("Estimate"=tauW,"Std.Error"=wSTD)
  lsthere[[paste(confidence*100,"% Confidence Interval", sep ="")]] <- c(wlow,wup)
  lsthere[["estimated sensitivity and estimated specificity"]] <- c(sens,speci)
  lsthere     
}

