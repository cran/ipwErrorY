#' Estimation of ATE with Two Replicates 
#'
#' Estimation of average treatment effect when misclassification probabilities are unknown but two independent replicates of the outcome are available
#'@param data The dataset to be analyzed in the form of R data frame without missing data
#'@param indA A column name indicating the binary treatment variable
#'@param indYerror A vector of two column names indicating replicates of the binary outcome variable
#'@param indX A vector of column names indicating the covariates included in the treatment model
#'@param constraint The constraint to be used; the default assumes sensitivity equals specificity
#'@param sensitivity The specified sensitivity between 0 and 1 when imposing the constraint that sensitivity is known, and the default is set to be NULL
#'@param specificity The specified specificity between 0 and 1 when imposing the constraint that specificity is known, and the default is set to be NULL
#'@param prevalence The specified prevalence between 0 and 1 when imposing the constraint that prevalence is known, and the default is set to be NULL
#'@param confidence The confidence level between 0 and 1; the default is 0.95 corresponding to a 95 per cent confidence interval
#'@return A list of the estimate of average treatment effect, sandwich-variance-based standard error, confidence interval, imposed constraint, and the information on sensitivity and specificity
#'@import nleqslv
#'@import stats
#'
#'@examples
#'#create a dataset with sensitivity=0.95 and specificity=0.85
#'set.seed(100)
#'X1=rnorm(2000)   
#'A=rbinom(2000,1,1/(1+exp(-0.2-X1)))
#'Y=rbinom(2000,1,1/(1+exp(-0.2-A-X1)))
#'y1=which(Y==1)
#'y0=which(Y==0) 
#'Yast1=Y
#'Yast1[y1]=rbinom(length(y1),1,0.95)
#'Yast1[y0]=rbinom(length(y0),1,0.15)
#'Yast2=Y
#'Yast2[y1]=rbinom(length(y1),1,0.95)  
#'Yast2[y0]=rbinom(length(y0),1,0.15)
#'da=data.frame(A=A,X1=X1,Yast1=Yast1,Yast2=Yast2)
#'head(da)
#'#apply the correction method assuming specificity=0.85
#'Est2Replicates(da,"A",c("Yast1","Yast2"),"X1","known specificity",NULL,0.85,NULL,0.95)
#'
#'@export
Est2Replicates<-function(data,indA,indYerror,indX, constraint=c("sensitivity equals specificity",
"known sensitivity", "known specificity","known prevalence"),sensitivity=NULL,specificity=NULL,prevalence=NULL,confidence=0.95)
  {
  
  if(sum(apply(data, 2, function(x) as.integer(any(is.na(x)))))>0) stop('invalid dataset with NAs (missing data detected)')
  
  n=nrow(data)
  A=data[,indA]
  YastMX=data[,indYerror]
  Yast1=YastMX[,1]
  Yast2=YastMX[,2]
  Yast=0.5*(Yast1+Yast2)
  nX=length(indX)+1
  trtX0=data[,indX]
  trtX=cbind(rep(1,n),trtX0)
  constraint<-match.arg(constraint)
  gg=glm(A~.,family="binomial",data=as.data.frame(cbind(A,trtX0)))
  ps=predict(gg, type = "response")
  nai=mean(A*Yast/ps)-mean((1-A)*Yast/(1-ps))

  if (constraint=="sensitivity equals specificity"){
  estqu<-function(b){
    b1=b[1]
    b2=b[2]
    yf=rep(0,2)
    yf[1]=mean((1-Yast1)*(1-Yast2)-b1*(1-b2)^2-(1-b1)*b2^2)
    yf[2]=mean(Yast1*(1-Yast2)+Yast2*(1-Yast1)-2*b2*(1-b2))
    yf
  }
  
  outpara=nleqslv(c(0.5,1), estqu)$x
  qq=outpara[1]
  p11=outpara[2]
  sens=p11
  cor=nai/(2*p11-1)
 AA=matrix(0,nX+3,nX+3)
 AA[1,nX+2]=(1-p11)^2-p11^2
 AA[1,nX+3]=2*p11-2*qq
 AA[2,nX+3]=2-4*p11
 AA[nX+3,1]=2*p11-1
 AA[nX+3,nX+3]=2*cor
  
  for(i in 3:(nX+2)){
    for(j in 2:(nX+1)){
      AA[i,j]=mean(ps*(1-ps)*trtX[,i-2]*trtX[,j-1])
    }
  }
  
  for(j in 2:(nX+1)){
    AA[nX+3,j]=mean((A*Yast*(1-ps)/ps+(1-A)*Yast*ps/(1-ps))*trtX[,j-1])
  }
 
 BB=matrix(0,nX+3,nX+3)
 B1=(1-Yast1)*(1-Yast2)-qq*(1-p11)^2-(1-qq)*p11^2
 B2=Yast1*(1-Yast2)+Yast2*(1-Yast1)-2*p11*(1-p11)
 B5=A*Yast/ps-(1-A)*Yast/(1-ps)-(2*p11-1)*cor
 BB[1,1]=mean(B1*B1)
 BB[1,2]=mean(B1*B2)
 BB[1,nX+3]=mean(B1*B5)
 BB[2,1]=mean(B2*B1)
 BB[2,2]=mean(B2*B2)
 BB[2,nX+3]=mean(B2*B5)
 BB[nX+3,1]=mean(B5*B1)
 BB[nX+3,2]=mean(B5*B2)
 BB[nX+3,nX+3]=mean(B5*B5)
 
 for(i in 3:(nX+2)){
   BB[i,1]=mean((A-ps)*trtX[,i-2]*B1)
   BB[i,2]=mean((A-ps)*trtX[,i-2]*B2)
   BB[i,nX+3]=mean((A-ps)*trtX[,i-2]*B5)
 }
 
 for(j in 3:(nX+2)){
   BB[1,j]=mean((A-ps)*trtX[,j-2]*B1)
   BB[2,j]=mean((A-ps)*trtX[,j-2]*B2)
   BB[nX+3,j]=mean((A-ps)*trtX[,j-2]*B5)
 }
 
 for(i in 3:(nX+2)){
   for(j in 3:(nX+2)){
     BB[i,j]=mean((A-ps)*trtX[,i-2]*(A-ps)*trtX[,j-2])
   }
 }
 
  mmm=((solve(AA)%*%BB)%*%t(solve(AA)))/n
  Vtau=mmm[1,1] 
  STD=sqrt(Vtau)
 low=cor-qnorm(1-(1-confidence)/2)*STD
 up=cor+qnorm(1-(1-confidence)/2)*STD
 
 outputHere=list("Estimate"=cor,"Std.Error"=STD)
 outputHere[[paste(confidence*100,"% Confidence Interval", sep ="")]] <- c(low,up)
 outputHere[["imposed constraint"]] <- "sensitivity equals specificity"
 outputHere[["estimated sensitivity and estimated specificity"]] <- c(sens,sens)

  }

 else if (constraint=="known sensitivity"){
 p11=sensitivity
 sens=p11
 estqu<-function(b){
   b1=b[1]
   b2=b[2]
   yf=rep(0,2)
   yf[1]=mean((1-Yast1)*(1-Yast2)-b1*(1-p11)^2-(1-b1)*(1-b2)^2)
   yf[2]=mean(Yast1*(1-Yast2)+Yast2*(1-Yast1)-2*b1*(1-p11)*p11-2*(1-b1)*b2*(1-b2))
   yf
 }
 
 outpara=nleqslv(c(0.5,0), estqu)$x
 qq=outpara[1]
 p10=outpara[2]
 speci=1-p10
 cor=nai/(p11-p10)
 AA=matrix(0,nX+3,nX+3)
 AA[1,2]=(1-p11)^2-(1-p10)^2 
 AA[1,3]=-2*(1-qq)*(1-2*p10) 
 AA[2,2]=2*(1-p11)*p11-2*(1-p10)*p10
 AA[2,3]=2*(1-qq)*(1-2*p10)
 AA[nX+3,1]=p11-p10
 AA[nX+3,3]=-cor
 
 for(i in 3:(nX+2)){
   for(j in 4:(nX+3)){
     AA[i,j]=mean(ps*(1-ps)*trtX[,i-2]*trtX[,j-3])
   }
 }
 
 for(j in 4:(nX+3)){
   AA[nX+3,j]=mean((A*Yast*(1-ps)/ps+(1-A)*Yast*ps/(1-ps))*trtX[,j-3])
 }
 
 BB=matrix(0,nX+3,nX+3)
 B1=(1-Yast1)*(1-Yast2)-qq*(1-p11)^2-(1-qq)*(1-p10)^2
 B2=Yast1*(1-Yast2)+Yast2*(1-Yast1)-2*qq*(1-p11)*p11-2*(1-qq)*(1-p10)*p10
 B5=A*Yast/ps-(1-A)*Yast/(1-ps)-(p11-p10)*cor
 BB[1,1]=mean(B1*B1)
 BB[1,2]=mean(B1*B2)
 BB[1,nX+3]=mean(B1*B5)
 BB[2,1]=mean(B2*B1)
 BB[2,2]=mean(B2*B2)
 BB[2,nX+3]=mean(B2*B5)
 BB[nX+3,1]=mean(B5*B1)
 BB[nX+3,2]=mean(B5*B2)
 BB[nX+3,nX+3]=mean(B5*B5)
 
 for(i in 3:(nX+2)){
   BB[i,1]=mean((A-ps)*trtX[,i-2]*B1)
   BB[i,2]=mean((A-ps)*trtX[,i-2]*B2)
   BB[i,nX+3]=mean((A-ps)*trtX[,i-2]*B5)
 }
 
 for(j in 3:(nX+2)){
   BB[1,j]=mean((A-ps)*trtX[,j-2]*B1)
   BB[2,j]=mean((A-ps)*trtX[,j-2]*B2)
   BB[nX+3,j]=mean((A-ps)*trtX[,j-2]*B5)
 }
 
 for(i in 3:(nX+2)){
   for(j in 3:(nX+2)){
     BB[i,j]=mean((A-ps)*trtX[,i-2]*(A-ps)*trtX[,j-2])
   }
 }
 
 mmm=((solve(AA)%*%BB)%*%t(solve(AA)))/n
 
 Vtau=mmm[1,1] 
 STD=sqrt(Vtau)


 low=cor-qnorm(1-(1-confidence)/2)*STD
 up=cor+qnorm(1-(1-confidence)/2)*STD
 
 outputHere=list("Estimate"=cor,"Std.Error"=STD)
 outputHere[[paste(confidence*100,"% Confidence Interval", sep ="")]] <- c(low,up)
 outputHere[["imposed constraint"]] <- "known sensitivity"
 outputHere[["assumed sensitivity and estimated specificity"]] <- c(sens,speci)
 
 }
 
 else if (constraint=="known specificity"){
 p00=specificity
 p10=1-p00
 speci=p00
 estqu<-function(b){
   b1=b[1]
   b2=b[2]
   yf=rep(0,2)
   yf[1]=mean((1-Yast1)*(1-Yast2)-b1*(1-b2)^2-(1-b1)*(1-p10)^2)
   yf[2]=mean(Yast1*(1-Yast2)+Yast2*(1-Yast1)-2*b1*(1-b2)*b2-2*(1-b1)*p10*(1-p10))
   yf
 }
 
 outpara=nleqslv(c(0.5,1), estqu)$x
 qq=outpara[1]
 p11=outpara[2]
 sens=p11
 cor=nai/(p11-p10)
 AA=matrix(0,nX+3,nX+3)
 AA[1,2]=(1-p11)^2-(1-p10)^2 
 AA[1,3]=-2*qq*(1-p11)
 AA[2,2]=2*(1-p11)*p11-2*(1-p10)*p10
 AA[2,3]=2*qq*(1-2*p11)
 AA[nX+3,1]=p11-p10
 AA[nX+3,3]=cor
 
 for(i in 3:(nX+2)){
   for(j in 4:(nX+3)){
     AA[i,j]=mean(ps*(1-ps)*trtX[,i-2]*trtX[,j-3])
   }
 }
 
 for(j in 4:(nX+3)){
   AA[nX+3,j]=mean((A*Yast*(1-ps)/ps+(1-A)*Yast*ps/(1-ps))*trtX[,j-3])
 }
 
 BB=matrix(0,nX+3,nX+3)
 B1=(1-Yast1)*(1-Yast2)-qq*(1-p11)^2-(1-qq)*(1-p10)^2
 B2=Yast1*(1-Yast2)+Yast2*(1-Yast1)-2*qq*(1-p11)*p11-2*(1-qq)*(1-p10)*p10
 B5=A*Yast/ps-(1-A)*Yast/(1-ps)-(p11-p10)*cor
 BB[1,1]=mean(B1*B1)
 BB[1,2]=mean(B1*B2)
 BB[1,nX+3]=mean(B1*B5)
 BB[2,1]=mean(B2*B1)
 BB[2,2]=mean(B2*B2)
 BB[2,nX+3]=mean(B2*B5)
 BB[nX+3,1]=mean(B5*B1)
 BB[nX+3,2]=mean(B5*B2)
 BB[nX+3,nX+3]=mean(B5*B5)
 
 for(i in 3:(nX+2)){
   BB[i,1]=mean((A-ps)*trtX[,i-2]*B1)
   BB[i,2]=mean((A-ps)*trtX[,i-2]*B2)
   BB[i,nX+3]=mean((A-ps)*trtX[,i-2]*B5)
 }
 
 for(j in 3:(nX+2)){
   BB[1,j]=mean((A-ps)*trtX[,j-2]*B1)
   BB[2,j]=mean((A-ps)*trtX[,j-2]*B2)
   BB[nX+3,j]=mean((A-ps)*trtX[,j-2]*B5)
 }
 
 for(i in 3:(nX+2)){
   for(j in 3:(nX+2)){
     BB[i,j]=mean((A-ps)*trtX[,i-2]*(A-ps)*trtX[,j-2])
   }
 }
 
 mmm=((solve(AA)%*%BB)%*%t(solve(AA)))/n
 Vtau=mmm[1,1] 
 STD=sqrt(Vtau)

 low=cor-qnorm(1-(1-confidence)/2)*STD
 up=cor+qnorm(1-(1-confidence)/2)*STD
 
 outputHere=list("Estimate"=cor,"Std.Error"=STD)
 outputHere[[paste(confidence*100,"% Confidence Interval", sep ="")]] <- c(low,up)
 outputHere[["imposed constraint"]] <- "known specificity"
 outputHere[["estimated sensitivity and assumed specificity"]] <- c(sens,speci)
 
 
 }

 else if (constraint=="known prevalence"){
qq=prevalence
 
 estqu<-function(b){
   b1=b[1]
   b2=b[2]
   yf=rep(0,2)
   yf[1]=mean((1-Yast1)*(1-Yast2)-qq*(1-b1)^2-(1-qq)*(1-b2)^2)
   yf[2]=mean(Yast1*(1-Yast2)+Yast2*(1-Yast1)-2*qq*(1-b1)*b1-2*(1-qq)*(1-b2)*b2)
   yf
 }
 
 outpara=nleqslv(c(1,0), estqu)$x
 p11=outpara[1]
 p10=outpara[2]
sens=p11
speci=1-p10
 cor=nai/(p11-p10)
 AA=matrix(0,nX+3,nX+3)
 AA[1,2]=-2*qq*(1-p11)
 AA[1,3]=-2*(1-qq)*(1-p10)
 AA[2,2]=2*qq*(1-2*p11)
 AA[2,3]=2*(1-qq)*(1-2*p10)
 AA[nX+3,1]=p11-p10
 AA[nX+3,2]=cor
 AA[nX+3,3]=-cor
 
for(i in 3:(nX+2)){
   for(j in 4:(nX+3)){
     AA[i,j]=mean(ps*(1-ps)*trtX[,i-2]*trtX[,j-3])
   }
 }
 
 for(j in 4:(nX+3)){
   AA[nX+3,j]=mean((A*Yast*(1-ps)/ps+(1-A)*Yast*ps/(1-ps))*trtX[,j-3])
 }
 
 BB=matrix(0,nX+3,nX+3)
 B1=(1-Yast1)*(1-Yast2)-qq*(1-p11)^2-(1-qq)*(1-p10)^2
 B2=Yast1*(1-Yast2)+Yast2*(1-Yast1)-2*qq*(1-p11)*p11-2*(1-qq)*(1-p10)*p10
 B5=A*Yast/ps-(1-A)*Yast/(1-ps)-(p11-p10)*cor
 BB[1,1]=mean(B1*B1)
 BB[1,2]=mean(B1*B2)
 BB[1,nX+3]=mean(B1*B5)
 BB[2,1]=mean(B2*B1)
 BB[2,2]=mean(B2*B2)
 BB[2,nX+3]=mean(B2*B5)
 BB[nX+3,1]=mean(B5*B1)
 BB[nX+3,2]=mean(B5*B2)
 BB[nX+3,nX+3]=mean(B5*B5)
 
 for(i in 3:(nX+2)){
   BB[i,1]=mean((A-ps)*trtX[,i-2]*B1)
   BB[i,2]=mean((A-ps)*trtX[,i-2]*B2)
   BB[i,nX+3]=mean((A-ps)*trtX[,i-2]*B5)
 }
 
 for(j in 3:(nX+2)){
   BB[1,j]=mean((A-ps)*trtX[,j-2]*B1)
   BB[2,j]=mean((A-ps)*trtX[,j-2]*B2)
   BB[nX+3,j]=mean((A-ps)*trtX[,j-2]*B5)
 }
 
 for(i in 3:(nX+2)){
   for(j in 3:(nX+2)){
     BB[i,j]=mean((A-ps)*trtX[,i-2]*(A-ps)*trtX[,j-2])
   }
 }
 
 mmm=((solve(AA)%*%BB)%*%t(solve(AA)))/n
 Vtau=mmm[1,1] 
 STD=sqrt(Vtau)
 
low=cor-qnorm(1-(1-confidence)/2)*STD
up=cor+qnorm(1-(1-confidence)/2)*STD

outputHere=list("Estimate"=cor,"Std.Error"=STD)
outputHere[[paste(confidence*100,"% Confidence Interval", sep ="")]] <- c(low,up)
outputHere[["imposed constraint"]] <- "known prevalence"
outputHere[["estimated sensitivity and estimated specificity"]] <- c(sens,speci)

 }
 
outputHere
}
