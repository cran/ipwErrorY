#' Estimation of ATE with Known Error
#'
#' Estimation of average treatment effect with known outcome misclassification probabilities, i.e., known sensitivity and specificity
#'@param data The dataset to be analyzed in the form of R data frame
#'@param indA A column name indicating the treatment variable
#'@param indYerror A column name indicating the misclassified outcome variable
#'@param indX A vector of column names indicating the covariates included in the treatment model
#'@param sensitivity The specified sensitivity between 0 and 1
#'@param specificity The specified specificity between 0 and 1
#'@param confidence The confidence level between 0 and 1; the default is 0.95 corresponding to a 95 per cent confidence interval
#'@return A list of the estimate of average treatment effect, sandwich standard error and confidence interval
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
#'Yast=Y
#'Yast[y1]=rbinom(length(y1),1,0.95)
#'Yast[y0]=rbinom(length(y0),1,0.15)
#'da=data.frame(X1=X1,A=A,Yast=Yast)
#'head(da)
#'#apply the correction method with sensitivity=0.95 and specificity=0.85
#'KnownError(da,"A","Yast","X1",0.95,0.85,0.95)
#'
#'@export
KnownError<-function(data,indA,indYerror,indX,sensitivity,specificity,confidence=0.95){
  p11=sensitivity
  p10=1-specificity 
  n=nrow(data)
  A=data[,indA]
  Y=data[,indYerror]
  nX=length(indX)+1
  trtX0=data[,indX]
  trtX=cbind(rep(1,n),trtX0)
  gg=glm(A~.,family="binomial",data=as.data.frame(cbind(A,trtX0)))
  psfit=predict(gg, type = "response")
  nai=mean(A*Y/psfit)-mean((1-A)*Y/(1-psfit))
  cor=nai/(p11-p10)
  A11=matrix(0,nX,nX)
  B11=A11
  
  for(i in 1:nX){
    for(j in 1:nX){
      A11[i,j]=mean((p11-p10)*psfit*(1-psfit)*trtX[,i]*trtX[,j])
      B11[i,j]=mean((A-psfit)^2*trtX[,i]*trtX[,j])
    }
  }
  
  A22=rep(0,nX)
  B12=A22
  
  for(i in 1:nX){
  A22[i]= mean((A*Y*(1-psfit)/psfit+(1-A)*Y*psfit/(1-psfit))*trtX[,i])
  B12[i]=mean((1-psfit)/psfit*A*Y*trtX[,i]+psfit/(1-psfit)*(1-A)*Y*trtX[,i])
  }
  
  AA=solve(A11)%*%A22
  B22=mean((A*Y/psfit-(1-A)*Y/(1-psfit)-(p11-p10)*cor)^2)
  vars=(t(AA)%*%(B11%*%AA-1/(p11-p10)*B12)-(t(B12)%*%AA-1/(p11-p10)*B22)/(p11-p10))/n
  STD=as.numeric(sqrt(vars))
  low=cor-qnorm(1-(1-confidence)/2)*STD
  up=cor+qnorm(1-(1-confidence)/2)*STD

  lsthere=list("Estimate"=cor,"Std.Error"=STD)
  lsthere[[paste(confidence*100,"% Confidence Interval", sep ="")]] <- c(low,up)
  lsthere   
}

