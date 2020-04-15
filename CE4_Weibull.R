library(rootSolve) # need this package for the median.surv function
library(survival)
library(eha)
library(msm)

weib <- function(data,varlist=0,age.mean=0,SevScaleBL.mean=0,smoke=0){
  # data contain the following columns: Y, Delta, Trt, M
  # Trt=1: treatment; Trt=0: control
  #p1=p.Aa/(p.Aa+p.aa)
  #p2=p.AA/(p.AA+p.Aa)
  p1 <- nrow(data[which(data$M==1),])/(nrow(data[which(data$M==1),])+nrow(data[which(data$M==2),]))
  p2 <- nrow(data[which(data$M==0),])/(nrow(data[which(data$M==1),])+nrow(data[which(data$M==0),]))
  if (varlist==0){
    fit.aftweib <- aftreg(Surv(Y,Delta) ~ as.factor(Trt)+as.factor(M)+as.factor(Trt)*as.factor(M), data=data)
    fit.LogCoeff <- fit.aftweib$coef
    fit.LogCoeff <- fit.LogCoeff[c(6,7,1:5)]
    fit.LogCoeff[3:7] <- fit.LogCoeff[3:7]*exp(fit.LogCoeff[2])
    fit.Coeff <- exp(fit.LogCoeff) # the coefficient in the "original" scale (exponentiated from the log scale)
    # change the order of the parameters (so that they correspond to "lambda","k","theta1","theta2","theta3")
    names(fit.Coeff)[1:2]<-c("Scale","Shape")
    ##deltamethod on the var-cov from aftweib
    fit.LogVar <- deltamethod(list(~exp(x7)*x1,~exp(x7)*x2,~exp(x7)*x3,~exp(x7)*x4,~exp(x7)*x5,~x6,~x7),mean=fit.aftweib$coef,cov=fit.aftweib$var,ses=FALSE)
    # change the order of the var parameters accordingly 
    tmp1 <- fit.LogVar[c(1:5),c(1:5)]
    tmp2 <- fit.LogVar[6:7,6:7]
    tmp3 <- fit.LogVar[c(1:5),6:7]
    
    fit.LogVar.2 <- rbind(cbind(tmp2,t(tmp3)),cbind(tmp3,tmp1))
    # delta method to compute the Variance matrix for the parameter estimates in the "original" scale
    fit.Var <- deltamethod(list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),~exp(x6),~exp(x7)),mean=fit.LogCoeff,cov=fit.LogVar.2,ses=FALSE)
  } else{
    fit.aftweib <- aftreg(Surv(Y,Delta) ~ as.factor(Trt)+as.factor(M)+as.factor(Trt)*as.factor(M)+enrollage+as.factor(smoke)+SevScaleBL, data=data)
    fit.LogCoeff <- fit.aftweib$coef
    fit.LogCoeff <- fit.LogCoeff[c(9,10,1:3,7:8,4:6)]
    fit.LogCoeff[3:10] <- fit.LogCoeff[3:10]*exp(fit.LogCoeff[2])
    fit.Coeff <- exp(fit.LogCoeff) # the coefficient in the "original" scale (exponentiated from the log scale)
    # change the order of the parameters (so that they correspond to "lambda","k","theta1","theta2","theta3")
    names(fit.Coeff)[1:2]<-c("Scale","Shape")
    # change the order of the var parameters accordingly 
    fit.LogVar <- deltamethod(list(~exp(x10)*x1,~exp(x10)*x2,~exp(x10)*x3,~exp(x10)*x4,~exp(x10)*x5,~exp(x10)*x6,~exp(x10)*x7,~exp(x10)*x8,~x9,~x10),mean=fit.aftweib$coef,cov=fit.aftweib$var,ses=FALSE)
    # change the order of the var parameters accordingly 
    tmp1 <- fit.LogVar[c(1:3,7:8,4:6),c(1:3,7:8,4:6)]
    tmp2 <- fit.LogVar[9:10,9:10]
    tmp3 <- fit.LogVar[c(1:3,7:8,4:6),9:10]
    
    fit.LogVar.2 <- rbind(cbind(tmp2,t(tmp3)),cbind(tmp3,tmp1))
    # delta method to compute the Variance matrix for the parameter estimates in the "original" scale
    fit.Var <- deltamethod(list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),~exp(x6),~exp(x7),~exp(x8*age.mean+x9*smoke+x10*SevScaleBL.mean)),mean=fit.LogCoeff,cov=fit.LogVar.2,ses=FALSE)
  }
  
  return(list(fit.aftweib=fit.aftweib,p1=p1,p2=p2,log.coef=fit.LogCoeff,log.var=fit.LogVar,coef=fit.Coeff,var=fit.Var))
}


# calculate the median survival time for each subgroup and the combined groups
median.surv <- function(coef,p1,p2,varlist=0,tau=0.5){
  lambda <- coef[1]
  k <- coef[2]
  if (varlist==0) {
    theta<-c(as.numeric(coef[3]),as.numeric(coef[4]),as.numeric(coef[5]),as.numeric(coef[6]),as.numeric(coef[7]),1)
  } else{
    theta<-c(as.numeric(coef[3]),as.numeric(coef[4]),as.numeric(coef[5]),as.numeric(coef[6]),as.numeric(coef[7]),as.numeric(exp(log(coef[8])*age.mean+log(coef[9])*0+log(coef[10])*SevScaleBL.mean)))
  }
  # lambda: Weibull distribution's scale parameter 
  # k: Weibull distribution's shape parameter 
  # theta=(theta1,theta2,theta3,theta4,theta5,theta6): =exp(beta) where beta is the slope parameter in the Cox PH model
  
  median.C.AA <- as.numeric(lambda*(log(1/tau)/theta[6])^(1/k))
  median.C.Aa <- as.numeric(lambda*(log(1/tau)/(theta[2]*theta[6]))^(1/k))
  median.C.aa <- as.numeric(lambda*(log(1/tau)/(theta[3]*theta[6]))^(1/k))
  median.Rx.AA <- as.numeric(lambda*(log(1/tau)/(theta[1]*theta[6]))^(1/k))  
  median.Rx.Aa <- as.numeric(lambda*(log(1/tau)/(theta[1]*theta[2]*theta[4]*theta[6]))^(1/k))
  median.Rx.aa <- as.numeric(lambda*(log(1/tau)/(theta[1]*theta[3]*theta[5]*theta[6]))^(1/k))
  
  f.Rx12 <- function(t){
    f=p1*exp(-theta[1]*theta[2]*theta[4]*theta[6]*(t/lambda)^k)+(1-p1)*exp(-theta[1]*theta[3]*theta[5]*theta[6]*(t/lambda)^k)-tau
    return(f)
  }
  f.C12 <- function(t){
    f=p1*exp(-theta[2]*theta[6]*(t/lambda)^k)+(1-p1)*exp(-theta[3]*theta[6]*(t/lambda)^k)-tau
    return(f)
  }
  f.Rx10 <- function(t){
    f=p2*exp(-theta[1]*theta[6]*(t/lambda)^k)+(1-p2)*exp(-theta[1]*theta[2]*theta[4]*theta[6]*(t/lambda)^k)-tau
    return(f)
  }
  f.C10 <- function(t){
    f=p2*exp(-theta[6]*(t/lambda)^k)+(1-p2)*exp(-theta[2]*theta[6]*(t/lambda)^k)-tau
    return(f)
  }
  median.Rx12 <- uniroot(f.Rx12,c(min(median.Rx.Aa,median.Rx.aa),max(median.Rx.Aa,median.Rx.aa)))$root
  median.C12 <- uniroot(f.C12,c(min(median.C.Aa,median.C.aa),max(median.C.Aa,median.C.aa)))$root
  median.Rx10 <- uniroot(f.Rx10,c(min(median.Rx.Aa,median.Rx.AA),max(median.Rx.Aa,median.Rx.AA)))$root
  median.C10 <- uniroot(f.C10,c(min(median.C.Aa,median.C.AA),max(median.C.Aa,median.C.AA)))$root
  
  return(list(median.Rx.AA=median.Rx.AA,median.C.AA=median.C.AA,
              median.Rx.Aa=median.Rx.Aa,median.C.Aa=median.C.Aa,
              median.Rx.aa=median.Rx.aa,median.C.aa=median.C.aa,
              median.Rx12=median.Rx12,median.C12=median.C12,
              median.Rx10=median.Rx10,median.C10=median.C10,
              r0=log(median.Rx.AA/median.C.AA),
              r1=log(median.Rx.Aa/median.C.Aa),
              r2=log(median.Rx.aa/median.C.aa),
              r12=log(median.Rx12/median.C12),
              r10=log(median.Rx10/median.C10)))
  
}
# check2<-median.surv(coef=check$coef,p1=check$p1,p2=check$p2,age=0)
# check3<-median.surv(coef=check1$coef,p1=check$p1,p2=check$p2,age=40)

# use the delta method for implicitly defined RVs
# reference paper: Jacques Benichou and Mitchell Gail (The American Statistician, 1989, Vol.43, No.1) 

var.median.surv.dt.dim7 <- function(median.result,p1,p2,coef,param.var,varlist=0,alpha){
  median.Rx12=median.result$median.Rx12
  median.C12=median.result$median.C12
  median.Rx10=median.result$median.Rx10
  median.C10=median.result$median.C10
  lambda <- coef[1]
  k <- coef[2]
  if (varlist==0) {
    theta<-c(as.numeric(coef[3]),as.numeric(coef[4]),as.numeric(coef[5]),as.numeric(coef[6]),as.numeric(coef[7]),1)
    param.var<-cbind(param.var,rep(0,7))
    param.var<-rbind(param.var,rep(0,8))
  } else{
    theta<-c(as.numeric(coef[3]),as.numeric(coef[4]),as.numeric(coef[5]),as.numeric(coef[6]),as.numeric(coef[7]),as.numeric(exp(log(coef[8])*age.mean+log(coef[9])*0+log(coef[10])*SevScaleBL.mean)))
  }
  
  # calculate the variance estimate for the median survival time estimate of all group
  # using delta method for implicitly defined random variables
  
  J.11<-J.22<-J.33<--1
  J.44<-(-k*theta[1]*theta[2]*theta[4]*theta[6]*p1*(median.Rx12/lambda)^(k-1)/lambda*exp(-theta[1]*theta[2]*theta[4]*theta[6]*(median.Rx12/lambda)^k))-
    (k*theta[1]*theta[3]*theta[5]*theta[6]*(1-p1)/lambda)*(median.Rx12/lambda)^(k-1)*exp(-theta[1]*theta[3]*theta[5]*theta[6]*(median.Rx12/lambda)^k)
  J.55<-(-k*theta[2]*theta[6]*p1*(median.C12/lambda)^(k-1)/lambda*exp(-theta[2]*theta[6]*(median.C12/lambda)^k))-
    (k*theta[3]*theta[6]*(1-p1)/lambda*(median.C12/lambda)^(k-1)*exp(-theta[3]*theta[6]*(median.C12/lambda)^k))
  J.66<-(-k*theta[1]*theta[6]*p2*(median.Rx10/lambda)^(k-1)/lambda*exp(-theta[1]*theta[6]*(median.Rx10/lambda)^k))-
    (k*theta[1]*theta[2]*theta[4]*theta[6]*(1-p2)/lambda*(median.Rx10/lambda)^(k-1)*exp(-theta[1]*theta[2]*theta[4]*theta[6]*(median.Rx10/lambda)^k))
  J.77<-(-k*theta[6]*p2*(median.C10/lambda)^(k-1)/lambda*exp(-theta[6]*(median.C10/lambda)^k))-
    (k*theta[2]*theta[6]*(1-p2)/lambda*(median.C10/lambda)^(k-1)*exp(-theta[2]*theta[6]*(median.C10/lambda)^k))
  
  J.mat <- diag(c(J.11, J.22, J.33, J.44, J.55, J.66, J.77))  
  J.mat.inv <- solve(J.mat)
  
  H.mat <- mat.or.vec(7,8)
  H.mat[1,] <- c(0,log(theta[1])/(k^2),-1/(k*theta[1]),0,0,0,0,0)
  H.mat[2,] <- c(0,log(theta[1]*theta[4])/(k^2),-1/(k*theta[1]),0,0,-1/(k*theta[4]),0,0)
  H.mat[3,] <- c(0,log(theta[1]*theta[5])/(k^2),-1/(k*theta[1]),0,0,0,-1/(k*theta[5]),0)
  H.mat[4,] <- c(p1*theta[1]*theta[2]*theta[4]*theta[6]*exp(-theta[1]*theta[2]*theta[4]*theta[6]*(median.Rx12/lambda)^k)*k*(median.Rx12/lambda)^(k-1)*(median.Rx12/(lambda^2))+
                   (1-p1)*theta[1]*theta[3]*theta[5]*theta[6]*exp(-theta[1]*theta[3]*theta[5]*theta[6]*(median.Rx12/lambda)^k)*k*(median.Rx12/lambda)^(k-1)*(median.Rx12/(lambda^2)), 
                 -p1*theta[1]*theta[2]*theta[4]*theta[6]*exp(-theta[1]*theta[2]*theta[4]*theta[6]*(median.Rx12/lambda)^k)*(median.Rx12/lambda)^k*log(median.Rx12/lambda)-
                   (1-p1)*theta[1]*theta[3]*theta[5]*theta[6]*exp(-theta[1]*theta[3]*theta[5]*theta[6]*(median.Rx12/lambda)^k)*(median.Rx12/lambda)^k*log(median.Rx12/lambda), 
                 -p1*exp(-theta[1]*theta[2]*theta[4]*theta[6]*(median.Rx12/lambda)^k)*theta[2]*theta[4]*theta[6]*(median.Rx12/lambda)^k-
                   (1-p1)*exp(-theta[1]*theta[3]*theta[5]*theta[6]*(median.Rx12/lambda)^k)*theta[3]*theta[5]*theta[6]*(median.Rx12/lambda)^k,
                 -p1*exp(-theta[1]*theta[2]*theta[4]*theta[6]*(median.Rx12/lambda)^k)*theta[1]*theta[4]*theta[6]*(median.Rx12/lambda)^k,
                 -(1-p1)*exp(-theta[1]*theta[3]*theta[5]*theta[6]*(median.Rx12/lambda)^k)*theta[1]*theta[5]*theta[6]*(median.Rx12/lambda)^k,
                 -p1*exp(-theta[1]*theta[2]*theta[4]*theta[6]*(median.Rx12/lambda)^k)*theta[1]*theta[2]*theta[6]*(median.Rx12/lambda)^k,
                 -(1-p1)*exp(-theta[1]*theta[3]*theta[5]*theta[6]*(median.Rx12/lambda)^k)*theta[1]*theta[3]*theta[6]*(median.Rx12/lambda)^k,
                 -p1*exp(-theta[1]*theta[2]*theta[4]*theta[6]*(median.Rx12/lambda)^k)*theta[1]*theta[2]*theta[4]*(median.Rx12/lambda)^k-
                   (1-p1)*exp(-theta[1]*theta[3]*theta[5]*theta[6]*(median.Rx12/lambda)^k)*theta[1]*theta[3]*theta[5]*(median.Rx12/lambda)^k)
  H.mat[5,] <- c(p1*theta[2]*theta[6]*exp(-theta[2]*theta[6]*(median.C12/lambda)^k)*k*(median.C12/lambda)^(k-1)*(median.C12/(lambda^2))+
                   (1-p1)*theta[3]*theta[6]*exp(-theta[3]*theta[6]*(median.C12/lambda)^k)*k*(median.C12/lambda)^(k-1)*(median.C12/(lambda^2)), 
                 -p1*theta[2]*theta[6]*exp(-theta[2]*theta[6]*(median.C12/lambda)^k)*(median.C12/lambda)^k*log(median.C12/lambda)-
                   (1-p1)*theta[3]*theta[6]*exp(-theta[3]*theta[6]*(median.C12/lambda)^k)*(median.C12/lambda)^k*log(median.C12/lambda), 
                 0,
                 -p1*exp(-theta[2]*theta[6]*(median.C12/lambda)^k)*theta[6]*(median.C12/lambda)^k,
                 -(1-p1)*exp(-theta[3]*theta[6]*(median.C12/lambda)^k)*theta[6]*(median.C12/lambda)^k,
                 0,
                 0,
                 -p1*exp(-theta[2]*theta[6]*(median.C12/lambda)^k)*theta[2]*(median.C12/lambda)^k-
                   (1-p1)*exp(-theta[3]*theta[6]*(median.C12/lambda)^k)*theta[3]*(median.C12/lambda)^k)
  H.mat[6,] <- c(p2*theta[1]*theta[6]*exp(-theta[1]*theta[6]*(median.Rx10/lambda)^k)*k*(median.Rx10/lambda)^(k-1)*(median.Rx10/(lambda^2))+
                   (1-p2)*theta[1]*theta[2]*theta[4]*theta[6]*exp(-theta[1]*theta[2]*theta[4]*theta[6]*(median.Rx10/lambda)^k)*k*(median.Rx10/lambda)^(k-1)*(median.Rx10/(lambda^2)), 
                 -p2*theta[1]*theta[6]*exp(-theta[1]*theta[6]*(median.Rx10/lambda)^k)*(median.Rx10/lambda)^k*log(median.Rx10/lambda)-
                   (1-p2)*theta[1]*theta[2]*theta[4]*theta[6]*exp(-theta[1]*theta[2]*theta[4]*theta[6]*(median.Rx10/lambda)^k)*(median.Rx10/lambda)^k*log(median.Rx10/lambda), 
                 -p2*exp(-theta[1]*theta[6]*(median.Rx10/lambda)^k)*theta[6]*(median.Rx10/lambda)^k-
                   (1-p2)*exp(-theta[1]*theta[2]*theta[4]*theta[6]*(median.Rx10/lambda)^k)*theta[2]*theta[4]*theta[6]*(median.Rx10/lambda)^k,
                 -(1-p2)*exp(-theta[1]*theta[2]*theta[4]*theta[6]*(median.Rx10/lambda)^k)*theta[1]*theta[4]*theta[6]*(median.Rx10/lambda)^k,
                 0,
                 -(1-p2)*exp(-theta[1]*theta[2]*theta[4]*theta[6]*(median.Rx10/lambda)^k)*theta[1]*theta[2]*theta[6]*(median.Rx10/lambda)^k,
                 0,
                 -p2*exp(-theta[1]*theta[6]*(median.Rx10/lambda)^k)*theta[1]*(median.Rx10/lambda)^k-
                   (1-p2)*exp(-theta[1]*theta[2]*theta[4]*theta[6]*(median.Rx10/lambda)^k)*theta[1]*theta[2]*theta[4]*(median.Rx10/lambda)^k)
  H.mat[7,] <- c(p2*theta[6]*exp(-theta[6]*(median.C10/lambda)^k)*k*(median.C10/lambda)^(k-1)*(median.C10/(lambda^2))+
                   (1-p2)*theta[2]*theta[6]*exp(-theta[2]*theta[6]*(median.C10/lambda)^k)*k*(median.C10/lambda)^(k-1)*(median.C10/(lambda^2)), 
                 -p2*theta[6]*exp(-theta[6]*(median.C10/lambda)^k)*(median.C10/lambda)^k*log(median.C10/lambda)-
                   (1-p2)*theta[2]*theta[6]*exp(-theta[2]*theta[6]*(median.C10/lambda)^k)*(median.C10/lambda)^k*log(median.C10/lambda), 
                 0,
                 -(1-p2)*exp(-theta[2]*theta[6]*(median.C10/lambda)^k)*theta[6]*(median.C10/lambda)^k,
                 0,
                 0,
                 0,
                 -p2*exp(-theta[6]*(median.C10/lambda)^k)*(median.C10/lambda)^k-
                   (1-p2)*exp(-theta[2]*theta[6]*(median.C10/lambda)^k)*theta[2]*(median.C10/lambda)^k)
  
  median.var <- J.mat.inv %*% H.mat %*% param.var %*% t(H.mat) %*% t(J.mat.inv)
  mean.CE4<-c(median.result$r12-median.result$r0,median.result$r2-median.result$r10,median.result$r1-median.result$r0,median.result$r2-median.result$r1)
  mean7<-c(median.result$r0,median.result$r1,median.result$r2,median.Rx12,median.C12,median.Rx10,median.C10)
  mean5<-c(mean7[1:3],log(mean7[4]/mean7[5]),log(mean7[6]/mean7[7]))
  p.mean5<-rep(NA,5)
  mean5.var <- deltamethod(list(~x1,~x2,~x3,~log(x4/x5),~log(x6/x7)),mean=mean7,cov=median.var,ses=FALSE)
  p.mean5[1] = 1-(pmvnorm(lower=rep(-abs(mean5[1]/sqrt(diag(mean5.var)[1])),5),upper=rep(abs(mean5[1]/sqrt(diag(mean5.var)[1])),5),corr=cov2cor(mean5.var))[1])
  p.mean5[2] = 1-(pmvnorm(lower=rep(-abs(mean5[2]/sqrt(diag(mean5.var)[2])),5),upper=rep(abs(mean5[2]/sqrt(diag(mean5.var)[2])),5),corr=cov2cor(mean5.var))[1])
  p.mean5[3] = 1-(pmvnorm(lower=rep(-abs(mean5[3]/sqrt(diag(mean5.var)[3])),5),upper=rep(abs(mean5[3]/sqrt(diag(mean5.var)[3])),5),corr=cov2cor(mean5.var))[1])
  p.mean5[4] = 1-(pmvnorm(lower=rep(-abs(mean5[4]/sqrt(diag(mean5.var)[4])),5),upper=rep(abs(mean5[4]/sqrt(diag(mean5.var)[4])),5),corr=cov2cor(mean5.var))[1])
  p.mean5[5] = 1-(pmvnorm(lower=rep(-abs(mean5[5]/sqrt(diag(mean5.var)[5])),5),upper=rep(abs(mean5[5]/sqrt(diag(mean5.var)[5])),5),corr=cov2cor(mean5.var))[1])
  
  
  CE4.var.full <- deltamethod(list(~x1,~x2,~x3,~(log(x4/x5)-x1),~(x3-log(x6/x7)),~(x2-x1),~(x3-x2)),mean=mean7,cov=median.var,ses=FALSE)
  CE4.var <- CE4.var.full[4:7,4:7]
  q<-qmvnorm(1-alpha,tail=("both.tails"),corr=cov2cor(CE4.var))$quantile
  CI.CE.4<-matrix(rep(NA,8),nrow=4)
  CI.CE.4[,1]<-mean.CE4-q*sqrt(diag(CE4.var))
  CI.CE.4[,2]<-mean.CE4+q*sqrt(diag(CE4.var))
  p.value<-p.value.ind<-rep(NA,4)
  p.value.ind[1]<-1-pnorm(abs(mean.CE4[1])/sqrt(diag(CE4.var)[1]))
  p.value.ind[2]<-1-pnorm(abs(mean.CE4[2])/sqrt(diag(CE4.var)[2]))
  p.value.ind[3]<-1-pnorm(abs(mean.CE4[3])/sqrt(diag(CE4.var)[3]))
  p.value.ind[4]<-1-pnorm(abs(mean.CE4[4])/sqrt(diag(CE4.var)[4]))
  p.value[1] = 1-(pmvnorm(lower=rep(-abs(mean.CE4[1]/sqrt(diag(CE4.var)[1])),4),upper=rep(abs(mean.CE4[1]/sqrt(diag(CE4.var)[1])),4),corr=cov2cor(CE4.var))[1])
  p.value[2] = 1-(pmvnorm(lower=rep(-abs(mean.CE4[2]/sqrt(diag(CE4.var)[2])),4),upper=rep(abs(mean.CE4[2]/sqrt(diag(CE4.var)[2])),4),corr=cov2cor(CE4.var))[1])
  p.value[3] = 1-(pmvnorm(lower=rep(-abs(mean.CE4[3]/sqrt(diag(CE4.var)[3])),4),upper=rep(abs(mean.CE4[3]/sqrt(diag(CE4.var)[3])),4),corr=cov2cor(CE4.var))[1])
  p.value[4] = 1-(pmvnorm(lower=rep(-abs(mean.CE4[4]/sqrt(diag(CE4.var)[4])),4),upper=rep(abs(mean.CE4[4]/sqrt(diag(CE4.var)[4])),4),corr=cov2cor(CE4.var))[1])
  p.maxN = 1-(pmvnorm(lower=rep(-max(abs(mean.CE4/sqrt(diag(CE4.var)))),4),upper=rep(max(abs(mean.CE4/sqrt(diag(CE4.var)))),4),corr=cov2cor(CE4.var))[1])
  return(list(mean5, p.mean5, mean7, median.var = median.var,median5.var=mean5.var,CE4.var = CE4.var,mean.CE4=mean.CE4,CI.CE.4=CI.CE.4,p.value=p.value,p.maxN=p.maxN))
}

