require(rootSolve) # need this package for the median.surv function
require(survival)
require(eha)
require(Hmisc)
require(rmeta)

source("CE4_Weibull.R")


fulldata_AREDS<-read.csv("Phenotype_AREDS.csv")

CFH_ARMS2<-c("rs7522681","rs412852","rs1061170","rs3766405","rs11200647","rs10490924")



alpha<-10/3837556

est.CE<-matrix(rep(NA,4*6),nrow=6)
p.value<-matrix(rep(NA,4*6),nrow=6)
CI.4<-matrix(rep(NA,8*6),nrow=6)
p.maxN<-matrix(rep(NA,6),nrow=6)

for(i in 1:6){
  data_i<-fulldata_AREDS[,c(which(names(fulldata_AREDS)==CFH_ARMS2[i]),2:9)]
  colnames(data_i)[c(1,9)]<-c("M","Delta")
  data_i$Y<-as.numeric(as.character(data_i$Y))
  data_i$enrollage<-as.numeric(as.character(data_i$enrollage))
  data_i$SevScaleBL<-as.numeric(as.character(data_i$SevScaleBL))
  data_i$Delta<-as.numeric(as.character(data_i$Delta))
  data_i$smoke<-as.numeric(as.character(data_i$smoke))
  age.mean<-mean(data_i$enrollage)
  SevScaleBL.mean<-mean(data_i$SevScaleBL)
  smoke=0
  weibfit<-weib(data_i,varlist=1,age.mean=age.mean,SevScaleBL.mean=SevScaleBL.mean,smoke=smoke)
  med.surv<-median.surv(coef=weibfit$coef,p1=weibfit$p1,p2=weibfit$p2,varlist=1,tau=0.75)
  CEresult<-var.median.surv.dt.dim7(median.result=med.surv,p1=weibfit$p1,p2=weibfit$p2,coef=weibfit$coef,param.var=weibfit$var,varlist=1,alpha=alpha)
  est.CE[i,]<-exp(CEresult$mean.CE4)
  CI.4[i,c(1,2)]<-exp(CEresult$CI.CE.4[1,])
  CI.4[i,c(3,4)]<-exp(CEresult$CI.CE.4[2,])
  CI.4[i,c(5,6)]<-exp(CEresult$CI.CE.4[3,])
  CI.4[i,c(7,8)]<-exp(CEresult$CI.CE.4[4,])
  p.value[i,]<-CEresult$p.value
  p.maxN[i,]<-CEresult$p.maxN
}

result<-cbind(CFH_ARMS2,p.maxN,est.CE[,1],CI.4[,1:2],p.value[,1],est.CE[,2],CI.4[,3:4],p.value[,2],est.CE[,3],CI.4[,5:6],p.value[,3],est.CE[,4],CI.4[,7:8],p.value[,4])
result<-as.data.frame(result)
names(result)<-c("SNP","p.maxN","theta_12,0","lower_12,0","higher_12,0","p_12,0","theta_2,01","lower_2,01","higher_2,01","p_2,01",
                 "theta_1,0","lower_1,0","higher_1,0","p_1,0","theta_2,1","lower_2,1","higher_2,1","p_2,1")
View(result)
