require(rootSolve) 
require(survival)
require(eha)
require(Hmisc)
require(rmeta)
require(plotrix)
source("CE4_Weibull.R")

fulldata_AREDS<-read.csv("Phenotype_AREDS.csv")

data_i<-fulldata_AREDS[,c("rs147106198","enrollage","SevScaleBL","smoke","Y","Trt","status")]
names(data_i)[c(1,7)]<-c("M","Delta")
data_i$Y<-as.numeric(as.character(data_i$Y))
data_i$enrollage<-as.numeric(as.character(data_i$enrollage))
data_i$SevScaleBL<-as.numeric(as.character(data_i$SevScaleBL))
data_i$Delta<-as.numeric(as.character(data_i$Delta))
data_i$smoke<-as.numeric(as.character(data_i$smoke))

age.mean<-mean(data_i$enrollage)
SevScaleBL.mean<-mean(data_i$SevScaleBL)
smoke=0

alpha<-10/3837556
weibfit<-weib(data_i,varlist=1,age.mean=age.mean,SevScaleBL.mean=SevScaleBL.mean,smoke=smoke)
med.surv<-median.surv(coef=weibfit$coef,p1=weibfit$p1,p2=weibfit$p2,varlist=1,tau=0.75)
CEresult<-var.median.surv.dt.dim7(median.result=med.surv,p1=weibfit$p1,p2=weibfit$p2,coef=weibfit$coef,param.var=weibfit$var,varlist=1,alpha=alpha)

mu.marker <- rep(0,3);
se.marker <- rep(0,3);

mu.marker[1] <- med.surv$r0
mu.marker[2] <- med.surv$r1
mu.marker[3] <- med.surv$r2
se.marker[1] <- sqrt(CEresult$median.var[1,1])
se.marker[2] <- sqrt(CEresult$median.var[2,2])
se.marker[3] <- sqrt(CEresult$median.var[3,3])

X11(display = "", width=1200, height=600)
par(mfrow=c(1,3))

######plot 1: profile in each biomarker group
errbar(x=c(0,1,2),y=mu.marker,
       yplus=mu.marker+se.marker,yminus=mu.marker-se.marker,
       xlab="",ylab="",
       xaxt = 'n',xlim=c(0,2),ylim=c(min(mu.marker-se.marker)-0.3,max(mu.marker+se.marker)+0.3),cex.axis=2.3);
axis(side=1,at=c(0,1,2),labels=c("AA","Aa","aa"),cex.axis=2.3);
lines(x=c(0,1,2),y=mu.marker, lwd=2);
abline(h=0,lty=2,col="blue",lwd=2)

#plot 2: CE4
CI.CE4.mat <- matrix(CEresult$CI.CE.4,4,2,byrow=FALSE)
mean.CE4 <- CEresult$mean.CE4
errbar(x=1:4,y=mean.CE4,
       yplus=CI.CE4.mat[,2],yminus=CI.CE4.mat[,1],
       xlab="",ylab="",
       xaxt = 'n',xlim=c(0.8,4.2), cex.axis=2.3,
       ylim=c(min(CI.CE4.mat[,1])-0.3,max(CI.CE4.mat[,2])+0.3),
       cex.lab=2.3)
labels.CE4 <- c(expression(kappa["(1,2):0"]),expression(kappa["2:(0,1)"]),
                expression(kappa["1:0"]),expression(kappa["2:1"]))
axis(side=1,at=1:4,labels=labels.CE4,font.axis=1.5,cex.axis=2.6,line=0.5,tick=FALSE)
abline(h=0,lty=2,col="blue",lwd=2)

#####plot 3: construct CI for targeted and non-targeted group, simultaneous 95% CP
alpha=0.05
CEresult<-var.median.surv.dt.dim7(median.result=med.surv,p1=weibfit$p1,p2=weibfit$p2,coef=weibfit$coef,param.var=weibfit$var,varlist=1,alpha=alpha)
newestimate<-c(med.surv$r12,med.surv$r0)
newvariance<-rbind(c(CEresult$median5.var[4,4],CEresult$median5.var[1,4]),c(CEresult$median5.var[4,1],CEresult$median5.var[1,1]))
q<-qmvnorm(1-alpha,tail=("both.tails"),corr=cov2cor(newvariance))$quantile
targetCI<-c(newestimate[1]-q*sqrt(diag(newvariance)[1]),newestimate[1]+q*sqrt(diag(newvariance)[1]))
nontargetCI<-c(newestimate[2]-q*sqrt(diag(newvariance)[2]),newestimate[2]+q*sqrt(diag(newvariance)[2]))

UB<-c(targetCI[2],nontargetCI[2])
LB<-c(targetCI[1],nontargetCI[1])

plotCI(1:2, newestimate, ui=UB, li=LB, xaxt = "n", xlab = "", ylab = "",
       xlim = c(0.5, 2.5), ylim = c(-1, 1), pch = 16, cex.axis = 1.5)
abline(h = 0)
axis(1, at=1:2, cex.axis = 1.5, labels=c("Targeted: {Aa,aa}", "Non-targeted: {AA}"), tick=FALSE)


