require(tableone)

fulldata_AREDS<-read.csv("Phenotype_AREDS.csv")
fulldata_AREDS$smoke<-as.factor(fulldata_AREDS$smoke)
fulldata_AREDS$status<-as.factor(fulldata_AREDS$status)
#####################################################################
##########                       Table 1                 ############
#####################################################################
table1_all<-CreateTableOne(data = fulldata_AREDS,vars=c("enrollage","smoke","Sex","SevScaleBL","status"))
table1_all$CatTable
table1_all$ContTable
summary(fulldata_AREDS$enrollage)
summary(fulldata_AREDS$SevScaleBL)
summary(fulldata_AREDS$Y)
table1_all<-CreateTableOne(data = fulldata_AREDS,vars=c("enrollage","smoke","Sex","SevScaleBL","status"),strata="Trt")
table1_all$CatTable
table1_all$ContTable
summary(fulldata_AREDS[which(fulldata_AREDS$Trt==0),]$enrollage)
summary(fulldata_AREDS[which(fulldata_AREDS$Trt==1),]$enrollage)
summary(fulldata_AREDS[which(fulldata_AREDS$Trt==0),]$SevScaleBL)
summary(fulldata_AREDS[which(fulldata_AREDS$Trt==1),]$SevScaleBL)


#####################################################################
##########                       Table 2                 ############
#####################################################################
fulldata_AREDS$target<-ifelse(fulldata_AREDS$rs147106198==0,0,1)
fulldata_AREDS$Trt<-as.factor(fulldata_AREDS$Trt)
table1_all<-CreateTableOne(data = fulldata_AREDS,vars=c("enrollage","smoke","Sex","Trt","SevScaleBL"),strata="target")
table1_all$CatTable
table1_all$ContTable
summary(fulldata_AREDS[which(fulldata_AREDS$target==0),]$enrollage)
summary(fulldata_AREDS[which(fulldata_AREDS$target==1),]$enrollage)
summary(fulldata_AREDS[which(fulldata_AREDS$target==0),]$SevScaleBL)
summary(fulldata_AREDS[which(fulldata_AREDS$target==1),]$SevScaleBL)


#####################################################################
##########                       Table 4                 ############
#####################################################################


#################AREDS: AREDS formulation arm
AREDS<-fulldata_AREDS[which(fulldata_AREDS$Trt==1),]

table1_all<-CreateTableOne(data = AREDS,vars=c("enrollage","smoke","Sex","SevScaleBL"),strata="target")
table1_all$CatTable
table1_all$ContTable
summary(AREDS[which(AREDS$target==0),]$enrollage)
summary(AREDS[which(AREDS$target==1),]$enrollage)
summary(AREDS[which(AREDS$target==0),]$SevScaleBL)
summary(AREDS[which(AREDS$target==1),]$SevScaleBL)


#################AREDS2: AREDS formulation arm
AREDS2_ctrl<-read.csv("Phenotype_AREDS2_trtAREDS.csv")
AREDS2_ctrl$target<-ifelse(AREDS2_ctrl$rs147106198==0,0,1)
AREDS2_ctrl$smoke<-as.factor(AREDS2_ctrl$smoke)

table1_all<-CreateTableOne(data = AREDS2_ctrl,vars=c("enrollage","smoke","Sex","SevScaleBL"),strata="target")
table1_all$CatTable
table1_all$ContTable
summary(AREDS2_ctrl[which(AREDS2_ctrl$target==0),]$enrollage)
summary(AREDS2_ctrl[which(AREDS2_ctrl$target==1),]$enrollage)
summary(AREDS2_ctrl[which(AREDS2_ctrl$target==0),]$SevScaleBL)
summary(AREDS2_ctrl[which(AREDS2_ctrl$target==1),]$SevScaleBL)
