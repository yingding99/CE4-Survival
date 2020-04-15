require(survival)
require(ggplot2)
require(survminer)
require(gridExtra)


fulldata_AREDS<-read.csv("Phenotype_AREDS.csv")
AREDS2_ctrl<-read.csv("Phenotype_AREDS2_trtAREDS.csv")

#####################K-M curves###############
####targeted population is {Aa,aa}, which means the rs147106198 coding is 1 or 2
snp1_Trt<-fulldata_AREDS[which(fulldata_AREDS$Trt==1),]
snp1_Trt$target<-ifelse(snp1_Trt$rs147106198==0,0,1)
fit1_t<-survfit(Surv(snp1_Trt$Y,as.numeric(as.character(snp1_Trt$status))) ~ snp1_Trt$target)

snp2_Trt<-AREDS2_ctrl
snp2_Trt$target<-ifelse(snp2_Trt$rs147106198==0,0,1)
fit2_t<-survfit(Surv(snp2_Trt$Y,as.numeric(as.character(snp2_Trt$status))) ~ snp2_Trt$target)

fit<-list(fit1_t,fit2_t)
data<-list(snp1_Trt,snp2_Trt)

labss<-c("Non-target","Target")

for (i in 1:2){
  assign(paste("p",i,sep=""),ggsurvplot(
    fit[[i]], 
    data = data[[i]], 
    censor = FALSE,
    #censor.shape=124,
    #censor.size=2,
    size = 1,                 # change line size
    #palette = cols1,# custom color palettes
    conf.int = F,          # Add confidence interval
    pval = FALSE,              # Add p-value
    risk.table = FALSE,        # Add risk table
    #legend.labs=labss, 
    legend="none",  # Change legend labels
    linetype=c(1,2),
    ylim=c(0,1),
    xlab="",
    ylab="",
    ggtheme = theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                    panel.background = element_blank(), axis.line = element_line(colour = "black"),
                    axis.text=element_text(size=22),
                    axis.title=element_text(size=22),
                    legend.title = element_text(size=20),
                    legend.text = element_text(size=20),plot.margin = unit(c(1,0.3,0.3,0.3), "cm"))
  )
  )
}
pall<-grid.arrange(p1$plot, p2$plot, nrow = 1, ncol=2)
plot(pall)
