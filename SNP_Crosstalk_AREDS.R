require(geometry)

fulldata_AREDS<-read.csv("Phenotype_AREDS.csv")
topSNPs<-read.csv("topSNP_m10_AREDS.csv")
results40<-topSNPs[which(topSNPs$GENE=="ESRRB,VASH1"|topSNPs$GENE=="VASH1"|topSNPs$GENE=="SPOCK2"|topSNPs$GENE=="C19orf44"|topSNPs$GENE=="C19orf44,CALR3"),]
SNP40<-fulldata_AREDS[,as.character(results40$SNP)]

##coordinates
##other SNP within selsected SNP: rs147106198
target_coord<-matrix(rep(0,9*ncol(SNP40)),ncol=9)
for (i in 1:ncol(SNP40)){
  table1<-table(SNP40[,i],SNP40[,which(colnames(SNP40)=="rs147106198")])
  perc<-prop.table(table1,2)
  target_coord[i,1]<-perc[1,1]
  target_coord[i,2]<-perc[2,1]
  target_coord[i,3]<-perc[3,1]
  target_coord[i,4]<-perc[1,2]
  target_coord[i,5]<-perc[2,2]
  target_coord[i,6]<-perc[3,2]
  target_coord[i,7]<-perc[1,3]
  target_coord[i,8]<-perc[2,3]
  target_coord[i,9]<-perc[3,3]
}
target_coord<-apply(target_coord,2,as.numeric)

#######This step is to create the visualization the coordinates on a ternary plot
y2contrast<-function(y)
{
  n <- nrow(y)
  if(is.matrix(y)) {
    v <- rep(0., n)
    h <- rep(0., n)
    for(i in 1.:n) {
      v[i] <- (2. * y[i, 3.] - (y[i, 1.] + y[i, 2.]))/sqrt(
        6.)
      h[i] <- (y[i, 2.] - y[i, 1.])/sqrt(2.)
    }
  }
  else {
    v <- (2. * y[3.] - (y[1.] + y[2.]))/sqrt(6.)
    h <- (y[2.] - y[1.])/sqrt(2.)
  }
  #	return(list(h=h, v=v))
  return(cbind(h, v))
}

topnum<-40
beta_sim<-matrix(rep(0,9*topnum),ncol=3)
coord<-as.matrix(cbind(colnames(SNP40),target_coord))
for (i in 1:topnum){
  beta_sim[c(3*i-2,3*i-1,3*i),]<-rbind(coord[i,c(2:4)],coord[i,c(5:7)],coord[i,c(8:10)])
}

beta_sim<-apply(beta_sim,2,as.numeric)
psim<-list()
####blue for SPOCK2, red for ESRRB/VASH1, and springgreen4 for C19orf44/CALR3
clist<-c(rep("blue",4),rep("red",30),rep("springgreen4",6))


par(mar=c(3.1,4.1,2.1,2.1))
tri.tips <- y2contrast(rbind(c(1,0,0),c(0,1,0),c(0,0,1)))
trimesh(rbind(1:3), tri.tips, asp=1)
text(c(-0.7571068,0.7571068,0), c(-0.4082483,-0.4082483,0.8464966), c("AA","Aa","aa"),cex=1.5) # Label tips
for (i in c(1:4,6:40)){
  psim[[i]]<-matrix(rep(0,6),nrow=3)
  psim[[i]]<-bary2cart(tri.tips, beta_sim[c(3*i-2,3*i-1,3*i),])
  lines(psim[[i]],col=clist[[i]],lty=1,lwd=2)
}

#####Causal SNP########
causalsnp<-bary2cart(tri.tips,beta_sim[c(13:15),])
lines(causalsnp,col="black",lty=2,lwd=3)


