taj<-read.table("~/analysis/angsd/taj",header=TRUE, sep=',')

taj2<-as.data.frame(taj[,4:10])

pops<-c("bbvb","bbpb","bbsjsp","bbbnp","bbsp","bbgb",
		"vbpb","vbsjsp","vbbnp","vbsp","vbgb",
		"pbsjsp","pbbnp","pbsp","pbgb",
		"sjspbnp","sjspsp","sjspgb",
		"bnpsp","bnpgb",
		"spgb")

k<-0
ddiff<-matrix(nrow=1027369,ncol=21)
colnames(ddiff)<-pops

for (i in 1:6)
{
	ii=i+1
	for (j in ii:7){
		k<-k+1
		ddiff[,k]<-taj2[,i]-taj2[,j]
		}
}


##Run All-pops-tajima.r and d_subtraction.r before this##
pops<-c("bbvb","bbpb","bbsjsp","bbbnp","bbsp","bbgb",
        "vbpb","vbsjsp","vbbnp","vbsp","vbgb",
        "pbsjsp","pbbnp","pbsp","pbgb",
        "sjspbnp","sjspsp","sjspgb",
        "bnpsp","bnpgb",
        "spgb")

dmeans<-colMeans(ddiff,na.rm=TRUE)
dmeans<-as.matrix(dmeans)

library(matrixStats)
dstdev<-colSds(ddiff,na.rm=TRUE)
dstdev<-as.matrix(dstdev)

ztaj<-matrix(nrow=1027369,ncol=21)
colnames(ztaj)<-pops


for(i in 1:21)
{
  for(k in 1:1027369)
  {
    ztaj[k,i]<-(ddiff[k,i]-dmeans[i])/dstdev[i]
  }
}


############################### don't need it since started doing per base
###keeping only the ones with over 20 rep
cov<-matrix(nrow=1027369,ncol=7)
cov<-cbind(bb[,5],vb[,5],pb[,5],sj[,5],bnp[,5],sp[,5],gb[,5])

nsnps<-cov[,1]
for (i in 1:7){
  nsnps <- nsnps + cov[,i]
}
nsnps <- nsnps/7

subw <- nsnps > 20

ztaj_keep<-cbind(taj[,1:3],ztaj)

write.table(ztaj_keep,"~/analysis/angsd/ztaj_keep_5kb",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")

#########
taj<-read.table("~/analysis/angsd/ztaj_keep_5kb",header=FALSE,sep='\t')
pops<-c("scaf","start","end","bbvb","bbpb","bbsjsp","bbbnp","bbsp","bbgb",
        "vbpb","vbsjsp","vbbnp","vbsp","vbgb",
        "pbsjsp","pbbnp","pbsp","pbgb",
        "sjspbnp","sjspsp","sjspgb",
        "bnpsp","bnpgb",
        "spgb")
colnames(taj)<-pops

#####PBS comparison in taj
distsp<-cbind(taj[,1:3],taj[,8],taj[,13],taj[,17],taj[,20],taj[,22])
distgb<-cbind(taj[,1:3],taj[,9],taj[,14],taj[,18],taj[,21],taj[,23])

distnames<-c("Scaf","start","end","bb","vb","pb","sj","bnp")
colnames(distsp)<-distnames
colnames(distgb)<-distnames

colnam<-names(distsp)[4:8]

total_dist<-cbind(taj[1:3],distsp[colnam]+distgb[match(distsp$Scaf,distgb$Scaf),colnam])

pbs_dist<-matrix(nrow=1027369,ncol=5)
for (i in 1:5)
{
  pbs_dist[,i]<-(total_dist[,i+3]-taj[,24])/2
}

pbs_dist<-cbind(taj[,1:3],pbs_dist)
colnames(pbs_dist)<-distnames


col<-c()
for (i in 1:5){
  col[i]<-quantile(pbs_dist[,i+3],prob=.99,na.rm=TRUE)
}

collow<-c()
for (i in 1:5){
  collow[i]<-quantile(pbs_dist[,i+3],prob=.01,na.rm=TRUE)
}


###DKGB
pbs_dist[grep("Scaffold1171",pbs_dist$Scaf),]

plot(pbs_dist[grep("Scaffold1171",pbs_dist$Scaf),4],pch=20,cex=.8,col="black",ylim=c(-10,10))
points(pbs_dist[grep("Scaffold1171",pbs_dist$Scaf),5],pch=20,cex=.8,col="grey")
points(pbs_dist[grep("Scaffold1171",pbs_dist$Scaf),6],pch=20,cex=.8,col="red")
points(pbs_dist[grep("Scaffold1171",pbs_dist$Scaf),7],pch=20,cex=.8,col="darkorange")
points(pbs_dist[grep("Scaffold1171",pbs_dist$Scaf),8],pch=20,cex=.8,col="gold")

abline(h=c(col),col=c("black","grey","red","darkorange","gold"),lwd=3,lty=3)
abline(h=c(collow),col=c("black","grey","red","darkorange","gold"),lwd=3,lty=3)

###Aquaporin
Scaffold24<-pbs_dist[grep("Scaffold24",pbs_dist$Scaf),]
plot(Scaffold24[2500:3000,4],pch=20,cex=.8,col="black",ylim=c(-10,10))
points(Scaffold24[2500:3000,5],pch=20,cex=.8,col="grey")
points(Scaffold24[2500:3000,6],pch=20,cex=.8,col="red")
points(Scaffold24[2500:3000,7],pch=20,cex=.8,col="darkorange")
points(Scaffold24[2500:3000,8],pch=20,cex=.8,col="gold")

abline(h=c(col),col=c("black","grey","red","darkorange","gold"),lwd=3,lty=3)
abline(h=c(collow),col=c("black","grey","red","darkorange","gold"),lwd=3,lty=3)

###ARNT
Scaffold0<-pbs_dist[grep("Scaffold0",pbs_dist$Scaf),]
plot(Scaffold0[600:900,4],pch=20,cex=.8,col="black",ylim=c(-10,10))
points(Scaffold0[600:900,5],pch=20,cex=.8,col="grey")
points(Scaffold0[600:900,6],pch=20,cex=.8,col="red")
points(Scaffold0[600:900,7],pch=20,cex=.8,col="darkorange")
points(Scaffold0[600:900,8],pch=20,cex=.8,col="gold")

###AHR1/2
pbs_dist[grep("Scaffold900",pbs_dist$Scaf),]

plot(pbs_dist[grep("Scaffold900",pbs_dist$Scaf),4],pch=20,cex=.8,col="black",ylim=c(-10,10))
points(pbs_dist[grep("Scaffold900",pbs_dist$Scaf),5],pch=20,cex=.8,col="grey")
points(pbs_dist[grep("Scaffold900",pbs_dist$Scaf),6],pch=20,cex=.8,col="red")
points(pbs_dist[grep("Scaffold900",pbs_dist$Scaf),7],pch=20,cex=.8,col="darkorange")
points(pbs_dist[grep("Scaffold900",pbs_dist$Scaf),8],pch=20,cex=.8,col="gold")

abline(h=c(col),col=c("black","grey","red","darkorange","gold"),lwd=3,lty=3)
abline(h=c(collow),col=c("black","grey","red","darkorange","gold"),lwd=3,lty=3)


###AHR1a

pbs_dist[grep("Scaffold1482",pbs_dist$Scaf),]

plot(pbs_dist[grep("Scaffold1482",pbs_dist$Scaf),4],pch=20,cex=1.2,col="black",ylim=c(-10,10),
     ylab="Population branch statistic compared to references")
points(pbs_dist[grep("Scaffold1482",pbs_dist$Scaf),5],pch=20,cex=1.2,col="grey")
points(pbs_dist[grep("Scaffold1482",pbs_dist$Scaf),6],pch=20,cex=1.2,col="red")
points(pbs_dist[grep("Scaffold1482",pbs_dist$Scaf),7],pch=20,cex=1.2,col="darkorange")
points(pbs_dist[grep("Scaffold1482",pbs_dist$Scaf),8],pch=20,cex=1.2,col="gold")

abline(h=c(col),col=c("black","grey","red","darkorange","gold"),lwd=3,lty=3)
abline(h=c(collow),col=c("black","grey","red","darkorange","gold"),lwd=3,lty=3)


######################

pdf(file="~/share/VCF/angsd_res/tajd_res.pdf")

par(mfrow=c(1,1))
plot(ddiff[,"bbgb"],pch=20,cex=.1)
quantile(ddiff[,"bbgb"],probs=.01,na.rm=TRUE)
abline(0.0266,0,col="red")
plot(ddiff[,"vbgb"],pch=20,cex=.1)
quantile(ddiff[,"vbgb"],probs=.01,na.rm=TRUE)
abline(0.035,0,col="red")
plot(ddiff[,"pbgb"],pch=20,cex=.1)
quantile(ddiff[,"pbgb"],probs=.01,na.rm=TRUE)
abline(0.020,0,col="red")




dev.off()
