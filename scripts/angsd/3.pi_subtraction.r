pi<-read.table("~/analysis/data/angsd/pi_neut_5kb",header=TRUE,sep=',')
pi2<-as.data.frame(pi[,4:10])

pops<-c("scaf","start","end","bbvb","bbpb","bbsjsp","bbbnp","bbsp","bbgb",
        "vbpb","vbsjsp","vbbnp","vbsp","vbgb",
        "pbsjsp","pbbnp","pbsp","pbgb",
        "sjspbnp","sjspsp","sjspgb",
        "bnpsp","bnpgb",
        "spgb","keep")

pidiff_temp<-matrix(nrow=1026857,ncol=21)

k<-0
for (i in 1:6)
{
  ii=i+1
  for (j in ii:7){
    k<-k+1
    pidiff_temp[,k]<-pi2[,i]-pi2[,j]
  }
}
#####################run until here for pi graphs


pimeans<-colMeans(pidiff_temp,na.rm=TRUE)
pimeans<-as.matrix(pimeans)

library(matrixStats)
pistdev<-colSds(pidiff_temp,na.rm=TRUE)
pistdev<-as.matrix(pistdev)

zpi_temp<-matrix(nrow=1026857,ncol=21)


for(i in 1:21)
{
  for(k in 1:1026857)
  {
    zpi_temp[k,i]<-(pidiff_temp[k,i]-pimeans[i])/pistdev[i]
  }
}

zpi<-cbind(pi[,1:3],zpi_temp,pi[,11])
colnames(zpi)<-pops

write.table(zpi,"~/analysis/data/angsd/zpi_keep_5kb",col.names=FALSE,row.names=FALSE,quote=FALSE,sep='\t')

###Run from here for pi graphs PBS
pops2<-c("scaf","start","end","bbvb","bbpb","bbsjsp","bbbnp","bbsp","bbgb",
        "vbpb","vbsjsp","vbbnp","vbsp","vbgb",
        "pbsjsp","pbbnp","pbsp","pbgb",
        "sjspbnp","sjspsp","sjspgb",
        "bnpsp","bnpgb",
        "spgb")

pidiff<-cbind(pi[,1:3],pidiff_temp)
colnames(pidiff)<-pops2

distsp<-cbind(pidiff[,1:3],pidiff[,8],pidiff[,13],pidiff[,17],pidiff[,20],pidiff[,22])
distgb<-cbind(pidiff[,1:3],pidiff[,9],pidiff[,14],pidiff[,18],pidiff[,21],pidiff[,23])

distnames<-c("Scaf","start","end","bb","vb","pb","sj","bnp")
colnames(distsp)<-distnames
colnames(distgb)<-distnames

colnam<-names(distsp)[4:8]

distsp2<-cbind(seq=seq(1:1026857),distsp)
distgb2<-cbind(seq=seq(1:1026857),distgb)

total_dist<-cbind(pidiff[1:3],distsp2[colnam]+distgb2[match(distsp2$seq,distgb2$seq),colnam])

pbs_dist<-matrix(nrow=1026857,ncol=5)
for (i in 1:5)
{
  pbs_dist[,i]<-(total_dist[,i+3]-pidiff[,24])/2
}

pbs_dist<-cbind(pidiff[,1:3],pbs_dist)
colnames(pbs_dist)<-distnames


col<-c()
for (i in 1:5){
  col[i]<-quantile(pbs_dist[,i+3],prob=.99,na.rm=TRUE)
}

collow<-c()
for (i in 1:5){
  collow[i]<-quantile(pbs_dist[,i+3],prob=.01,na.rm=TRUE)
}

###chr1
pbs_dist[grep("chr1",pbs_dist$Scaf),]

plot(pbs_dist[grep("chr1\\b",pbs_dist$Scaf),4],pch=20,cex=.8,col="black",ylim=c(-.1,.1),
     ylab="delta pi (lower means resistant populations have lower diversity)",
     xlab="chr1 position")
points(pbs_dist[grep("chr1\\b",pbs_dist$Scaf),5],pch=20,cex=.8,col="grey")
points(pbs_dist[grep("chr1\\b",pbs_dist$Scaf),6],pch=20,cex=.8,col="red")
points(pbs_dist[grep("chr1\\b",pbs_dist$Scaf),7],pch=20,cex=.8,col="darkorange")
points(pbs_dist[grep("chr1\\b",pbs_dist$Scaf),8],pch=20,cex=.8,col="gold")

abline(h=c(collow),col=c("black","grey","red","darkorange","gold"),lwd=3,lty=3)

legend(x=c(300000,350000),y=c(.05,.1),c("BB","VB","PB","SJSP","BNP"),pch=20,cex=1.2,
       col=c("black","grey","red","darkorange","gold"), x.intersp=.3,y.intersp=.6)

#chr1 pi
plot(pi[grep("chr1\\b",pi$scaf),4],pch=20,cex=.5,col="black",ylim=c(-.1,.4))
points(pi[grep("chr1\\b",pi$scaf),5],pch=20,cex=.5,col="grey",ylim=c(-.1,.4))
points(pi[grep("chr1\\b",pi$scaf),6],pch=20,cex=.5,col="red",ylim=c(-.1,.4))
points(pi[grep("chr1\\b",pi$scaf),7],pch=20,cex=.5,col="darkorange",ylim=c(-.1,.4))
points(pi[grep("chr1\\b",pi$scaf),8],pch=20,cex=.5,col="gold",ylim=c(-.1,.4))
points(pi[grep("chr1\\b",pi$scaf),9],pch=20,cex=.5,col="cyan",ylim=c(-.1,.4))
points(pi[grep("chr1\\b",pi$scaf),10],pch=20,cex=.5,col="blue",ylim=c(-.1,.4))
