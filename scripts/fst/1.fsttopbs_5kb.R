fs <- list.files("~/analysis/fst/raw/", "*5kb1kb.bed",full.names=TRUE)

fst <- list()

for (i in 1:21){
	fst[[i]] <- read.table(fs[i],stringsAsFactors=FALSE)
	fst[[i]][,4] <- as.numeric(fst[[i]][,4])
}
print("done reading files")
nfs <- gsub(".*\\/","",fs)
nfs <- gsub(".fst.*","",nfs)
names(fst)<-nfs

pbs <- function(t1,t2,c12){
  
  t1 <- -log(1-t1)
  t2 <- -log(1-t2)
  c12 <- -log(1-c12)
  
  stat <- (t1 + t2 - c12)/2
  return(stat)
}

nsnps <-fst[[1]][,5]

for (i in 2:21){
  
  nsnps <- nsnps + fst[[i]][,5]
}

nsnps <- nsnps/21

subw <- nsnps > 20

BBpbs <- pbs(fst[["BB.GB"]][,4],fst[["BB.SP"]][,4],fst[["GB.SP"]][,4])

VBpbs <- pbs(fst[["VB.GB"]][,4],fst[["VB.SP"]][,4],fst[["GB.SP"]][,4])

PBpbs <- pbs(fst[["PB.GB"]][,4],fst[["PB.SP"]][,4],fst[["GB.SP"]][,4])

SJpbs <- pbs(fst[["SJ.GB"]][,4],fst[["SJ.SP"]][,4],fst[["GB.SP"]][,4])

BNPpbs <- pbs(fst[["BNP.GB"]][,4],fst[["BNP.SP"]][,4],fst[["GB.SP"]][,4])
print("done reading files2")

Allpbs <- cbind(
  fst[[1]][,1:3],
  BB = BBpbs,
  VB = VBpbs,
  PB = PBpbs,
  SJ = SJpbs,
  BNP = BNPpbs,
  keep = as.numeric(subw))

write.table(Allpbs,
            file="~/analysis/fst/allpbs5kb",
            sep = "\t", 
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)


plot(BBpbs[subw],pch=20,cex=.5, col=factor(Allpbs[,1]))
quantile(BBpbs[subw],probs=.99,na.rm=TRUE)
abline(h=.1804205)


####Gene comparisons

pbs<-read.table("~/analysis/data/fst/allpbs5kb",header=FALSE)
pbsname<-c("Scaf","start","end","BBpbs","VBpbs","PBpbs","SJpbs","BNPpbs","keep")
colnames(pbs)<-pbsname

col<-c()
for (i in 1:5){
  col[i]<-quantile(pbs[,i+3],prob=.99,na.rm=TRUE)
}


###Outlier regions

chr2<-pbs[grep("chr2\\b",pbs$Scaf),]
plot(chr2[27000:28000,4],pch=20,cex=.5,col="black",ylim=c(-.05,.3),
     ylab="chr2_distance",
     xlab="Chr2:27000000:28000000 (in 5kb chunks)")
points(chr2[27000:28000,5],pch=20,cex=.5,col="grey")
points(chr2[27000:28000,6],pch=20,cex=.5,col="red")
points(chr2[27000:28000,7],pch=20,cex=.5,col="darkorange")
points(chr2[27000:28000,8],pch=20,cex=.5,col="gold")
legend("topright",y=c(.05,.1),c("BB","VB","PB","SJSP","BNP"),pch=20,cex=1.2,
       col=c("black","grey","red","darkorange","gold"), x.intersp=.3,y.intersp=.6)
abline(h=0,col="purple",lty=2)


################################################AHR2a

###finding neutral regions to do admixture on

pbs<-read.table("~/analysis/data/fst/allpbs5kb",header=FALSE)
pbsname<-c("Scaf","start","end","BBpbs","VBpbs","PBpbs","SJpbs","BNPpbs","keep")
colnames(pbs)<-pbsname

npbs<-c()
for (i in 1:1026857){
  npbs[i]<-rowSums((pbs[i,4:8])/5)
}
subw2<-npbs<0.02
subw3<-subw+subw2
subw4<-subw3>1

plot(pbsalt[subw4,4])

pbsalt<- cbind(
  fst[[1]][,1:3],
  BB = BBpbs,
  VB = VBpbs,
  PB = PBpbs,
  SJ = SJpbs,
  BNP = BNPpbs,
  keep=as.numeric(subw4))


write.table(pbsalt,"~/share/fst/pbsneut",
            quote=FALSE,row.names = FALSE,col.names = FALSE,sep='\t')
table(pbsalt[,9])



###OLD OUTLIER GRAPHS

# scaf1195<-Allpbs[grep("Scaffold1195",Allpbs$Scaf),]
# 
# plot(scaf1195[1:51,4],pch=20,cex=1.2,col="black",ylim=c(-.2,2.5),xlim=c(0,500))
# points(scaf1195[1:51,5],pch=20,cex=1.2,col="grey")
# points(scaf1195[1:51,6],pch=20,cex=1.2,col="red")
# points(scaf1195[1:51,7],pch=20,cex=1.2,col="darkorange")
# points(scaf1195[1:51,8],pch=20,cex=1.2,col="gold")
# 
# abline(h=c(col),col=c("black","grey","red","darkorange","gold"),lwd=3,lty=3)
# 
# scaf8639<-Allpbs[grep("Scaffold8639",Allpbs$Scaf),]
# 
# plot(scaf8639[,4],pch=20,cex=1.2,col="black",ylim=c(-.2,2.5),xlim=c(0,500))
# points(scaf8639[,5],pch=20,cex=1.2,col="grey")
# points(scaf8639[,6],pch=20,cex=1.2,col="red")
# points(scaf8639[,7],pch=20,cex=1.2,col="darkorange")
# points(scaf8639[,8],pch=20,cex=1.2,col="gold")
# 
# abline(h=c(col),col=c("black","grey","red","darkorange","gold"),lwd=3,lty=3)
# 
# scaf1650<-Allpbs[grep("Scaffold1650",Allpbs$Scaf),]
# 
# plot(scaf1650[100:114,4],pch=20,cex=1.2,col="black",ylim=c(-.2,2.5),xlim=c(0,500))
# points(scaf1650[100:114,5],pch=20,cex=1.2,col="grey")
# points(scaf1650[100:114,6],pch=20,cex=1.2,col="red")
# points(scaf1650[100:114,7],pch=20,cex=1.2,col="darkorange")
# points(scaf1650[100:114,8],pch=20,cex=1.2,col="gold")
# 
# abline(h=c(col),col=c("black","grey","red","darkorange","gold"),lwd=3,lty=3)
# 
# ###Aquaporin
# Scaffold24<-Allpbs[grep("Scaffold24",Allpbs$Scaf),]
# plot(Scaffold24[2500:3000,4],pch=20,cex=.8,col="black",ylim=c(-.5,2.5))
# points(Scaffold24[2500:3000,5],pch=20,cex=.8,col="grey")
# points(Scaffold24[2500:3000,6],pch=20,cex=.8,col="red")
# points(Scaffold24[2500:3000,7],pch=20,cex=.8,col="darkorange")
# points(Scaffold24[2500:3000,8],pch=20,cex=.8,col="gold")
# 
# abline(h=c(col),col=c("black","grey","red","darkorange","gold"),lwd=3,lty=3)
# 
# abline(v=c(200,206))
# abline(v=c(257,272))
# 
# 
# ###ARNT
# Scaffold0<-Allpbs[grep("Scaffold0",Allpbs$Scaf),]
# plot(Scaffold0[600:1000,4],pch=20,cex=.8,col="black",ylim=c(-.5,2.5))
# points(Scaffold0[600:1000,5],pch=20,cex=.8,col="grey")
# points(Scaffold0[600:1000,6],pch=20,cex=.8,col="red")
# points(Scaffold0[600:1000,7],pch=20,cex=.8,col="darkorange")
# points(Scaffold0[600:1000,8],pch=20,cex=.8,col="gold")
# 
# abline(h=c(col),col=c("black","grey","red","darkorange","gold"),lwd=3,lty=3)
# 
# ###ARNT-full
# Scaffold1305<-Allpbs[grep("Scaffold1305",Allpbs$Scaf),]
# plot(Scaffold1305[,4],pch=20,cex=.8,col="black",ylim=c(-.5,2.5))
# points(Scaffold1305[,5],pch=20,cex=.8,col="grey")
# points(Scaffold1305[,6],pch=20,cex=.8,col="red")
# points(Scaffold1305[,7],pch=20,cex=.8,col="darkorange")
# points(Scaffold1305[,8],pch=20,cex=.8,col="gold")
# 
# abline(h=c(col),col=c("black","grey","red","darkorange","gold"),lwd=3,lty=3)
# 
# ###AHR1/2
# scaf900<-Allpbs[grep("Scaffold900",Allpbs$Scaf),]
# 
# plot(Allpbs[grep("Scaffold900",Allpbs$Scaf),4],pch=20,cex=1.2,col="black",ylim=c(-.2,2.5))
# points(Allpbs[grep("Scaffold900",Allpbs$Scaf),5],pch=20,cex=1.2,col="grey")
# points(Allpbs[grep("Scaffold900",Allpbs$Scaf),6],pch=20,cex=1.2,col="red")
# points(Allpbs[grep("Scaffold900",Allpbs$Scaf),7],pch=20,cex=1.2,col="darkorange")
# points(Allpbs[grep("Scaffold900",Allpbs$Scaf),8],pch=20,cex=1.2,col="gold")
# 
# abline(h=c(col),col=c("black","grey","red","darkorange","gold"),lwd=3,lty=3)
# abline(v=c(212,218))
# abline(v=c(265,345))
# 
# ####GCR
# 
# Allpbs[grep("Scaffold2181",Allpbs$Scaf),]
# 
# plot(Allpbs[grep("Scaffold2181",Allpbs$Scaf),4],pch=20,cex=1.2,col="black",ylim=c(-.2,2.5))
# points(Allpbs[grep("Scaffold2181",Allpbs$Scaf),5],pch=20,cex=1.2,col="grey")
# points(Allpbs[grep("Scaffold2181",Allpbs$Scaf),6],pch=20,cex=1.2,col="red")
# points(Allpbs[grep("Scaffold2181",Allpbs$Scaf),7],pch=20,cex=1.2,col="darkorange")
# points(Allpbs[grep("Scaffold2181",Allpbs$Scaf),8],pch=20,cex=1.2,col="gold")
# 
# abline(h=c(col),col=c("black","grey","red","darkorange","gold"),lwd=3,lty=3)
# 
# 
# 
# ###AHR1a
# 
# scaf1650<-Allpbs[grep("Scaffold1650",Allpbs$Scaf),]
# 
# plot(scaf1650[60:75,4],pch=20,cex=1.2,col="black",ylim=c(-.2,2.5),xlim=c(0,500),
#      ylab="Population branch statistic compared to references")
# points(scaf1650[60:75,5],pch=20,cex=1.2,col="grey")
# points(scaf1650[60:75,6],pch=20,cex=1.2,col="red")
# points(scaf1650[60:75,7],pch=20,cex=1.2,col="darkorange")
# points(scaf1650[60:75,8],pch=20,cex=1.2,col="gold")
# 
# abline(h=c(col),col=c("black","grey","red","darkorange","gold"),lwd=3,lty=3)
# 
# 
# scaf1482<-Allpbs[grep("Scaffold1482",Allpbs$Scaf),]
# 
# plot(scaf1482[99:147,4],pch=20,cex=1.2,col="black",ylim=c(-.2,2.5),xlim=c(0,500),
#      ylab="Population branch statistic compared to references")
# points(scaf1482[99:147,5],pch=20,cex=1.2,col="grey")
# points(scaf1482[99:147,6],pch=20,cex=1.2,col="red")
# points(scaf1482[99:147,7],pch=20,cex=1.2,col="darkorange")
# points(scaf1482[99:147,8],pch=20,cex=1.2,col="gold")
# 
# abline(h=c(col),col=c("black","grey","red","darkorange","gold"),lwd=3,lty=3)
# 
# 
# 
# 
