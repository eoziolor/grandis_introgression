library("RColorBrewer")
library("lattice")
library("gplots")
#Loaidng list of fst files into a list object ----
fs <- list.files("~/analysis/data/fst/raw/", "*fst.1kb.bed",full.names=TRUE)

fst <- list()

for (i in 1:21){
  fst[[i]] <- read.table(fs[i],stringsAsFactors=FALSE)
  fst[[i]][,4] <- as.numeric(fst[[i]][,4])
}

nfs <- gsub(".*\\/","",fs)
nfs <- gsub(".fst.*","",nfs)
names(fst)<-nfs

#loading the pbs function----
pbs <- function(t1,t2,c12){
  
  t1 <- -log(1-t1)
  t2 <- -log(1-t2)
  c12 <- -log(1-c12)
  
  stat <- (t1 + t2 - c12)/2
  return(stat)
}

#selecting sites that ahave a minimum representation of 200 snps per region----
nsnps <-fst[[1]][,5]

for (i in 2:21){
  
  nsnps <- nsnps + fst[[i]][,5]
}

nsnps <- nsnps/21

subw <- nsnps > 4

#calculating FST for all

fstl<-c()
for(i in 1:21){
  fstl[i]<-mean(fst[[i]][subw,4],na.rm=TRUE)
}

names(fstl)<-names(fst)

fsth<-matrix(nrow = 7,ncol=7)
colnames(fsth)<-c("BB","VB","PB","SJ","BNP","SP","GB")
rownames(fsth)<-c("BB","VB","PB","SJ","BNP","SP","GB")
fsth[1,]<-c(0,fstl["BB.VB"],fstl["BB.PB"],fstl["BB.SJ"],fstl["BB.BNP"],fstl["BB.SP"],fstl["BB.GB"])
fsth[2,]<-c(fstl["BB.VB"],0,fstl["VB.PB"],fstl["VB.SJ"],fstl["VB.BNP"],fstl["VB.SP"],fstl["VB.GB"])
fsth[3,]<-c(fstl["BB.PB"],fstl["VB.PB"],0,fstl["PB.SJ"],fstl["PB.BNP"],fstl["PB.SP"],fstl["PB.GB"])
fsth[4,]<-c(fstl["BB.SJ"],fstl["VB.SJ"],fstl["PB.SJ"],0,fstl["SJ.BNP"],fstl["SJ.SP"],fstl["SJ.GB"])
fsth[5,]<-c(fstl["BB.BNP"],fstl["VB.BNP"],fstl["PB.BNP"],fstl["SJ.BNP"],0,fstl["BNP.SP"],fstl["BNP.GB"])
fsth[6,]<-c(fstl["BB.SP"],fstl["VB.SP"],fstl["PB.SP"],fstl["SJ.SP"],fstl["BNP.SP"],0,fstl["GB.SP"])
fsth[7,]<-c(fstl["BB.GB"],fstl["VB.GB"],fstl["PB.GB"],fstl["SJ.GB"],fstl["BNP.GB"],fstl["GB.SP"],0)

levelplot(fsth,aspect="iso",col.regions=brewer.pal(9,"YlOrRd"),scale=list(x=list(rot=45)),cuts=8)

#Doing pbs on all pops----

BBpbs <- pbs(fst[["BB.GB"]][,4],fst[["BB.SP"]][,4],fst[["GB.SP"]][,4])

VBpbs <- pbs(fst[["VB.GB"]][,4],fst[["VB.SP"]][,4],fst[["GB.SP"]][,4])

PBpbs <- pbs(fst[["PB.GB"]][,4],fst[["PB.SP"]][,4],fst[["GB.SP"]][,4])

SJpbs <- pbs(fst[["SJ.GB"]][,4],fst[["SJ.SP"]][,4],fst[["GB.SP"]][,4])

BNPpbs <- pbs(fst[["BNP.GB"]][,4],fst[["BNP.SP"]][,4],fst[["GB.SP"]][,4])


Allpbs <- cbind(
  fst[[1]][,1:3],
  BB = BBpbs,
  VB = VBpbs,
  PB = PBpbs,
  SJ = SJpbs,
  BNP = BNPpbs,
  keep = as.numeric(subw))

write.table(Allpbs,
            file="~/analysis/data/fst/allpbs1kb",
            sep = "\t", 
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)


plot(BBpbs[subw],pch=20,cex=.5, col=factor(Allpbs[,1]))
quantile(BBpbs[subw],probs=.99,na.rm=TRUE)
abline(h=.1804205)


####Gene comparisons

pbs<-read.table("~/analysis/data/fst/allpbs1kb",header=FALSE)
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

###Chr10 arnt

chr10<-pbs[grep("chr10\\b",pbs$Scaf),]
plot(chr10[,4],pch=20,cex=.5,col="black",ylim=c(-1,1),ylab="PBS divergence stat",xlab="chromosome 10")
points(chr10[,5],pch=20,cex=.5,col="grey",ylim=c(-1,1))
points(chr10[,6],pch=20,cex=.5,col="red",ylim=c(-1,1))
points(chr10[,7],pch=20,cex=.5,col="darkorange",ylim=c(-1,1))
points(chr10[,8],pch=20,cex=.5,col="gold",ylim=c(-1,1))

legend("topright",y=c(.05,.1),legend=c("BB","VB","PB","SJ","BNP"),pch=20,cex=1.2,
       col=c("black","grey","red","darkorange","gold"), x.intersp=.8,y.intersp=.8)

###Chr18 arnt

chr18<-pbs[grep("chr18\\b",pbs$Scaf),]
plot(chr18[15000:20793,4],pch=20,cex=.5,col="black",ylim=c(-1,1),ylab="PBS divergence stat",xlab="tail end of chromosome 18")
points(chr18[15000:20793,5],pch=20,cex=.5,col="grey",ylim=c(-1,1))
points(chr18[15000:20793,6],pch=20,cex=.5,col="red",ylim=c(-1,1))
points(chr18[15000:20793,7],pch=20,cex=.5,col="darkorange",ylim=c(-1,1))
points(chr18[15000:20793,8],pch=20,cex=.5,col="gold",ylim=c(-1,1))

legend("topleft",y=c(.05,.1),legend=c("BB","VB","PB","SJ","BNP"),pch=20,cex=1.2,
       col=c("black","grey","red","darkorange","gold"), x.intersp=.8,y.intersp=.8)

################################################AHR2a

###finding neutral regions to do admixture on----

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
