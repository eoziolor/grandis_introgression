library(stringr)
library(dplyr)
library(gtools)

z<-read.table("~/analysis/data/dfst/zmerge_fst_pi_5kb",stringsAsFactors = FALSE)

names<- c("Scaf","start","end","bbvb","bbpb","bbsjsp","bbbnp","bbsp","bbgb",
          "vbpb","vbsjsp","vbbnp","vbsp","vbgb",
          "pbsjsp","pbbnp","pbsp","pbgb",
          "sjspbnp","sjspsp","sjspgb",
          "bnpsp","bnpgb",
          "spgb","keep")
colnames(z)<-names

###Calculating total distance from SP and GB for each pop for merged z statistic.
subw<-z[,25]>1
distsp<-cbind(z[subw,1:3],z[subw,8],z[subw,13],z[subw,17],z[subw,20],z[subw,22])
distgb<-cbind(z[subw,1:3],z[subw,9],z[subw,14],z[subw,18],z[subw,21],z[subw,23])

distnames<-c("Scaf","start","end","bb","vb","pb","sj","bnp")
colnames(distsp)<-distnames
colnames(distgb)<-distnames

colnam<-names(distsp)[4:8]
distsp2<-cbind(seq=seq(1:848015),distsp)
distgb2<-cbind(seq=seq(1:848015),distgb)

total_dist<-cbind(z[subw,1:3],distsp2[colnam]+distgb2[match(distsp2$seq,distgb2$seq),colnam])

pbs_dist<-matrix(nrow=848015,ncol=5)
for(i in 1:5)
{
  pbs_dist[,i]<-total_dist[,i+3]/2
}

pbs_dist<-cbind(z[subw,1:3],pbs_dist)
colnames(pbs_dist)<-distnames

col<-c()
for (i in 1:5){
  col[i]<-quantile(pbs_dist[,i+3],prob=.99,na.rm=TRUE)
}

###Overall outliers

all<-pbs_dist[,4]>col[1] & pbs_dist[,5]>col[2] & pbs_dist[,6]>col[3] & pbs_dist[,7]>col[4] & pbs_dist[,8]>col[5]
res<-pbs_dist[,4]>col[1] & pbs_dist[,5]>col[2] & pbs_dist[,6]>col[3] & pbs_dist[,7]<col[4] & pbs_dist[,8]<col[5]
interm<-pbs_dist[,4]<col[1] & pbs_dist[,5]<col[2] & pbs_dist[,6]<col[3] & pbs_dist[,7]>col[4] & pbs_dist[,8]>col[5]


write.table(na.omit(pbs_dist[all,1:3]),"~/analysis/data/dfst/zshared_outliers_all",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(na.omit(pbs_dist[res,1:3]),"~/analysis/data/dfst/zres_outliers_all",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(na.omit(pbs_dist[interm,1:3]),"~/analysis/data/dfst/zinterm_outliers_all",row.names=FALSE,col.names=FALSE,quote=FALSE)

####Chr arranged outliers

pbsd_chr<-pbs_dist %>% filter(str_detect(Scaf,"chr"))

all<-pbsd_chr[,4]>col[1] & pbsd_chr[,5]>col[2] & pbsd_chr[,6]>col[3] & pbsd_chr[,7]>col[4] & pbsd_chr[,8]>col[5]
res<-pbsd_chr[,4]>col[1] & pbsd_chr[,5]>col[2] & pbsd_chr[,6]>col[3] & pbsd_chr[,7]<col[4] & pbsd_chr[,8]<col[5]
interm<-pbsd_chr[,4]<col[1] & pbsd_chr[,5]<col[2] & pbsd_chr[,6]<col[3] & pbsd_chr[,7]>col[4] & pbsd_chr[,8]>col[5]


palette(c("black","grey","grey30","grey50"))
par(mfrow=c(5,1),mar=c(0,5,0,0))
plot(pbsd_chr[,4],pch=20,cex=.8,
     col=ifelse((all),"red",
                ifelse((res),"darkorange",
                       ifelse((interm),"gold",factor(pbsd_chr[,1])))),
     xlab="",xaxt='n',ylab="BB z score",cex.lab=1.7,cex.axis=.8,bty="n",ylim=c(-16,23))

plot(pbsd_chr[,5],pch=20,cex=.8,
     col=ifelse((all),"red",
                ifelse((res),"darkorange",
                       ifelse((interm),"gold",factor(pbsd_chr[,1])))),
     xlab="",xaxt='n',ylab="VB z score",cex.lab=1.7,cex.axis=.8,bty="n",ylim=c(-16,23))

plot(pbsd_chr[,6],pch=20,cex=.8,
     col=ifelse((all),"red",
                ifelse((res),"darkorange",
                       ifelse((interm),"gold",factor(pbsd_chr[,1])))),
     xlab="",xaxt='n',ylab="PB z score",cex.lab=1.7,cex.axis=.8,bty="n",ylim=c(-16,23))

plot(pbsd_chr[,7],pch=20,cex=.8,
     col=ifelse((all),"red",
                ifelse((res),"darkorange",
                       ifelse((interm),"gold",factor(pbsd_chr[,1])))),
     xlab="",xaxt='n',ylab="SJ z score",cex.lab=1.7,cex.axis=.8,bty="n",ylim=c(-16,23))

plot(pbsd_chr[,8],pch=20,cex=.8,
     col=ifelse((all),"red",
                ifelse((res),"darkorange",
                       ifelse((interm),"gold",factor(pbsd_chr[,1])))),
     xlab="Chromosomes 1-24",ylab="BNP z score",cex.lab=1.7,cex.axis=.8,bty="n",ylim=c(-16,23))

legend("topleft",legend=c("Shared","Resistant only","Intermediate only"),
       col=c("red","darkorange","gold"),pch=20,cex=1,y.intersp=.5)


#Regions

write.table(na.omit(pbsd_chr[all,1:3]),"~/analysis/data/dfst/zshared_outliers",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(na.omit(pbsd_chr[res,1:3]),"~/analysis/data/dfst/zres_outliers",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(na.omit(pbsd_chr[interm,1:3]),"~/analysis/data/dfst/zinterm_outliers",row.names=FALSE,col.names=FALSE,quote=FALSE)

#AHR b's

chr18<-pbs_dist[grep("chr18\\b",pbs_dist$Scaf),]
plot(chr18[13000:17499,4],pch=20,cex=.5,col="black",ylim=c(-10,20),ylab="PBS divergence stat",xlab="tail end of chromosome 18")
points(chr18[13000:17499,5],pch=20,cex=.5,col="grey",ylim=c(-1,1))
points(chr18[13000:17499,6],pch=20,cex=.5,col="red",ylim=c(-1,1))
points(chr18[13000:17499,7],pch=20,cex=.5,col="darkorange",ylim=c(-1,1))
points(chr18[13000:17499,8],pch=20,cex=.5,col="gold",ylim=c(-1,1))

legend("topleft",y=c(.05,.1),legend=c("BB","VB","PB","SJ","BNP"),pch=20,cex=1.2,
       col=c("black","grey","red","darkorange","gold"), x.intersp=.8,y.intersp=.8)

#arnt

chr8<-pbs_dist[grep("chr8\\b",pbs_dist$Scaf),]
plot(chr8[13000:17499,4],pch=20,cex=.5,col="black",ylim=c(-10,20),ylab="PBS divergence stat",xlab="tail end of chromosome 8")
points(chr8[13000:17499,5],pch=20,cex=.5,col="grey",ylim=c(-1,1))
points(chr8[13000:17499,6],pch=20,cex=.5,col="red",ylim=c(-1,1))
points(chr8[13000:17499,7],pch=20,cex=.5,col="darkorange",ylim=c(-1,1))
points(chr8[13000:17499,8],pch=20,cex=.5,col="gold",ylim=c(-1,1))

legend("topleft",y=c(.05,.1),legend=c("BB","VB","PB","SJ","BNP"),pch=20,cex=1.2,
       col=c("black","grey","red","darkorange","gold"), x.intersp=.8,y.intersp=.8)