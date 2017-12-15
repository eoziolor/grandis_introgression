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

pbsc<-pbs_dist %>% filter(str_detect(Scaf,"chr"))

all<-pbsc[,4]>col[1] & pbsc[,5]>col[2] & pbsc[,6]>col[3] & pbsc[,7]>col[4] & pbsc[,8]>col[5]
res<-pbsc[,4]>col[1] & pbsc[,5]>col[2] & pbsc[,6]>col[3] & pbsc[,7]<col[4] & pbsc[,8]<col[5]
interm<-pbsc[,4]<col[1] & pbsc[,5]<col[2] & pbsc[,6]<col[3] & pbsc[,7]>col[4] & pbsc[,8]>col[5]
bbu<-pbsc[,4]>col[1] & pbsc[,5]<col[2] & pbsc[,6]<col[3] & pbsc[,7]<col[4] & pbsc[,8]<col[5]
vbu<-pbsc[,4]<col[1] & pbsc[,5]>col[2] & pbsc[,6]<col[3] & pbsc[,7]<col[4] & pbsc[,8]<col[5]
pbu<-pbsc[,4]<col[1] & pbsc[,5]<col[2] & pbsc[,6]>col[3] & pbsc[,7]<col[4] & pbsc[,8]<col[5]
sju<-pbsc[,4]<col[1] & pbsc[,5]<col[2] & pbsc[,6]<col[3] & pbsc[,7]>col[4] & pbsc[,8]<col[5]
bnpu<-pbsc[,4]<col[1] & pbsc[,5]<col[2] & pbsc[,6]<col[3] & pbsc[,7]<col[4] & pbsc[,8]>col[5]


palette(c("grey50","grey70"))
par(mfrow=c(5,1),mar=c(0,3,0,0))
plot(pbsc[,4],pch=20,cex=1.2,
     col=ifelse((all),"purple",
                ifelse((res),"black",
                       ifelse((interm),"firebrick2",
                              ifelse((bbu),"gold2",
                                     ifelse(pbsc[,4]>col[1],"green2",sort(as.factor(pbsc[,1]))))))),
     xlab="",xaxt='n',ylab="BB (PBS)",cex.lab=1,cex.axis=2.2,bty="n",ylim=c(-16,23),xaxs="i",yaxs="i")

# legend("topright",legend=c("Shared by all adapted","Resistant only","Intermediate only","Shared (group non-specific)","Local"),
#        col=c("purple","black","firebrick2","green2","gold2"),pch=20,cex=1.8,y.intersp=.5,x.intersp=.8,bty='n')

plot(pbsc[,5],pch=20,cex=1.2,
     col=ifelse((all),"purple",
                ifelse((res),"black",
                       ifelse((interm),"firebrick2",
                              ifelse((vbu),"gold2",
                                     ifelse(pbsc[,5]>col[2],"green2",sort(as.factor(pbsc[,1]))))))),
     xlab="",xaxt='n',ylab="VB (PBS)",cex.lab=1,cex.axis=2.2,bty="n",ylim=c(-16,23),xaxs="i",yaxs="i")

plot(pbsc[,6],pch=20,cex=1.2,
     col=ifelse((all),"purple",
                ifelse((res),"black",
                       ifelse((interm),"firebrick2",
                              ifelse((pbu),"gold2",
                                     ifelse(pbsc[,6]>col[3],"green2",sort(as.factor(pbsc[,1]))))))),
     xlab="",xaxt='n',ylab="PB (PBS)",cex.lab=1,cex.axis=2.2,bty="n",ylim=c(-16,23),xaxs="i",yaxs="i")

plot(pbsc[,7],pch=20,cex=1.2,
     col=ifelse((all),"purple",
                ifelse((res),"black",
                       ifelse((interm),"firebrick2",
                              ifelse((sju),"gold2",
                                     ifelse(pbsc[,7]>col[4],"green2",sort(as.factor(pbsc[,1]))))))),
     xlab="",xaxt='n',ylab="SJ (PBS)",cex.lab=1,cex.axis=2.2,bty="n",ylim=c(-16,23),xaxs="i",yaxs="i")

plot(pbsc[,8],pch=20,cex=1.2,
     col=ifelse((all),"purple",
                ifelse((res),"black",
                       ifelse((interm),"firebrick2",
                              ifelse((bnpu),"gold2",
                                     ifelse(pbsc[,8]>col[5],"green2",sort(as.factor(pbsc[,1]))))))),
     xlab="",xaxt='n',ylab="BNP (PBS)",cex.lab=1,cex.axis=2.2,bty="n",ylim=c(-16,23),xaxs="i",yaxs="i")


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