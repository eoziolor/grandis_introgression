pbs<-read.table("~/analysis/data/fst/allpbs5kb",header=FALSE,stringsAsFactors = FALSE)
pbsname<-c("Scaf","start","end","BBpbs","VBpbs","PBpbs","SJpbs","BNPpbs","keep")
colnames(pbs)<-pbsname

col<-c()
for (i in 1:5){
  col[i]<-quantile(pbs[,i+3],prob=.99,na.rm=TRUE)
}

nsnps <-pbs[,"keep"]
subw<-nsnps>0
#chr<-read.table("~/analysis/fst/scripts/chr_colors",stringsAsFactors=FALSE,sep="\t")

library(stringr)
library(dplyr)
library(gtools)

#plot(pbs[subw,4],pch=20,cex=.5,col=factor(pbs[subw,1]))
#plot(pbs[subw,4],pch=20,cex=.5,col=chr[pbs[subw,1],2])
#legend('topright',legend=levels(mixedsort(pbs[,1])),col=1:2,cex=.5,pch=1)

pbsc<-pbs %>% filter(str_detect(Scaf,"chr"))

nsnps <-pbsc[,"keep"]
subwc<-nsnps>0

palette(c("black","grey","grey30","grey50"))
par(mfrow=c(5,1),mar=c(0,5,0,0))
plot(pbsc[subwc,4],pch=20,cex=1,col=ifelse((pbsc[subw,4]<col[1]),factor(pbsc[subwc,1]),"red"),
     xlab="",xaxt='n',ylab="BB (PBS)",cex.lab=1.7,cex.axis=.8,bty="n",
     ylim=c(-.5,4))
plot(pbsc[subwc,5],pch=20,cex=1,col=ifelse((pbsc[subw,5]<col[2]),factor(pbsc[subwc,1]),"red"),
     xlab="",xaxt='n',ylab="VB (PBS)",cex.lab=1.7,cex.axis=.8,bty="n",
     ylim=c(-.5,4))
plot(pbsc[subwc,6],pch=20,cex=1,col=ifelse((pbsc[subw,6]<col[3]),factor(pbsc[subwc,1]),"red"),
     xlab="",xaxt='n',ylab="PB (PBS)",cex.lab=1.7,cex.axis=.8,bty="n",
     ylim=c(-.5,4))
plot(pbsc[subwc,7],pch=20,cex=1,col=ifelse((pbsc[subw,7]<col[4]),factor(pbsc[subwc,1]),"red"),
     xlab="",xaxt='n',ylab="SJ (PBS)",cex.lab=1.7,cex.axis=.8,bty="n",
     ylim=c(-.5,4))
plot(pbsc[subwc,8],pch=20,cex=1,col=ifelse((pbsc[subw,8]<col[5]),factor(pbsc[subwc,1]),"red"),
     xlab="Chromosomes 1-24",ylab="BNP (PBS)",cex.lab=1.7,cex.axis=.8,bty="n",
     ylim=c(-.5,4))

#####Plotting common regions

all<-pbsc[subw,4]>col[1] & pbsc[subw,5]>col[2] & pbsc[subw,6]>col[3] & pbsc[subw,7]>col[4] & pbsc[subw,8]>col[5]
res<-pbsc[subw,4]>col[1] & pbsc[subw,5]>col[2] & pbsc[subw,6]>col[3] & pbsc[subw,7]<col[4] & pbsc[subw,8]<col[5]
interm<-pbsc[subw,4]<col[1] & pbsc[subw,5]<col[2] & pbsc[subw,6]<col[3] & pbsc[subw,7]>col[4] & pbsc[subw,8]>col[5]
bbu<-pbsc[subw,4]>col[1] & pbsc[subw,5]<col[2] & pbsc[subw,6]<col[3] & pbsc[subw,7]<col[4] & pbsc[subw,8]<col[5]
vbu<-pbsc[subw,4]<col[1] & pbsc[subw,5]>col[2] & pbsc[subw,6]<col[3] & pbsc[subw,7]<col[4] & pbsc[subw,8]<col[5]
pbu<-pbsc[subw,4]<col[1] & pbsc[subw,5]<col[2] & pbsc[subw,6]>col[3] & pbsc[subw,7]<col[4] & pbsc[subw,8]<col[5]
sju<-pbsc[subw,4]<col[1] & pbsc[subw,5]<col[2] & pbsc[subw,6]<col[3] & pbsc[subw,7]>col[4] & pbsc[subw,8]<col[5]
bnpu<-pbsc[subw,4]<col[1] & pbsc[subw,5]<col[2] & pbsc[subw,6]<col[3] & pbsc[subw,7]<col[4] & pbsc[subw,8]>col[5]


palette(c("grey50","grey70"))
par(mfrow=c(5,1),mar=c(0,3,0,0))
plot(pbsc[subwc,4],pch=20,cex=1.2,
     col=ifelse((all),"purple",
                ifelse((res),"black",
                       ifelse((interm),"firebrick2",
                              ifelse((bbu),"gold2",
                                     ifelse(pbsc[subwc,4]>col[1],"green2",sort(as.factor(pbsc[subwc,1]))))))),
     xlab="",xaxt='n',ylab="BB (PBS)",cex.lab=1,cex.axis=2.2,bty="n",ylim=c(-.5,3.8),xaxs="i",yaxs="i")

# legend("topright",legend=c("Shared by all adapted","Resistant only","Intermediate only","Shared (group non-specific)","Local"),
#        col=c("purple","black","firebrick2","green2","gold2"),pch=20,cex=1.8,y.intersp=.5,x.intersp=.8,bty='n')

plot(pbsc[subwc,5],pch=20,cex=1.2,
     col=ifelse((all),"purple",
                ifelse((res),"black",
                       ifelse((interm),"firebrick2",
                              ifelse((vbu),"gold2",
                                     ifelse(pbsc[subwc,5]>col[2],"green2",sort(as.factor(pbsc[subwc,1]))))))),
     xlab="",xaxt='n',ylab="VB (PBS)",cex.lab=1,cex.axis=2.2,bty="n",ylim=c(-.5,3.8),xaxs="i",yaxs="i")

plot(pbsc[subwc,6],pch=20,cex=1.2,
     col=ifelse((all),"purple",
                ifelse((res),"black",
                       ifelse((interm),"firebrick2",
                              ifelse((pbu),"gold2",
                                     ifelse(pbsc[subwc,6]>col[3],"green2",sort(as.factor(pbsc[subwc,1]))))))),
     xlab="",xaxt='n',ylab="PB (PBS)",cex.lab=1,cex.axis=2.2,bty="n",ylim=c(-.5,3.8),xaxs="i",yaxs="i")

plot(pbsc[subwc,7],pch=20,cex=1.2,
     col=ifelse((all),"purple",
                ifelse((res),"black",
                       ifelse((interm),"firebrick2",
                              ifelse((sju),"gold2",
                                     ifelse(pbsc[subwc,7]>col[4],"green2",sort(as.factor(pbsc[subwc,1]))))))),
     xlab="",xaxt='n',ylab="SJ (PBS)",cex.lab=1,cex.axis=2.2,bty="n",ylim=c(-.5,3.8),xaxs="i",yaxs="i")

plot(pbsc[subwc,8],pch=20,cex=1.2,
     col=ifelse((all),"purple",
                ifelse((res),"black",
                       ifelse((interm),"firebrick2",
                              ifelse((bnpu),"gold2",
                                     ifelse(pbsc[subwc,8]>col[5],"green2",sort(as.factor(pbsc[subwc,1]))))))),
     xlab="",xaxt='n',ylab="BNP (PBS)",cex.lab=1,cex.axis=2.2,bty="n",ylim=c(-.5,3.8),xaxs="i",yaxs="i")


###chr8
# 
# chr2<-pbsc[grep("chr2\\b",pbsc$Scaf),]
# plot(chr2[,4],pch=20,cex=.5)
# 
# chr24<-pbsc[grep("chr24\\b",pbsc$Scaf),]
# plot(chr24[,4],pch=20,cex=.5)
