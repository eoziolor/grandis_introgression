##########START LOWER FOR GRAPHS
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

total_dist<-cbind(pidiff[,1:3],distsp2[colnam]+distgb2[match(distsp2$seq,distgb2$seq),colnam])

# pbs_dist<-matrix(nrow=1026857,ncol=5)
# for (i in 1:5)
# {
#   pbs_dist[,i]<-(total_dist[,i+3]-pidiff[,24])/2
# }

pbs_dist<-matrix(nrow=1026857,ncol=5)
for(i in 1:5)
{
  pbs_dist[,i]<-total_dist[,i+3]/2
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
#pbs_dist[grep("chr1",pbs_dist$Scaf),]
pbs_dist<-na.omit(pbs_dist)

write.table(pbs_dist,"~/analysis/data/angsd/pbs_pi_all5kb", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


##############run from here for graphs
pbs_dist<-read.table("~/analysis/data/angsd/pbs_pi_all5kb",header=FALSE,stringsAsFactors = FALSE)
distnames<-c("Scaf","start","end","bb","vb","pb","sj","bnp")
colnames(pbs_dist)<-distnames
#quantile info (0.01)
#BB=-0.007121807
#VB=-0.005270508
#PB=-0.006081543
#SJ=-0.007564454 
#BNP=-0.004136112

cat ~/analysis/data/angsd/pbs_pi_all5kb | \
awk '$4<-0.007121807 || $5<-0.005270508 || $6<-0.006081543 || $7<-0.007564454 || $8<-0.004136112' | \
~/program/bedtools2/bin/bedtools merge -i stdin -d 50000 \
-c 4,4,5,5,6,6,7,7,8,8 \
-o min,count,min,count,min,count,min,count,min,count \
-g <(cut -f 1-2 ~/analysis/data/genome/unsplit_merge.fasta.fai) > ~/analysis/data/angsd/PBSoutliers_50kb_all.bed

pbs_out_temp<-read.table("~/analysis/data/angsd/PBSoutliers_50kb_all.bed",stringsAsFactors = FALSE) #loads a pbs vector with windows merged within 50kb of each other and with max and windows count statistics
names<-c("Scaf","start","end","BBmax","BBcount","VBmax","VBcount","PBmax","PBcount","SJmax","SJcount","BNPmax","BNPcount")
colnames(pbs_out_temp)<-names

pbs_out<-pbs_out_temp %>% filter(str_detect(Scaf,"chr"))

all<-pbs_out[,4]<collow[1] & pbs_out[,6]<collow[2] & pbs_out[,8]<collow[3] & pbs_out[,10]<collow[4] & pbs_out[,12]<collow[5]
res<-pbs_out[,4]<collow[1] & pbs_out[,6]<collow[2] & pbs_out[,8]<collow[3] & pbs_out[,10]>collow[4] & pbs_out[,12]>collow[5]
interm<-pbs_out[,4]>collow[1] & pbs_out[,6]>collow[2] & pbs_out[,8]>collow[3] & pbs_out[,10]<collow[4] & pbs_out[,12]<collow[5]
bbu<-pbs_out[,4]<collow[1] & pbs_out[,6]>collow[2] & pbs_out[,8]>collow[3] & pbs_out[,10]>collow[4] & pbs_out[,12]>collow[5]
vbu<-pbs_out[,4]>collow[1] & pbs_out[,6]<collow[2] & pbs_out[,8]>collow[3] & pbs_out[,10]>collow[4] & pbs_out[,12]>collow[5]
pbu<-pbs_out[,4]>collow[1] & pbs_out[,6]>collow[2] & pbs_out[,8]<collow[3] & pbs_out[,10]>collow[4] & pbs_out[,12]>collow[5]
sju<-pbs_out[,4]>collow[1] & pbs_out[,6]>collow[2] & pbs_out[,8]>collow[3] & pbs_out[,10]<collow[4] & pbs_out[,12]>collow[5]
bnpu<-pbs_out[,4]>collow[1] & pbs_out[,6]>collow[2] & pbs_out[,8]>collow[3] & pbs_out[,10]>collow[4] & pbs_out[,12]<collow[5]

write.table(pbs_dist[,1:3],"~/analysis/data/angsd/PBS_keep_5kb.bed",row.names = FALSE,col.names = FALSE,quote=FALSE)
write.table(pbs_out[all,1:3],"~/analysis/data/angsd/pbs_regions_sharedall.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(pbs_out[res,1:3],"~/analysis/data/angsd/pbs_regions_sharedres.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(pbs_out[interm,1:3],"~/analysis/data/angsd/pbs_regions_sharedinterm.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(pbs_out[bbu,1:3],"~/analysis/data/angsd/pbs_regions_sharedbbu.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(pbs_out[vbu,1:3],"~/analysis/data/angsd/pbs_regions_sharedvbu.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(pbs_out[pbu,1:3],"~/analysis/data/angsd/pbs_regions_sharedpbu.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(pbs_out[sju,1:3],"~/analysis/data/angsd/pbs_regions_sharedsju.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(pbs_out[bnpu,1:3],"~/analysis/data/angsd/pbs_regions_sharedbnpu.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)

source("http://bioconductor.org/biocLite.R")
biocLite()
library("rtracklayer")

bed1=import("~/analysis/data/angsd/PBS_keep_5kb.bed")

bedall=import("~/analysis/data/angsd/pbs_regions_sharedall.bed")
bed1overlall=bed1[bed1 %over% bedall]
hitsall<-findOverlaps(bedall,bed1)
allhit<-subjectHits(hitsall)

bedres=import("~/analysis/data/angsd/pbs_regions_sharedres.bed")
bed1overlres=bed1[bed1 %over% bedres]
hitsres<-findOverlaps(bedres,bed1)
reshit<-subjectHits(hitsres)

bedinterm=import("~/analysis/data/angsd/pbs_regions_sharedinterm.bed")
bed1overlinterm=bed1[bed1 %over% bedinterm]
hitsinterm<-findOverlaps(bedinterm,bed1)
intermhit<-subjectHits(hitsinterm)

bedbbu=import("~/analysis/data/angsd/pbs_regions_sharedbbu.bed")
bed1overlbbu=bed1[bed1 %over% bedbbu]
hitsbbu<-findOverlaps(bedbbu,bed1)
bbuhit<-subjectHits(hitsbbu)

bedvbu=import("~/analysis/data/angsd/pbs_regions_sharedvbu.bed")
bed1overlvbu=bed1[bed1 %over% bedvbu]
hitsvbu<-findOverlaps(bedvbu,bed1)
vbuhit<-subjectHits(hitsvbu)

bedpbu=import("~/analysis/data/angsd/pbs_regions_sharedpbu.bed")
bed1overlpbu=bed1[bed1 %over% bedpbu]
hitspbu<-findOverlaps(bedpbu,bed1)
pbuhit<-subjectHits(hitspbu)

bedsju=import("~/analysis/data/angsd/pbs_regions_sharedsju.bed")
bed1overlsju=bed1[bed1 %over% bedsju]
hitssju<-findOverlaps(bedsju,bed1)
sjuhit<-subjectHits(hitssju)

bedbnpu=import("~/analysis/data/angsd/pbs_regions_sharedbnpu.bed")
bed1overlbnpu=bed1[bed1 %over% bedbnpu]
hitsbnpu<-findOverlaps(bedbnpu,bed1)
bnpuhit<-subjectHits(hitsbnpu)

pbs_dist<-cbind(pbs_dist,0,0,0,0,0,0,0,0)
newn<-c("Scaf","start","end","BB","VB","PB","SJ","BNP","all","res","interm","bbu","vbu","pbu","sju","bnpu")
colnames(pbs_dist)<-newn
pbs_dist[allhit,"all"]<-pbs_dist[allhit,"all"]+1
pbs_dist[reshit,"res"]<-pbs_dist[reshit,"res"]+1
pbs_dist[intermhit,"interm"]<-pbs_dist[intermhit,"interm"]+1
pbs_dist[bbuhit,"bbu"]<-pbs_dist[bbuhit,"bbu"]+1
pbs_dist[vbuhit,"vbu"]<-pbs_dist[vbuhit,"vbu"]+1
pbs_dist[pbuhit,"pbu"]<-pbs_dist[pbuhit,"pbu"]+1
pbs_dist[sjuhit,"sju"]<-pbs_dist[sjuhit,"sju"]+1
pbs_dist[bnpuhit,"bnpu"]<-pbs_dist[bnpuhit,"bnpu"]+1

palette(c("grey50","grey70","black","grey30"))
par(mfrow=c(5,1),mar=c(0,3,0,0))
plot(pbs_dist[grep("chr1\\b",pbs_dist$Scaf),4],pch=20,cex=1.2,
     col=ifelse(pbs_dist[grep("chr1\\b",pbs_dist$Scaf),"all"]>0,"purple",
                ifelse(pbs_dist[grep("chr1\\b",pbs_dist$Scaf),"res"]>0,"black",
                       ifelse(pbs_dist[grep("chr1\\b",pbs_dist$Scaf),"interm"]>0,"firebrick2",
                              ifelse(pbs_dist[grep("chr1\\b",pbs_dist$Scaf),"bbu"]>0,"gold2",
                                     ifelse(pbs_dist[grep("chr1\\b",pbs_dist$Scaf),4]<collow[1],"green2",as.factor(pbs_dist[,1])))))),
     xlab="",xaxt='n',cex.lab=1,cex.axis=1.2,bty="n",yaxs="i",ylim=c(-.1,.1))

plot(pbs_dist[grep("chr1\\b",pbs_dist$Scaf),5],pch=20,cex=1.2,
     col=ifelse(pbs_dist[grep("chr1\\b",pbs_dist$Scaf),"all"]>0,"purple",
                ifelse(pbs_dist[grep("chr1\\b",pbs_dist$Scaf),"res"]>0,"black",
                       ifelse(pbs_dist[grep("chr1\\b",pbs_dist$Scaf),"interm"]>0,"firebrick2",
                              ifelse(pbs_dist[grep("chr1\\b",pbs_dist$Scaf),"vbu"]>0,"gold2",
                                     ifelse(pbs_dist[grep("chr1\\b",pbs_dist$Scaf),5]<collow[2],"green2",as.factor(pbs_dist[grep("chr1\\b",pbs_dist$Scaf),1])))))),
     xlab="",xaxt='n',cex.lab=1,cex.axis=1.2,bty="n",yaxs="i",ylim=c(-.1,.1))

plot(pbs_dist[grep("chr1\\b",pbs_dist$Scaf),6],pch=20,cex=1.2,
     col=ifelse(pbs_dist[grep("chr1\\b",pbs_dist$Scaf),"all"]>0,"purple",
                ifelse(pbs_dist[grep("chr1\\b",pbs_dist$Scaf),"res"]>0,"black",
                       ifelse(pbs_dist[grep("chr1\\b",pbs_dist$Scaf),"interm"]>0,"firebrick2",
                              ifelse(pbs_dist[grep("chr1\\b",pbs_dist$Scaf),"pbu"]>0,"gold2",
                                     ifelse(pbs_dist[grep("chr1\\b",pbs_dist$Scaf),6]<collow[3],"green2",as.factor(pbs_dist[grep("chr1\\b",pbs_dist$Scaf),1])))))),
     xlab="",xaxt='n',cex.lab=1,cex.axis=1.2,bty="n",yaxs="i",ylim=c(-.1,.1))

plot(pbs_dist[grep("chr1\\b",pbs_dist$Scaf),7],pch=20,cex=1.2,
     col=ifelse(pbs_dist[grep("chr1\\b",pbs_dist$Scaf),"all"]>0,"purple",
                ifelse(pbs_dist[grep("chr1\\b",pbs_dist$Scaf),"res"]>0,"black",
                       ifelse(pbs_dist[grep("chr1\\b",pbs_dist$Scaf),"interm"]>0,"firebrick2",
                              ifelse(pbs_dist[grep("chr1\\b",pbs_dist$Scaf),"sju"]>0,"gold2",
                                     ifelse(pbs_dist[grep("chr1\\b",pbs_dist$Scaf),7]<collow[4],"green2",as.factor(pbs_dist[grep("chr1\\b",pbs_dist$Scaf),1])))))),
     xlab="",xaxt='n',cex.lab=1,cex.axis=1.2,bty="n",yaxs="i",ylim=c(-.1,.1))

plot(pbs_dist[grep("chr1\\b",pbs_dist$Scaf),8],pch=20,cex=1.2,
     col=ifelse(pbs_dist[grep("chr1\\b",pbs_dist$Scaf),"all"]>0,"purple",
                ifelse(pbs_dist[grep("chr1\\b",pbs_dist$Scaf),"res"]>0,"black",
                       ifelse(pbs_dist[grep("chr1\\b",pbs_dist$Scaf),"interm"]>0,"firebrick2",
                              ifelse(pbs_dist[grep("chr1\\b",pbs_dist$Scaf),"bnpu"]>0,"gold2",
                                     ifelse(pbs_dist[grep("chr1\\b",pbs_dist$Scaf),8]<collow[5],"green2",as.factor(pbs_dist[grep("chr1\\b",pbs_dist$Scaf),1])))))),
     xlab="",xaxt='n',cex.lab=1,cex.axis=1.2,bty="n",yaxs="i",ylim=c(-.1,.1))



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

###AIP region (chr2:27525885-27626885)

chr2<-pbs_dist[grep("chr2\\b",pbs_dist$Scaf),]
plot(chr2[27500:27600,4],pch=20,cex=.5,col="black",ylim=c(-.1,.1))
points(chr2[27500:27600,5],pch=20,cex=.5,col="grey")
points(chr2[27500:27600,6],pch=20,cex=.5,col="red")
points(chr2[27500:27600,7],pch=20,cex=.5,col="darkorange")
points(chr2[27500:27600,8],pch=20,cex=.5,col="gold")

plot(chr2[27000:28000,4],pch=20,cex=.5,col="black",ylim=c(-.01,.005),
     xlab="Chr2:27000000:28000000 (in 5kb chunks)",
     ylab="average difference between res pi and references pi")
points(chr2[27000:28000,5],pch=20,cex=.5,col="grey")
points(chr2[27000:28000,6],pch=20,cex=.5,col="red")
points(chr2[27000:28000,7],pch=20,cex=.5,col="darkorange")
points(chr2[27000:28000,8],pch=20,cex=.5,col="gold")

legend("topright",y=c(.05,.1),legend=c("BB","VB","PB","SJSP","BNP"),pch=20,cex=1.2,
       col=c("black","grey","red","darkorange","gold"), x.intersp=.3,y.intersp=.6)
abline(h=0,col="purple",lty=2)


###Chromosome 10 ARNT region
chr10<-pbs_dist[grep("chr10\\b",pbs_dist$Scaf),]
plot(chr10[,4],pch=20,cex=.5,col="black",ylim=c(-.05,.05),ylab="pi difference (lower means lower pi in resistant)",xlab="chromosome 10")
points(chr10[,5],pch=20,cex=.5,col="grey")
points(chr10[,6],pch=20,cex=.5,col="red")

legend("topright",y=c(.05,.1),legend=c("BB","VB","PB"),pch=20,cex=1.2,
       col=c("black","grey","red"), x.intersp=.8,y.intersp=.8)


#Chromosome 18 AHR b's
chr18<-pbs_dist[grep("chr18\\b",pbs_dist$Scaf),]
plot(chr18[15000:20793,4],pch=20,cex=.5,col="black",ylim=c(-.05,.05),ylab="pi difference (lower means lower pi in resistant)",
     xlab="tail end of chromosome 18")
points(chr18[15000:20793,5],pch=20,cex=.5,col="grey")
points(chr18[15000:20793,6],pch=20,cex=.5,col="red")
points(chr18[15000:20793,7],pch=20,cex=.5,col="darkorange")
points(chr18[15000:20793,8],pch=20,cex=.5,col="gold")

legend("topleft",y=c(.05,.1),legend=c("BB","VB","PB","SJ","BNP"),pch=20,cex=1.2,
       col=c("black","grey","red","darkorange","gold"), x.intersp=.8,y.intersp=.8)

