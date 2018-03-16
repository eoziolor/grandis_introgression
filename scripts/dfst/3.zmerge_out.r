library(stringr)
library(dplyr)
library(gtools)

#reading in z values table----
z<-read.table("~/analysis/data/dfst/zmerge_fst_pi_5kb",stringsAsFactors = FALSE)

names<- c("Scaf","start","end","bbvb","bbpb","bbsjsp","bbbnp","bbsp","bbgb",
          "vbpb","vbsjsp","vbbnp","vbsp","vbgb",
          "pbsjsp","pbbnp","pbsp","pbgb",
          "sjspbnp","sjspsp","sjspgb",
          "bnpsp","bnpgb",
          "spgb","keep")
colnames(z)<-names

###Calculating total distance from SP and GB for each pop for merged z statistic----
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

write.table(pbs_dist,"~/analysis/data/dfst/zpbs_5kb.bed",row.names = FALSE,col.names = FALSE,quote = FALSE,sep="\t")

###Overall outliers----

all<-pbs_dist[,4]>col[1] & pbs_dist[,5]>col[2] & pbs_dist[,6]>col[3] & pbs_dist[,7]>col[4] & pbs_dist[,8]>col[5]
res<-pbs_dist[,4]>col[1] & pbs_dist[,5]>col[2] & pbs_dist[,6]>col[3] & pbs_dist[,7]<col[4] & pbs_dist[,8]<col[5]
interm<-pbs_dist[,4]<col[1] & pbs_dist[,5]<col[2] & pbs_dist[,6]<col[3] & pbs_dist[,7]>col[4] & pbs_dist[,8]>col[5]


write.table(na.omit(pbs_dist[all,1:3]),"~/analysis/data/dfst/zshared_outliers_all",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(na.omit(pbs_dist[res,1:3]),"~/analysis/data/dfst/zres_outliers_all",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(na.omit(pbs_dist[interm,1:3]),"~/analysis/data/dfst/zinterm_outliers_all",row.names=FALSE,col.names=FALSE,quote=FALSE)

####Chr arranged outliers----

pbs_out_temp<-read.table("~/analysis/data/dfst/zregions_max.bed",stringsAsFactors = FALSE) #loads a pbs vector with windows merged within 50kb of each other and with max and windows count statistics
names<-c("Scaf","start","end","BBmax","BBcount","VBmax","VBcount","PBmax","PBcount","SJmax","SJcount","BNPmax","BNPcount")
colnames(pbs_out_temp)<-names

pbs_out<-pbs_out_temp %>% filter(str_detect(Scaf,"chr"))

pbsc<-pbs_dist %>% filter(str_detect(Scaf,"chr"))

all<-pbs_out[,4]>col[1] & pbs_out[,6]>col[2] & pbs_out[,8]>col[3] & pbs_out[,10]>col[4] & pbs_out[,12]>col[5]
res<-pbs_out[,4]>col[1] & pbs_out[,6]>col[2] & pbs_out[,8]>col[3] & pbs_out[,10]<col[4] & pbs_out[,12]<col[5]
interm<-pbs_out[,4]<col[1] & pbs_out[,6]<col[2] & pbs_out[,8]<col[3] & pbs_out[,10]>col[4] & pbs_out[,12]>col[5]
bbu<-pbs_out[,4]>col[1] & pbs_out[,6]<col[2] & pbs_out[,8]<col[3] & pbs_out[,10]<col[4] & pbs_out[,12]<col[5]
vbu<-pbs_out[,4]<col[1] & pbs_out[,6]>col[2] & pbs_out[,8]<col[3] & pbs_out[,10]<col[4] & pbs_out[,12]<col[5]
pbu<-pbs_out[,4]<col[1] & pbs_out[,6]<col[2] & pbs_out[,8]>col[3] & pbs_out[,10]<col[4] & pbs_out[,12]<col[5]
sju<-pbs_out[,4]<col[1] & pbs_out[,6]<col[2] & pbs_out[,8]<col[3] & pbs_out[,10]>col[4] & pbs_out[,12]<col[5]
bnpu<-pbs_out[,4]<col[1] & pbs_out[,6]<col[2] & pbs_out[,8]<col[3] & pbs_out[,10]<col[4] & pbs_out[,12]>col[5]

write.table(pbsc[,1:3],"~/analysis/data/dfst/PBS_keep_5kb.bed",row.names = FALSE,col.names = FALSE,quote=FALSE)
write.table(pbs_out[all,1:3],"~/analysis/data/dfst/pbs_regions_sharedall.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(pbs_out[res,1:3],"~/analysis/data/dfst/pbs_regions_sharedres.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(pbs_out[interm,1:3],"~/analysis/data/dfst/pbs_regions_sharedinterm.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(pbs_out[bbu,1:3],"~/analysis/data/dfst/pbs_regions_sharedbbu.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(pbs_out[vbu,1:3],"~/analysis/data/dfst/pbs_regions_sharedvbu.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(pbs_out[pbu,1:3],"~/analysis/data/dfst/pbs_regions_sharedpbu.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(pbs_out[sju,1:3],"~/analysis/data/dfst/pbs_regions_sharedsju.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(pbs_out[bnpu,1:3],"~/analysis/data/dfst/pbs_regions_sharedbnpu.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)

#Finding regions of overlap----

# source("http://bioconductor.org/biocLite.R")
# biocLite()
library("rtracklayer")


bed1=import("~/analysis/data/dfst/PBS_keep_5kb.bed")

bedall=import("~/analysis/data/dfst/pbs_regions_sharedall.bed")
bed1overlall=bed1[bed1 %over% bedall]
hitsall<-findOverlaps(bedall,bed1)
allhit<-subjectHits(hitsall)

bedres=import("~/analysis/data/dfst/pbs_regions_sharedres.bed")
bed1overlres=bed1[bed1 %over% bedres]
hitsres<-findOverlaps(bedres,bed1)
reshit<-subjectHits(hitsres)

bedinterm=import("~/analysis/data/dfst/pbs_regions_sharedinterm.bed")
bed1overlinterm=bed1[bed1 %over% bedinterm]
hitsinterm<-findOverlaps(bedinterm,bed1)
intermhit<-subjectHits(hitsinterm)

bedbbu=import("~/analysis/data/dfst/pbs_regions_sharedbbu.bed")
bed1overlbbu=bed1[bed1 %over% bedbbu]
hitsbbu<-findOverlaps(bedbbu,bed1)
bbuhit<-subjectHits(hitsbbu)

bedvbu=import("~/analysis/data/dfst/pbs_regions_sharedvbu.bed")
bed1overlvbu=bed1[bed1 %over% bedvbu]
hitsvbu<-findOverlaps(bedvbu,bed1)
vbuhit<-subjectHits(hitsvbu)

bedpbu=import("~/analysis/data/dfst/pbs_regions_sharedpbu.bed")
bed1overlpbu=bed1[bed1 %over% bedpbu]
hitspbu<-findOverlaps(bedpbu,bed1)
pbuhit<-subjectHits(hitspbu)

bedsju=import("~/analysis/data/dfst/pbs_regions_sharedsju.bed")
bed1overlsju=bed1[bed1 %over% bedsju]
hitssju<-findOverlaps(bedsju,bed1)
sjuhit<-subjectHits(hitssju)

bedbnpu=import("~/analysis/data/dfst/pbs_regions_sharedbnpu.bed")
bed1overlbnpu=bed1[bed1 %over% bedbnpu]
hitsbnpu<-findOverlaps(bedbnpu,bed1)
bnpuhit<-subjectHits(hitsbnpu)

#putting those regions into a dataframe----
pbsc<-cbind(pbsc,0,0,0,0,0,0,0,0)
newn<-c("Scaf","start","end","BB","VB","PB","SJ","BNP","all","res","interm","bbu","vbu","pbu","sju","bnpu")
colnames(pbsc)<-newn
pbsc[allhit,"all"]<-pbsc[allhit,"all"]+1
pbsc[reshit,"res"]<-pbsc[reshit,"res"]+1
pbsc[intermhit,"interm"]<-pbsc[intermhit,"interm"]+1
pbsc[bbuhit,"bbu"]<-pbsc[bbuhit,"bbu"]+1
pbsc[vbuhit,"vbu"]<-pbsc[vbuhit,"vbu"]+1
pbsc[pbuhit,"pbu"]<-pbsc[pbuhit,"pbu"]+1
pbsc[sjuhit,"sju"]<-pbsc[sjuhit,"sju"]+1
pbsc[bnpuhit,"bnpu"]<-pbsc[bnpuhit,"bnpu"]+1

#plotting those results by using the pbs_out vector. Have to find a way to intersect it with a region----
palette(c("grey50","grey70","black","grey30"))
par(mfrow=c(5,1),mar=c(0,3,0,0))
plot(pbsc[,4],pch=20,cex=1.2,
     col=ifelse(pbsc[,"all"]>0,"purple",
                ifelse(pbsc[,"res"]>0,"black",
                       ifelse(pbsc[,"interm"]>0,"firebrick2",
                              ifelse(pbsc[,"bbu"]>0,"gold2",
                                     ifelse(pbsc[,4]>col[1],"green2",sort(as.factor(pbsc[,1]))))))),
     xlab="",xaxt='n',cex.lab=1,cex.axis=2.2,bty="n",ylim=c(-16,23),xaxs="i",yaxs="i")

plot(pbsc[,5],pch=20,cex=1.2,
     col=ifelse(pbsc[,"all"]>0,"purple",
                ifelse(pbsc[,"res"]>0,"black",
                       ifelse(pbsc[,"interm"]>0,"firebrick2",
                              ifelse(pbsc[,"vbu"]>0,"gold2",
                                     ifelse(pbsc[,5]>col[2],"green2",sort(as.factor(pbsc[,1]))))))),
     xlab="",xaxt='n',cex.lab=1,cex.axis=2.2,bty="n",ylim=c(-16,23),xaxs="i",yaxs="i")

plot(pbsc[,6],pch=20,cex=1.2,
     col=ifelse(pbsc[,"all"]>0,"purple",
                ifelse(pbsc[,"res"]>0,"black",
                       ifelse(pbsc[,"interm"]>0,"firebrick2",
                              ifelse(pbsc[,"pbu"]>0,"gold2",
                                     ifelse(pbsc[,6]>col[3],"green2",sort(as.factor(pbsc[,1]))))))),
     xlab="",xaxt='n',cex.lab=1,cex.axis=2.2,bty="n",ylim=c(-16,23),xaxs="i",yaxs="i")

plot(pbsc[,7],pch=20,cex=1.2,
     col=ifelse(pbsc[,"all"]>0,"purple",
                ifelse(pbsc[,"res"]>0,"black",
                       ifelse(pbsc[,"interm"]>0,"firebrick2",
                              ifelse(pbsc[,"sju"]>0,"gold2",
                                     ifelse(pbsc[,7]>col[4],"green2",sort(as.factor(pbsc[,1]))))))),
     xlab="",xaxt='n',cex.lab=1,cex.axis=2.2,bty="n",ylim=c(-16,23),xaxs="i",yaxs="i")

plot(pbsc[,8],pch=20,cex=1.2,
     col=ifelse(pbsc[,"all"]>0,"purple",
                ifelse(pbsc[,"res"]>0,"black",
                       ifelse(pbsc[,"interm"]>0,"firebrick2",
                              ifelse(pbsc[,"bnpu"]>0,"gold2",
                                     ifelse(pbsc[,8]>col[5],"green2",sort(as.factor(pbsc[,1]))))))),
     xlab="",xaxt='n',cex.lab=1,cex.axis=2.2,bty="n",ylim=c(-16,23),xaxs="i",yaxs="i")


write.table(pbsc[allhit,1:3],"~/analysis/data/dfst/zregions_split5kb.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)

###Cleaner, by removing local regions----

palette(c("grey40","grey60","grey80","grey20"))
par(mfrow=c(5,1),mar=c(0,3,0,0))
plot(pbsc[,4],pch=20,cex=1.2,
     col=ifelse(pbsc[,"all"]>0,"purple",
                ifelse(pbsc[,"res"]>0,"black",
                       ifelse(pbsc[,"interm"]>0,"firebrick2",sort(as.factor(pbsc[,1]))))),
     xlab="",xaxt='n',cex.lab=1,cex.axis=2.2,bty="n",ylim=c(-16,23),yaxs="i")

plot(pbsc[,5],pch=20,cex=1.2,
     col=ifelse(pbsc[,"all"]>0,"purple",
                ifelse(pbsc[,"res"]>0,"black",
                       ifelse(pbsc[,"interm"]>0,"firebrick2",sort(as.factor(pbsc[,1]))))),
     xlab="",xaxt='n',cex.lab=1,cex.axis=2.2,bty="n",ylim=c(-16,23),yaxs="i")

plot(pbsc[,6],pch=20,cex=1.2,
     col=ifelse(pbsc[,"all"]>0,"purple",
                ifelse(pbsc[,"res"]>0,"black",
                       ifelse(pbsc[,"interm"]>0,"firebrick2",sort(as.factor(pbsc[,1]))))),
     xlab="",xaxt='n',cex.lab=1,cex.axis=2.2,bty="n",ylim=c(-16,23),yaxs="i")

plot(pbsc[,7],pch=20,cex=1.2,
     col=ifelse(pbsc[,"all"]>0,"purple",
                ifelse(pbsc[,"res"]>0,"black",
                       ifelse(pbsc[,"interm"]>0,"firebrick2",sort(as.factor(pbsc[,1]))))),
     xlab="",xaxt='n',cex.lab=1,cex.axis=2.2,bty="n",ylim=c(-16,23),yaxs="i")

plot(pbsc[,8],pch=20,cex=1.2,
     col=ifelse(pbsc[,"all"]>0,"purple",
                ifelse(pbsc[,"res"]>0,"black",
                       ifelse(pbsc[,"interm"]>0,"firebrick2",sort(as.factor(pbsc[,1]))))),
     xlab="",xaxt='n',cex.lab=1,cex.axis=2.2,bty="n",ylim=c(-16,23),yaxs="i")

###plotting distribution of outliers----
pbsc<-pbs_dist %>% filter(str_detect(Scaf,"chr"))
pbsc<-cbind(pbsc,0,0,0,0,0,0,0,0)
newn<-c("Scaf","start","end","BB","VB","PB","SJ","BNP","all","res","interm","bbu","vbu","pbu","sju","bnpu")
colnames(pbsc)<-newn
pbsc<-na.omit(pbsc)

all<-pbsc[,4]>col[1] & pbsc[,5]>col[2] & pbsc[,6]>col[3] & pbsc[,7]>col[4] & pbsc[,8]>col[5]
res<-pbsc[,4]>col[1] & pbsc[,5]>col[2] & pbsc[,6]>col[3] & pbsc[,7]<col[4] & pbsc[,8]<col[5]
interm<-pbsc[,4]<col[1] & pbsc[,5]<col[2] & pbsc[,6]<col[3] & pbsc[,7]>col[4] & pbsc[,8]>col[5]
bbu<-pbsc[,4]>col[1] & pbsc[,5]<col[2] & pbsc[,6]<col[3] & pbsc[,7]<col[4] & pbsc[,8]<col[5]
vbu<-pbsc[,4]<col[1] & pbsc[,5]>col[2] & pbsc[,6]<col[3] & pbsc[,7]<col[4] & pbsc[,8]<col[5]
pbu<-pbsc[,4]<col[1] & pbsc[,5]<col[2] & pbsc[,6]>col[3] & pbsc[,7]<col[4] & pbsc[,8]<col[5]
sju<-pbsc[,4]<col[1] & pbsc[,5]<col[2] & pbsc[,6]<col[3] & pbsc[,7]>col[4] & pbsc[,8]<col[5]
bnpu<-pbsc[,4]<col[1] & pbsc[,5]<col[2] & pbsc[,6]<col[3] & pbsc[,7]<col[4] & pbsc[,8]>col[5]

write.table(pbsc[,1:3],"~/analysis/data/dfst/PBS_keep_5kb.bed",row.names = FALSE,col.names = FALSE,quote=FALSE)
write.table(na.omit(pbsc[all,1:3]),"~/analysis/data/dfst/pbs_regions_sharedall.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(na.omit(pbsc[res,1:3]),"~/analysis/data/dfst/pbs_regions_sharedres.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(na.omit(pbsc[interm,1:3]),"~/analysis/data/dfst/pbs_regions_sharedinterm.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(na.omit(pbsc[bbu,1:3]),"~/analysis/data/dfst/pbs_regions_sharedbbu.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(na.omit(pbsc[vbu,1:3]),"~/analysis/data/dfst/pbs_regions_sharedvbu.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(na.omit(pbsc[pbu,1:3]),"~/analysis/data/dfst/pbs_regions_sharedpbu.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(na.omit(pbsc[sju,1:3]),"~/analysis/data/dfst/pbs_regions_sharedsju.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(na.omit(pbsc[bnpu,1:3]),"~/analysis/data/dfst/pbs_regions_sharedbnpu.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)

library("rtracklayer")


bed1=import("~/analysis/data/dfst/PBS_keep_5kb.bed")

bedall=import("~/analysis/data/dfst/pbs_regions_sharedall.bed")
bed1overlall=bed1[bed1 %over% bedall]
hitsall<-findOverlaps(bedall,bed1,type="equal")
allhit<-subjectHits(hitsall)

bedres=import("~/analysis/data/dfst/pbs_regions_sharedres.bed")
bed1overlres=bed1[bed1 %over% bedres]
hitsres<-findOverlaps(bedres,bed1,type="equal")
reshit<-subjectHits(hitsres)

bedinterm=import("~/analysis/data/dfst/pbs_regions_sharedinterm.bed")
bed1overlinterm=bed1[bed1 %over% bedinterm]
hitsinterm<-findOverlaps(bedinterm,bed1,type="equal")
intermhit<-subjectHits(hitsinterm)

bedbbu=import("~/analysis/data/dfst/pbs_regions_sharedbbu.bed")
bed1overlbbu=bed1[bed1 %over% bedbbu]
hitsbbu<-findOverlaps(bedbbu,bed1,type="equal")
bbuhit<-subjectHits(hitsbbu)

bedvbu=import("~/analysis/data/dfst/pbs_regions_sharedvbu.bed")
bed1overlvbu=bed1[bed1 %over% bedvbu]
hitsvbu<-findOverlaps(bedvbu,bed1,type="equal")
vbuhit<-subjectHits(hitsvbu)

bedpbu=import("~/analysis/data/dfst/pbs_regions_sharedpbu.bed")
bed1overlpbu=bed1[bed1 %over% bedpbu]
hitspbu<-findOverlaps(bedpbu,bed1,type="equal")
pbuhit<-subjectHits(hitspbu)

bedsju=import("~/analysis/data/dfst/pbs_regions_sharedsju.bed")
bed1overlsju=bed1[bed1 %over% bedsju]
hitssju<-findOverlaps(bedsju,bed1,type="equal")
sjuhit<-subjectHits(hitssju)

bedbnpu=import("~/analysis/data/dfst/pbs_regions_sharedbnpu.bed")
bed1overlbnpu=bed1[bed1 %over% bedbnpu]
hitsbnpu<-findOverlaps(bedbnpu,bed1,type="equal")
bnpuhit<-subjectHits(hitsbnpu)

pbsc<-cbind(pbsc,0,0,0,0,0,0,0,0)
newn<-c("Scaf","start","end","BB","VB","PB","SJ","BNP","all","res","interm","bbu","vbu","pbu","sju","bnpu")
colnames(pbsc)<-newn
pbsc[allhit,"all"]<-pbsc[allhit,"all"]+1
pbsc[reshit,"res"]<-pbsc[reshit,"res"]+1
pbsc[intermhit,"interm"]<-pbsc[intermhit,"interm"]+1
pbsc[bbuhit,"bbu"]<-pbsc[bbuhit,"bbu"]+1
pbsc[vbuhit,"vbu"]<-pbsc[vbuhit,"vbu"]+1
pbsc[pbuhit,"pbu"]<-pbsc[pbuhit,"pbu"]+1
pbsc[sjuhit,"sju"]<-pbsc[sjuhit,"sju"]+1
pbsc[bnpuhit,"bnpu"]<-pbsc[bnpuhit,"bnpu"]+1

par(mfrow=c(2,3),mar=c(4,4,0,0))
plot(pbsc[,"BB"],pbsc[,"VB"],pch=20,cex=.7,
     col=ifelse(pbsc[,"all"]>0,"purple",
                ifelse(pbsc[,"interm"]>0,"red",
                       ifelse(pbsc[,"res"],"black",NA))),bty='l',
     xlab="VB z values",ylab="BB z values")
abline(h=0,v=0)

plot(pbsc[,"BB"],pbsc[,"PB"],pch=20,cex=.7,
     col=ifelse(pbsc[,"all"]>0,"purple",
                ifelse(pbsc[,"interm"]>0,"red",
                       ifelse(pbsc[,"res"],"black",NA))),bty='l',
     xlab="PB z values",ylab="BB z values")
abline(h=0,v=0)

plot(pbsc[,"VB"],pbsc[,"PB"],pch=20,cex=.7,
     col=ifelse(pbsc[,"all"]>0,"purple",
                ifelse(pbsc[,"interm"]>0,"red",
                       ifelse(pbsc[,"res"],"black",NA))),bty='l',
     xlab="PB z values",ylab="VB z values")
abline(h=0,v=0)

plot(pbsc[,"BNP"],pbsc[,"PB"],pch=20,cex=.7,
     col=ifelse(pbsc[,"all"]>0,"purple",
                ifelse(pbsc[,"interm"]>0,"red",
                       ifelse(pbsc[,"res"],"black",NA))),bty='l',
     xlab="PB z values",ylab="BNP z values")
abline(h=0,v=0)

plot(pbsc[,"SJ"],pbsc[,"PB"],pch=20,cex=.7,
     col=ifelse(pbsc[,"all"]>0,"purple",
                ifelse(pbsc[,"interm"]>0,"red",
                       ifelse(pbsc[,"res"],"black",NA))),bty='l',
     xlab="PB z values",ylab="SJ z values")
abline(h=0,v=0)

plot(pbsc[,"BNP"],pbsc[,"SJ"],pch=20,cex=.7,
     col=ifelse(pbsc[,"all"]>0,"purple",
                ifelse(pbsc[,"interm"]>0,"red",
                       ifelse(pbsc[,"res"],"black",NA))),bty='l',
     xlab="SJ z values",ylab="BNP z values")
abline(h=0,v=0)

###Trying to figure out if this can be attributed to a real increase in z for resistant pops----

intermeans<-c()
for(i in 1:5){
  intermeans[i]<-mean(pbsc[interm,i+3])
}

resmeans<-c()
for(i in 1:5){
  resmeans[i]<-mean(pbsc[res,i+3])
}

allmeans<-c()
for(i in 1:5){
  allmeans[i]<-mean(pbsc[all,i+3])
}

#plotting histogram for intermediate regions----
rimeans<-c()
b<-c()
for(i in 1:5){
  for(j in 1:1000){
    b[j]<-mean(sample(pbsc[,i+3],size=388,replace=FALSE))
  }
  rimeans<-cbind(rimeans,b)
}

nam<-c("BB","VB","PB","SJ","BNP")
colnames(rimeans)<-nam
cols<-c("black","black","black","firebrick2","firebrick2")

par(mfrow=c(2,3))
for(i in 1:length(nam)){
  hist(rimeans[,i],main='',breaks=30,xlim=c(range(rimeans[,1])[[1]],intermeans[i]+.5),
       bty='l',col=cols[i],border=cols[i],xlab=nam[i])
  abline(v=intermeans[i],lwd=3,col="green")
}

#plotting histogram for resistant only regions----

rrmeans<-c()
b<-c()
for(i in 1:5){
  for(j in 1:1000){
    b[j]<-mean(sample(pbsc[,i+3],size=2435,replace=FALSE))
  }
  rrmeans<-cbind(rrmeans,b)
}

nam<-c("BB","VB","PB","SJ","BNP")
colnames(rrmeans)<-nam
cols<-c("black","black","black","firebrick2","firebrick2")

par(mfrow=c(2,3))
for(i in 1:length(nam)){
  hist(rrmeans[,i],main='',breaks=30,xlim=c(range(rrmeans[,1])[[1]],resmeans[i]+.5),
       bty='l',col=cols[i],border=cols[i],xlab=nam[i])
  abline(v=resmeans[i],lwd=3,col="green")
}

###plotting histogram for shared regions----

rameans<-c()
b<-c()
for(i in 1:5){
  for(j in 1:1000){
    b[j]<-mean(sample(pbsc[,i+3],size=252,replace=FALSE))
  }
  rameans<-cbind(rameans,b)
}

nam<-c("BB","VB","PB","SJ","BNP")
colnames(rameans)<-nam
cols<-c("black","black","black","firebrick2","firebrick2")

par(mfrow=c(2,3))
for(i in 1:length(nam)){
  hist(rameans[,i],main='',breaks=30,xlim=c(range(rameans[,1])[[1]],allmeans[i]+.5),
       bty='l',col=cols[i],border=cols[i],xlab=nam[i])
  abline(v=allmeans[i],lwd=3,col="green")
}


# palette(c("grey50","grey70"))
# par(mfrow=c(5,1),mar=c(0,3,0,0))
# plot(pbsc[,4],pch=20,cex=1.2,
#      col=ifelse((all),"purple",
#                 ifelse((res),"black",
#                        ifelse((interm),"firebrick2",
#                               ifelse((bbu),"gold2",
#                                      ifelse(pbsc[,4]>col[1],"green2",sort(as.factor(pbsc[,1]))))))),
#      xlab="",xaxt='n',ylab="BB (PBS)",cex.lab=1,cex.axis=2.2,bty="n",ylim=c(-16,23),xaxs="i",yaxs="i")
# 
# # legend("topright",legend=c("Shared by all adapted","Resistant only","Intermediate only","Shared (group non-specific)","Local"),
# #        col=c("purple","black","firebrick2","green2","gold2"),pch=20,cex=1.8,y.intersp=.5,x.intersp=.8,bty='n')
# 
# plot(pbsc[,5],pch=20,cex=1.2,
#      col=ifelse((all),"purple",
#                 ifelse((res),"black",
#                        ifelse((interm),"firebrick2",
#                               ifelse((vbu),"gold2",
#                                      ifelse(pbsc[,5]>col[2],"green2",sort(as.factor(pbsc[,1]))))))),
#      xlab="",xaxt='n',ylab="VB (PBS)",cex.lab=1,cex.axis=2.2,bty="n",ylim=c(-16,23),xaxs="i",yaxs="i")
# 
# plot(pbsc[,6],pch=20,cex=1.2,
#      col=ifelse((all),"purple",
#                 ifelse((res),"black",
#                        ifelse((interm),"firebrick2",
#                               ifelse((pbu),"gold2",
#                                      ifelse(pbsc[,6]>col[3],"green2",sort(as.factor(pbsc[,1]))))))),
#      xlab="",xaxt='n',ylab="PB (PBS)",cex.lab=1,cex.axis=2.2,bty="n",ylim=c(-16,23),xaxs="i",yaxs="i")
# 
# plot(pbsc[,7],pch=20,cex=1.2,
#      col=ifelse((all),"purple",
#                 ifelse((res),"black",
#                        ifelse((interm),"firebrick2",
#                               ifelse((sju),"gold2",
#                                      ifelse(pbsc[,7]>col[4],"green2",sort(as.factor(pbsc[,1]))))))),
#      xlab="",xaxt='n',ylab="SJ (PBS)",cex.lab=1,cex.axis=2.2,bty="n",ylim=c(-16,23),xaxs="i",yaxs="i")
# 
# plot(pbsc[,8],pch=20,cex=1.2,
#      col=ifelse((all),"purple",
#                 ifelse((res),"black",
#                        ifelse((interm),"firebrick2",
#                               ifelse((bnpu),"gold2",
#                                      ifelse(pbsc[,8]>col[5],"green2",sort(as.factor(pbsc[,1]))))))),
#      xlab="",xaxt='n',ylab="BNP (PBS)",cex.lab=1,cex.axis=2.2,bty="n",ylim=c(-16,23),xaxs="i",yaxs="i")
# 
# 
# #Regions
# 
# write.table(na.omit(pbsd_chr[all,1:3]),"~/analysis/data/dfst/zshared_outliers",row.names=FALSE,col.names=FALSE,quote=FALSE)
# write.table(na.omit(pbsd_chr[res,1:3]),"~/analysis/data/dfst/zres_outliers",row.names=FALSE,col.names=FALSE,quote=FALSE)
# write.table(na.omit(pbsd_chr[interm,1:3]),"~/analysis/data/dfst/zinterm_outliers",row.names=FALSE,col.names=FALSE,quote=FALSE)
# 
# #AHR b's
# 
# chr18<-pbs_dist[grep("chr18\\b",pbs_dist$Scaf),]
# plot(chr18[13000:17499,4],pch=20,cex=.5,col="black",ylim=c(-10,20),ylab="PBS divergence stat",xlab="tail end of chromosome 18")
# points(chr18[13000:17499,5],pch=20,cex=.5,col="grey",ylim=c(-1,1))
# points(chr18[13000:17499,6],pch=20,cex=.5,col="red",ylim=c(-1,1))
# points(chr18[13000:17499,7],pch=20,cex=.5,col="darkorange",ylim=c(-1,1))
# points(chr18[13000:17499,8],pch=20,cex=.5,col="gold",ylim=c(-1,1))
# 
# legend("topleft",y=c(.05,.1),legend=c("BB","VB","PB","SJ","BNP"),pch=20,cex=1.2,
#        col=c("black","grey","red","darkorange","gold"), x.intersp=.8,y.intersp=.8)
# 
# #arnt
# 
# chr8<-pbs_dist[grep("chr8\\b",pbs_dist$Scaf),]
# plot(chr8[13000:17499,4],pch=20,cex=.5,col="black",ylim=c(-10,20),ylab="PBS divergence stat",xlab="tail end of chromosome 8")
# points(chr8[13000:17499,5],pch=20,cex=.5,col="grey",ylim=c(-1,1))
# points(chr8[13000:17499,6],pch=20,cex=.5,col="red",ylim=c(-1,1))
# points(chr8[13000:17499,7],pch=20,cex=.5,col="darkorange",ylim=c(-1,1))
# points(chr8[13000:17499,8],pch=20,cex=.5,col="gold",ylim=c(-1,1))
# 
# legend("topleft",y=c(.05,.1),legend=c("BB","VB","PB","SJ","BNP"),pch=20,cex=1.2,
#        col=c("black","grey","red","darkorange","gold"), x.intersp=.8,y.intersp=.8)