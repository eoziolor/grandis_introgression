# install.packages('stringr')
# install.packages('dplyr')
# install.packages('gtools')
#install.packages('magrittr')
#install.packages("XML")
library(XML)
library(magrittr)
library(stringr)
library(dplyr)
library(gtools)

#Reading in table and getting quantiles----
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
#plot(pbs[subw,4],pch=20,cex=.5,col=factor(pbs[subw,1]))
#plot(pbs[subw,4],pch=20,cex=.5,col=chr[pbs[subw,1],2])
#legend('topright',legend=levels(mixedsort(pbs[,1])),col=1:2,cex=.5,pch=1)

#Subsampling into chromosome only------
pbsct<-pbs %>% filter(str_detect(Scaf,"chr"))

nsnps <-pbsct[,"keep"]
subwc<-nsnps>0


pbsctt<-pbsct[subwc,]
rownames(pbsctt)<-seq(1:dim(pbsctt[subwc,])[1])
pbsctt<-pbsctt[,1:8]

#removing the crappy scaffold that has the first 141 windows of chromosome 16 (mismapped from chr1); this is discovered in dxy.r script in introgression folder
chr16<-str_detect(pbsctt$Scaf,"chr16") #grab chr16
ord<-order(pbsctt[chr16,"BBpbs"],decreasing=TRUE)
pbsc16<-pbsctt[chr16,]
pbsc16[ord,]
chr16rows<-as.numeric(rownames(pbsctt[chr16,])) #grab rownames for it
crappyrows<-chr16rows[1:300] #get the first 141 rows which contain scaffold "crappy"
pbsc<-pbsctt[-c(crappyrows),] #remove thos rows from total
chr16.2<-str_detect(pbsc$Scaf,"chr16")
head(pbsc[chr16,])

###Doing this on merged windows to avoid patchyness of peak coloration------

pbs_out_temp<-read.table("~/analysis/data/fst/PBSoutliers_5kb_all_max.bed",stringsAsFactors = FALSE) #loads a pbs vector with windows merged within 50kb of each other and with max and windows count statistics
names<-c("Scaf","start","end","BBmax","BBcount","VBmax","VBcount","PBmax","PBcount","SJmax","SJcount","BNPmax","BNPcount")
colnames(pbs_out_temp)<-names

pbs_out<-pbs_out_temp %>% filter(str_detect(Scaf,"chr")) #selecting only chromosome mapped scaffolds

#checking for whether those are outliers in different groups--------
all<-pbs_out[,4]>col[1] & pbs_out[,6]>col[2] & pbs_out[,8]>col[3] & pbs_out[,10]>col[4] & pbs_out[,12]>col[5]
res<-pbs_out[,4]>col[1] & pbs_out[,6]>col[2] & pbs_out[,8]>col[3] & pbs_out[,10]<col[4] & pbs_out[,12]<col[5]
interm<-pbs_out[,4]<col[1] & pbs_out[,6]<col[2] & pbs_out[,8]<col[3] & pbs_out[,10]>col[4] & pbs_out[,12]>col[5]
bbu<-pbs_out[,4]>col[1] & pbs_out[,6]<col[2] & pbs_out[,8]<col[3] & pbs_out[,10]<col[4] & pbs_out[,12]<col[5]
vbu<-pbs_out[,4]<col[1] & pbs_out[,6]>col[2] & pbs_out[,8]<col[3] & pbs_out[,10]<col[4] & pbs_out[,12]<col[5]
pbu<-pbs_out[,4]<col[1] & pbs_out[,6]<col[2] & pbs_out[,8]>col[3] & pbs_out[,10]<col[4] & pbs_out[,12]<col[5]
sju<-pbs_out[,4]<col[1] & pbs_out[,6]<col[2] & pbs_out[,8]<col[3] & pbs_out[,10]>col[4] & pbs_out[,12]<col[5]
bnpu<-pbs_out[,4]<col[1] & pbs_out[,6]<col[2] & pbs_out[,8]<col[3] & pbs_out[,10]<col[4] & pbs_out[,12]>col[5]


write.table(pbsc[,1:3],"~/analysis/data/fst/PBS_keep_5kb.bed",row.names = FALSE,col.names = FALSE,quote=FALSE)
write.table(pbs_out[all,1:3],"~/analysis/data/fst/pbs_regions_sharedall.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(pbs_out[res,1:3],"~/analysis/data/fst/pbs_regions_sharedres.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(pbs_out[interm,1:3],"~/analysis/data/fst/pbs_regions_sharedinterm.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(pbs_out[bbu,1:3],"~/analysis/data/fst/pbs_regions_sharedbbu.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(pbs_out[vbu,1:3],"~/analysis/data/fst/pbs_regions_sharedvbu.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(pbs_out[pbu,1:3],"~/analysis/data/fst/pbs_regions_sharedpbu.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(pbs_out[sju,1:3],"~/analysis/data/fst/pbs_regions_sharedsju.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(pbs_out[bnpu,1:3],"~/analysis/data/fst/pbs_regions_sharedbnpu.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)

#source("http://bioconductor.org/biocLite.R")
#biocLite()
#biocLite('rtracklayer')

#Finding the overlaps in full data---------
library("rtracklayer")

bed1=import("~/analysis/data/fst/PBS_keep_5kb.bed")

bedall=import("~/analysis/data/fst/pbs_regions_sharedall.bed")
bed1overlall=bed1[bed1 %over% bedall]
hitsall<-findOverlaps(bedall,bed1)
allhit<-subjectHits(hitsall)

bedres=import("~/analysis/data/fst/pbs_regions_sharedres.bed")
bed1overlres=bed1[bed1 %over% bedres]
hitsres<-findOverlaps(bedres,bed1)
reshit<-subjectHits(hitsres)

bedinterm=import("~/analysis/data/fst/pbs_regions_sharedinterm.bed")
bed1overlinterm=bed1[bed1 %over% bedinterm]
hitsinterm<-findOverlaps(bedinterm,bed1)
intermhit<-subjectHits(hitsinterm)

bedbbu=import("~/analysis/data/fst/pbs_regions_sharedbbu.bed")
bed1overlbbu=bed1[bed1 %over% bedbbu]
hitsbbu<-findOverlaps(bedbbu,bed1)
bbuhit<-subjectHits(hitsbbu)

bedvbu=import("~/analysis/data/fst/pbs_regions_sharedvbu.bed")
bed1overlvbu=bed1[bed1 %over% bedvbu]
hitsvbu<-findOverlaps(bedvbu,bed1)
vbuhit<-subjectHits(hitsvbu)

bedpbu=import("~/analysis/data/fst/pbs_regions_sharedpbu.bed")
bed1overlpbu=bed1[bed1 %over% bedpbu]
hitspbu<-findOverlaps(bedpbu,bed1)
pbuhit<-subjectHits(hitspbu)

bedsju=import("~/analysis/data/fst/pbs_regions_sharedsju.bed")
bed1overlsju=bed1[bed1 %over% bedsju]
hitssju<-findOverlaps(bedsju,bed1)
sjuhit<-subjectHits(hitssju)

bedbnpu=import("~/analysis/data/fst/pbs_regions_sharedbnpu.bed")
bed1overlbnpu=bed1[bed1 %over% bedbnpu]
hitsbnpu<-findOverlaps(bedbnpu,bed1)
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


#plotting those results by using the pbs_out vector-------------
palette(c("grey40","grey60","grey80","grey20"))
par(mfrow=c(5,1),mar=c(0,3,0,0))
plot(pbsc[,4],pch=20,cex=1.2,
     col=ifelse(pbsc[,"all"]>0,"purple",
                ifelse(pbsc[,"res"]>0,"black",
                       ifelse(pbsc[,"interm"]>0,"firebrick2",
                              ifelse(pbsc[,"bbu"]>0,"gold2",
                                     ifelse(pbsc[,4]>col[1],"green2",as.factor(pbsc[,1])))))),
     xlab="",xaxt='n',cex.lab=1,cex.axis=2.2,bty="n",ylim=c(-.5,3.8),yaxs="i")

plot(pbsc[,5],pch=20,cex=1.2,
     col=ifelse(pbsc[,"all"]>0,"purple",
                ifelse(pbsc[,"res"]>0,"black",
                       ifelse(pbsc[,"interm"]>0,"firebrick2",
                              ifelse(pbsc[,"vbu"]>0,"gold2",
                                     ifelse(pbsc[,5]>col[2],"green2",as.factor(pbsc[,1])))))),
     xlab="",xaxt='n',cex.lab=1,cex.axis=2.2,bty="n",ylim=c(-.5,3.8),yaxs="i")

plot(pbsc[,6],pch=20,cex=1.2,
     col=ifelse(pbsc[,"all"]>0,"purple",
                ifelse(pbsc[,"res"]>0,"black",
                       ifelse(pbsc[,"interm"]>0,"firebrick2",
                              ifelse(pbsc[,"pbu"]>0,"gold2",
                                     ifelse(pbsc[,6]>col[3],"green2",as.factor(pbsc[,1])))))),
     xlab="",xaxt='n',cex.lab=1,cex.axis=2.2,bty="n",ylim=c(-.5,3.8),yaxs="i")

plot(pbsc[,7],pch=20,cex=1.2,
     col=ifelse(pbsc[,"all"]>0,"purple",
                ifelse(pbsc[,"res"]>0,"black",
                       ifelse(pbsc[,"interm"]>0,"firebrick2",
                              ifelse(pbsc[,"sju"]>0,"gold2",
                                     ifelse(pbsc[,7]>col[4],"green2",as.factor(pbsc[,1])))))),
     xlab="",xaxt='n',cex.lab=1,cex.axis=2.2,bty="n",ylim=c(-.5,3.8),yaxs="i")

plot(pbsc[,8],pch=20,cex=1.2,
     col=ifelse(pbsc[,"all"]>0,"purple",
                ifelse(pbsc[,"res"]>0,"black",
                       ifelse(pbsc[,"interm"]>0,"firebrick2",
                              ifelse(pbsc[,"bnpu"]>0,"gold2",
                                     ifelse(pbsc[,8]>col[5],"green2",as.factor(pbsc[,1])))))),
     xlab="",xaxt='n',cex.lab=1,cex.axis=2.2,bty="n",ylim=c(-.5,3.8),yaxs="i")


###Plotting outliers vs each other-----------
pbsct<-pbs %>% filter(str_detect(Scaf,"chr"))

nsnps <-pbsct[,"keep"]
subwc<-nsnps>0


pbsc<-pbsct[subwc,]
rownames(pbsc)<-seq(1:dim(pbsc[subwc,])[1])
pbsc<-pbsc[,1:8]

#checking for whether those are outliers in different groups
all<-pbsc[,4]>col[1] & pbsc[,5]>col[2] & pbsc[,6]>col[3] & pbsc[,7]>col[4] & pbsc[,8]>col[5]
res<-pbsc[,4]>col[1] & pbsc[,5]>col[2] & pbsc[,6]>col[3] & pbsc[,7]<col[4] & pbsc[,8]<col[5]
interm<-pbsc[,4]<col[1] & pbsc[,5]<col[2] & pbsc[,6]<col[3] & pbsc[,7]>col[4] & pbsc[,8]>col[5]
bbu<-pbsc[,4]>col[1] & pbsc[,5]<col[2] & pbsc[,6]<col[3] & pbsc[,7]<col[4] & pbsc[,8]<col[5]
vbu<-pbsc[,4]<col[1] & pbsc[,5]>col[2] & pbsc[,6]<col[3] & pbsc[,7]<col[4] & pbsc[,8]<col[5]
pbu<-pbsc[,4]<col[1] & pbsc[,5]<col[2] & pbsc[,6]>col[3] & pbsc[,7]<col[4] & pbsc[,8]<col[5]
sju<-pbsc[,4]<col[1] & pbsc[,5]<col[2] & pbsc[,6]<col[3] & pbsc[,7]>col[4] & pbsc[,8]<col[5]
bnpu<-pbsc[,4]<col[1] & pbsc[,5]<col[2] & pbsc[,6]<col[3] & pbsc[,7]<col[4] & pbsc[,8]>col[5]

write.table(pbsc[,1:3],"~/analysis/data/fst/PBS_keep_5kb.bed",row.names = FALSE,col.names = FALSE,quote=FALSE)
write.table(na.omit(pbsc[all,1:3]),"~/analysis/data/fst/pbs_regions_sharedall.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(na.omit(pbsc[res,1:3]),"~/analysis/data/fst/pbs_regions_sharedres.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(na.omit(pbsc[interm,1:3]),"~/analysis/data/fst/pbs_regions_sharedinterm.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(na.omit(pbsc[bbu,1:3]),"~/analysis/data/fst/pbs_regions_sharedbbu.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(na.omit(pbsc[vbu,1:3]),"~/analysis/data/fst/pbs_regions_sharedvbu.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(na.omit(pbsc[pbu,1:3]),"~/analysis/data/fst/pbs_regions_sharedpbu.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(na.omit(pbsc[sju,1:3]),"~/analysis/data/fst/pbs_regions_sharedsju.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(na.omit(pbsc[bnpu,1:3]),"~/analysis/data/fst/pbs_regions_sharedbnpu.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)

#source("http://bioconductor.org/biocLite.R")
#biocLite()
library("rtracklayer")

bed1=import("~/analysis/data/fst/PBS_keep_5kb.bed")

bedall=import("~/analysis/data/fst/pbs_regions_sharedall.bed")
bed1overlall=bed1[bed1 %over% bedall]
hitsall<-findOverlaps(bedall,bed1,type="equal")
allhit<-subjectHits(hitsall)

bedres=import("~/analysis/data/fst/pbs_regions_sharedres.bed")
bed1overlres=bed1[bed1 %over% bedres]
hitsres<-findOverlaps(bedres,bed1,type="equal")
reshit<-subjectHits(hitsres)

bedinterm=import("~/analysis/data/fst/pbs_regions_sharedinterm.bed")
bed1overlinterm=bed1[bed1 %over% bedinterm]
hitsinterm<-findOverlaps(bedinterm,bed1,type="equal")
intermhit<-subjectHits(hitsinterm)

bedbbu=import("~/analysis/data/fst/pbs_regions_sharedbbu.bed")
bed1overlbbu=bed1[bed1 %over% bedbbu]
hitsbbu<-findOverlaps(bedbbu,bed1,type="equal")
bbuhit<-subjectHits(hitsbbu)

bedvbu=import("~/analysis/data/fst/pbs_regions_sharedvbu.bed")
bed1overlvbu=bed1[bed1 %over% bedvbu]
hitsvbu<-findOverlaps(bedvbu,bed1,type="equal")
vbuhit<-subjectHits(hitsvbu)

bedpbu=import("~/analysis/data/fst/pbs_regions_sharedpbu.bed")
bed1overlpbu=bed1[bed1 %over% bedpbu]
hitspbu<-findOverlaps(bedpbu,bed1,type="equal")
pbuhit<-subjectHits(hitspbu)

bedsju=import("~/analysis/data/fst/pbs_regions_sharedsju.bed")
bed1overlsju=bed1[bed1 %over% bedsju]
hitssju<-findOverlaps(bedsju,bed1,type="equal")
sjuhit<-subjectHits(hitssju)

bedbnpu=import("~/analysis/data/fst/pbs_regions_sharedbnpu.bed")
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

plot(pbsc[,"SJ"],pbsc[,"PB"],pch=20,cex=.7,
     col=ifelse(pbsc[,"all"]>0,"purple",
                ifelse(pbsc[,"interm"]>0,"red",
                       ifelse(pbsc[,"res"],"black",NA))),bty='l',
     xlab="PB z values",ylab="SJ z values")
abline(h=0,v=0)

plot(pbsc[,"BNP"],pbsc[,"PB"],pch=20,cex=.7,
     col=ifelse(pbsc[,"all"]>0,"purple",
                ifelse(pbsc[,"interm"]>0,"red",
                       ifelse(pbsc[,"res"],"black",NA))),bty='l',
     xlab="PB z values",ylab="BNP z values")
abline(h=0,v=0)

plot(pbsc[,"BNP"],pbsc[,"SJ"],pch=20,cex=.7,
     col=ifelse(pbsc[,"all"]>0,"purple",
                ifelse(pbsc[,"interm"]>0,"red",
                       ifelse(pbsc[,"res"],"black",NA))),bty='l',
     xlab="SJ z values",ylab="BNP z values")
abline(h=0,v=0)

#Plotting outliers in ggplot----
require(gridExtra)
pbsw<-cbind(pbsc[,1:8],0)
newn<-c("Scaf","start","end","BB","VB","PB","SJ","BNP","outliers")
colnames(pbsw)<-newn

pbsw[allhit,"outliers"]<-"all"
pbsw[reshit,"outliers"]<-"res"
pbsw[intermhit,"outliers"]<-"int"

pbsws<-subset(pbsw,outliers!=0)
a<-ggplot(pbsws,
       aes(BB,VB,col=outliers))+
  geom_point() +
  scale_color_manual(values=c("purple","firebrick2","black")) +
  theme_classic()

b<-ggplot(pbsws,
          aes(BB,PB,col=outliers))+
  geom_point() +
  scale_color_manual(values=c("purple","firebrick2","black")) +
  theme_classic()

c<-ggplot(pbsws,
       aes(VB,PB,col=outliers))+
  geom_point() +
  scale_color_manual(values=c("purple","firebrick2","black")) +
  theme_classic()

d<-ggplot(pbsws,
          aes(SJ,BB,col=outliers))+
  geom_point() +
  scale_color_manual(values=c("purple","firebrick2","black")) +
  theme_classic()

e<-ggplot(pbsws,
          aes(BNP,BB,col=outliers))+
  geom_point() +
  scale_color_manual(values=c("purple","firebrick2","black")) +
  theme_classic()

f<-ggplot(pbsws,
          aes(BNP,SJ,col=outliers))+
  geom_point() +
  scale_color_manual(values=c("purple","firebrick2","black")) +
  theme_classic()

grid.arrange(a,b,c,d,e,f,ncol=3,nrow=2)

# ggplot(pbsw,
#        aes(SJ,BB,col=outliers))+
#   geom_point() +
#   scale_color_manual(values=c("grey","purple","firebrick2","black")) +
#   theme_classic()

###Trying to figure out if this can be attributed to a real increase in z for resistant pops----

intermeans<-c()
for(i in 1:5){
  intermeans[i]<-mean(pbsc[interm,i+3],na.rm=TRUE)
}

resmeans<-c()
for(i in 1:5){
  resmeans[i]<-mean(pbsc[res,i+3],na.rm=TRUE)
}

allmeans<-c()
for(i in 1:5){
  allmeans[i]<-mean(pbsc[all,i+3],na.rm=TRUE)
}

#plotting histogram for intermediate regions----
rimeans<-c()
b<-c()
for(i in 1:5){
  for(j in 1:1000){
    b[j]<-mean(sample(pbsc[,i+3],size=124,replace=FALSE))
  }
  rimeans<-cbind(rimeans,b)
}

nam<-c("BB","VB","PB","SJ","BNP")
colnames(rimeans)<-nam
cols<-c("black","black","black","firebrick2","firebrick2")

par(mfrow=c(2,3),mar=c(4,4,2,2))
for(i in 1:length(nam)){
  hist(rimeans[,i],main='',breaks=30,xlim=c(range(rimeans[,i]-.1,na.rm=TRUE)[[1]],intermeans[i]+.5),
       bty='l',col=cols[i],border=cols[i],xlab=nam[i],cex.axis=3,ylab='')
  abline(v=intermeans[i],lwd=3,col="green")
  box(bty='l',lwd=3)
}

#plotting histogram for resistant only regions----

rrmeans<-c()
b<-c()
for(i in 1:5){
  for(j in 1:1000){
    b[j]<-mean(sample(pbsc[,i+3],size=2549,replace=FALSE))
  }
  rrmeans<-cbind(rrmeans,b)
}

nam<-c("BB","VB","PB","SJ","BNP")
colnames(rrmeans)<-nam
cols<-c("black","black","black","firebrick2","firebrick2")

par(mfrow=c(2,3))
for(i in 1:length(nam)){
  hist(rrmeans[,i],main='',breaks=30,xlim=c(range(rrmeans[,i]-.1,na.rm=TRUE)[[1]],resmeans[i]+.5),
       bty='l',col=cols[i],border=cols[i],xlab=nam[i])
  abline(v=resmeans[i],lwd=3,col="green")
}

###plotting histogram for shared regions----

rameans<-c()
b<-c()
for(i in 1:5){
  for(j in 1:1000){
    b[j]<-mean(sample(pbsc[,i+3],size=259,replace=FALSE))
  }
  rameans<-cbind(rameans,b)
}

nam<-c("BB","VB","PB","SJ","BNP")
colnames(rameans)<-nam
cols<-c("black","black","black","firebrick2","firebrick2")

par(mfrow=c(2,3))
for(i in 1:length(nam)){
  hist(rameans[,i],main='',breaks=30,xlim=c(range(rameans[,i]-.1,na.rm=TRUE)[[1]],allmeans[i]+.5),
       bty='l',col=cols[i],border=cols[i],xlab=nam[i])
  abline(v=allmeans[i],lwd=3,col="green")
}

#smoothing funciton----
subsmooth <- function(vec,by=10,width=11){
  
  len <- length(vec)
  subl <- seq(from=by,to=len,by=by)
  submax <- length(subl)
  width <- width/2
  test <- vec[subl]
  
  for(i in 1:submax){
    
    j <- i - width
    k <- i + width
    if(j < 1) {j <- 1}
    if(k > submax) {k <- submax}
    test[i] <- mean(test[j:k],na.rm=TRUE)
  }
  
  return(test)
  
}


###Plotting CHR1/AHR region------------
pbsc1<-pbsc %>% filter(str_detect(Scaf,"\\bchr1\\b"))

par(mfrow=c(5,1),mar=c(3,3,0,0),mgp=c(1,1,0))

for(i in 1:5){
  plot(subsmooth(pbsc1[1:5000,i+3]),pch=20,cex=.5,ylim=c(0,1.8),bty='l',cex.axis=2,ylab='',xlab='')
  box(bty='l',lwd=3)
  abline(v=c(50,60),col="red",lty=2,lwd=1.5)
}

#plotting ARNT chr8----
pbsc8<-pbsc %>% filter(str_detect(Scaf,"\\bchr8\\b"))

par(mfrow=c(5,1),mar=c(3,3,0,0),mgp=c(1,1,0))

for(i in 1:5){
  plot(subsmooth(pbsc8[14000:18000,i+3]),pch=20,cex=.5,ylim=c(0,.5),bty='l',cex.axis=2,ylab='',xlab='')
  box(bty='l',lwd=3)
  abline(v=c(170,173),col="red",lty=2,lwd=1.5)
}


#plotting AIP chr2----
pbsc2<-pbsc %>% filter(str_detect(Scaf,"\\bchr2\\b"))

par(mfrow=c(5,1),mar=c(3,3,0,0),mgp=c(1,1,0))

for(i in 1:5){
  plot(subsmooth(pbsc2[23000:26000,i+3]),pch=20,cex=.5,ylim=c(0,.4),bty='l',cex.axis=2,ylab='',xlab='')
  box(bty='l',lwd=3)
  abline(v=c(140,142),col="red",lty=2,lwd=1.5)
}

#plotting AQP3

pbsc24<-pbsc %>% filter(str_detect(Scaf,"\\bchr24\\b"))

par(mfrow=c(5,1),mar=c(3,3,0,0),mgp=c(1,1,0))

for(i in 1:5){
  plot(subsmooth(pbsc24[22000:26000,i+3]),pch=20,cex=.5,ylim=c(0,.5),bty='l',cex.axis=2,ylab='',xlab='')
  box(bty='l',lwd=3)
  abline(v=c(159,162),col="red",lty=2,lwd=1.5)
}

#plotting AQP3

pbsc11<-pbsc %>% filter(str_detect(Scaf,"\\bchr11\\b"))

par(mfrow=c(5,1),mar=c(3,3,0,0),mgp=c(1,1,0))

for(i in 1:5){
  plot(subsmooth(pbsc11[27000:28940,i+3]),pch=20,cex=.5,ylim=c(0,.5),bty='l',cex.axis=2,ylab='',xlab='')
  box(bty='l',lwd=3)
  abline(v=c(181,188),col="red",lty=2,lwd=1.5)
}


# #old chr1 plotting----
# 
# plot(pbsc1[,4],pch=20,cex=1.2,
#      col=ifelse(pbsc1[,"all"]>0,"purple",
#                 ifelse(pbsc1[,"res"]>0,"black",
#                        ifelse(pbsc1[,"interm"]>0,"firebrick2",
#                               ifelse(pbsc1[,"bbu"]>0,"gold2",
#                                      ifelse(pbsc1[,4]>col[1],"green2",sort(as.factor(pbsc1[,1]))))))),
#      xlab="",xaxt='n',cex.lab=1,cex.axis=2.2,bty="n",ylim=c(-.5,3.8),yaxs="i")
# abline(v=c(547,555),lty=2,col="red")
# 
# plot(pbsc1[,5],pch=20,cex=1.2,
#      col=ifelse(pbsc1[,"all"]>0,"purple",
#                 ifelse(pbsc1[,"res"]>0,"black",
#                        ifelse(pbsc1[,"interm"]>0,"firebrick2",
#                               ifelse(pbsc1[,"vbu"]>0,"gold2",
#                                      ifelse(pbsc1[,5]>col[2],"green2",sort(as.factor(pbsc1[,1]))))))),
#      xlab="",xaxt='n',cex.lab=1,cex.axis=2.2,bty="n",ylim=c(-.5,3.8),yaxs="i")
# abline(v=c(547,555),lty=2,col="red")
# 
# plot(pbsc1[,6],pch=20,cex=1.2,
#      col=ifelse(pbsc1[,"all"]>0,"purple",
#                 ifelse(pbsc1[,"res"]>0,"black",
#                        ifelse(pbsc1[,"interm"]>0,"firebrick2",
#                               ifelse(pbsc1[,"pbu"]>0,"gold2",
#                                      ifelse(pbsc1[,6]>col[3],"green2",sort(as.factor(pbsc1[,1]))))))),
#      xlab="",xaxt='n',cex.lab=1,cex.axis=2.2,bty="n",ylim=c(-.5,3.8),yaxs="i")
# abline(v=c(547,555),lty=2,col="red")
# 
# plot(pbsc1[,7],pch=20,cex=1.2,
#      col=ifelse(pbsc1[,"all"]>0,"purple",
#                 ifelse(pbsc1[,"res"]>0,"black",
#                        ifelse(pbsc1[,"interm"]>0,"firebrick2",
#                               ifelse(pbsc1[,"sju"]>0,"gold2",
#                                      ifelse(pbsc1[,7]>col[4],"green2",sort(as.factor(pbsc1[,1]))))))),
#      xlab="",xaxt='n',cex.lab=1,cex.axis=2.2,bty="n",ylim=c(-.5,3.8),yaxs="i")
# abline(v=c(547,555),lty=2,col="red")
# 
# plot(pbsc1[,8],pch=20,cex=1.2,
#      col=ifelse(pbsc1[,"all"]>0,"purple",
#                 ifelse(pbsc1[,"res"]>0,"black",
#                        ifelse(pbsc1[,"interm"]>0,"firebrick2",
#                               ifelse(pbsc1[,"bnpu"]>0,"gold2",
#                                      ifelse(pbsc1[,8]>col[5],"green2",sort(as.factor(pbsc1[,1]))))))),
#      xlab="",xaxt='n',cex.lab=1,cex.axis=2.2,bty="n",ylim=c(-.5,3.8),yaxs="i")
# abline(v=c(547,555),lty=2,col="red")
# 

###plotting in 5kb windows
#####Plotting common regions
# 
# all<-pbsc[,4]>col[1] & pbsc[,5]>col[2] & pbsc[,6]>col[3] & pbsc[,7]>col[4] & pbsc[,8]>col[5]
# res<-pbsc[,4]>col[1] & pbsc[,5]>col[2] & pbsc[,6]>col[3] & pbsc[,7]<col[4] & pbsc[,8]<col[5]
# interm<-pbsc[,4]<col[1] & pbsc[,5]<col[2] & pbsc[,6]<col[3] & pbsc[,7]>col[4] & pbsc[,8]>col[5]
# bbu<-pbsc[,4]>col[1] & pbsc[,5]<col[2] & pbsc[,6]<col[3] & pbsc[,7]<col[4] & pbsc[,8]<col[5]
# vbu<-pbsc[,4]<col[1] & pbsc[,5]>col[2] & pbsc[,6]<col[3] & pbsc[,7]<col[4] & pbsc[,8]<col[5]
# pbu<-pbsc[,4]<col[1] & pbsc[,5]<col[2] & pbsc[,6]>col[3] & pbsc[,7]<col[4] & pbsc[,8]<col[5]
# sju<-pbsc[,4]<col[1] & pbsc[,5]<col[2] & pbsc[,6]<col[3] & pbsc[,7]>col[4] & pbsc[,8]<col[5]
# bnpu<-pbsc[,4]<col[1] & pbsc[,5]<col[2] & pbsc[,6]<col[3] & pbsc[,7]<col[4] & pbsc[,8]>col[5]
# 
# 
# palette(c("grey50","grey70"))
# par(mfrow=c(5,1),mar=c(0,3,0,0))
# plot(pbsc[,4],pch=20,cex=1.2,
#      col=ifelse((all),"purple",
#                 ifelse((res),"black",
#                        ifelse((interm),"firebrick2",
#                               ifelse((bbu),"gold2",
#                                      ifelse(pbsc[,4]>col[1],"green2",sort(as.factor(pbsc[,1]))))))),
#      xlab="",xaxt='n',ylab="BB (PBS)",cex.lab=1,cex.axis=2.2,bty="n",ylim=c(-.5,3.8),xaxs="i",yaxs="i")
# 
# # legend("topright",legend=c("Shared by all adapted","Resistant only","Intermediate only","Shared (group non-specific)","Local"),
# #        col=c("purple","black","firebrick2","green2","gold2"),pch=20,cex=1.8,y.intersp=.5,x.intersp=.8,bty='n')
# 
# plot(pbsc[,5],pch=20,cex=1.2,
#      col=ifelse((all),"purple",
#                 ifelse((res),"black",
#                        ifelse((interm),"firebrickas.numeric(rownames(pbsc))==allhit2",
#                               ifelse((vbu),"gold2",
#                                      ifelse(pbsc[,5]>col[2],"green2",sort(as.factor(pbsc[,1]))))))),
#      xlab="",xaxt='n',ylab="VB (PBS)",cex.lab=1,cex.axis=2.2,bty="n",ylim=c(-.5,3.8),xaxs="i",yaxs="i")
# 
# plot(pbsc[,6],pch=20,cex=1.2,
#      col=ifelse((all),"purple",
#                 ifelse((res),"black",
#                        ifelse((interm),"firebrick2",
#                               ifelse((pbu),"gold2",
#                                      ifelse(pbsc[,6]>col[3],"green2",sort(as.factor(pbsc[,1]))))))),
#      xlab="",xaxt='n',ylab="PB (PBS)",cex.lab=1,cex.axis=2.2,bty="n",ylim=c(-.5,3.8),xaxs="i",yaxs="i")
# 
# plot(pbsc[,7],pch=20,cex=1.2,
#      col=ifelse((all),"purple",
#                 ifelse((res),"black",
#                        ifelse((interm),"firebrick2",
#                               ifelse((sju),"gold2",
#                                      ifelse(pbsc[,7]>col[4],"green2",sort(as.factor(pbsc[,1]))))))),
#      xlab="",xaxt='n',ylab="SJ (PBS)",cex.lab=1,cex.axis=2.2,bty="n",ylim=c(-.5,3.8),xaxs="i",yaxs="i")
# 
# plot(pbsc[,8],pch=20,cex=1.2,
#      col=ifelse((all),"purple",
#                 ifelse((res),"black",
#                        ifelse((interm),"firebrick2",
#                               ifelse((bnpu),"gold2",
#                                      ifelse(pbsc[,8]>col[5],"green2",sort(as.factor(pbsc[,1]))))))),
#      xlab="",xaxt='n',ylab="BNP (PBS)",cex.lab=1,cex.axis=2.2,bty="n",ylim=c(-.5,3.8),xaxs="i",yaxs="i")


###chr8
# 
# chr2<-pbsc[grep("chr2\\b",pbsc$Scaf),]
# plot(chr2[,4],pch=20,cex=.5)
# 
# chr24<-pbsc[grep("chr24\\b",pbsc$Scaf),]
# plot(chr24[,4],pch=20,cex=.5)

