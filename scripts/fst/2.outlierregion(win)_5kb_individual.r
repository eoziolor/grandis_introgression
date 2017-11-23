#Open fsttopbs first and load the PBS files

pbs<-read.table("~/analysis/data/fst/allpbs5kb",header=FALSE,stringsAsFactors = FALSE)
pbsname<-c("Scaf","start","end","BBpbs","VBpbs","PBpbs","SJpbs","BNPpbs","keep")
colnames(pbs)<-pbsname

quantile(pbs[,4],probs=.99,na.rm=TRUE)
quantile(pbs[,5],probs=.99,na.rm=TRUE)
quantile(pbs[,6],probs=.99,na.rm=TRUE)
quantile(pbs[,7],probs=.99,na.rm=TRUE)
quantile(pbs[,8],probs=.99,na.rm=TRUE)


Quantile info
BBpbs - 0.2533448
VBpbs - 0.1461146 
PBpbs - 0.1765786 
SJpbs - 0.2164984
BNP - 0.05471878

###for each individual pop, take the 1% outliers and save them in separate files

##BB

cat ~/analysis/data/fst/allpbs5kb | awk '$9>0' | \
grep -v NA | \
awk '$4>0.2533448' | \
~/program/bedtools2/bin/bedtools merge -i stdin -d 5000 \
-c 4,4 \
-o sum,count \
-g <(cut -f 1-2 ~/analysis/data/genome/unsplit_merge.fasta.fai) > ~/analysis/data/fst/individual_pbs/BB_pbs_5kb.bed

#VB

cat ~/analysis/data/fst/allpbs5kb | awk '$9>0' | \
grep -v NA | \
awk '$5>0.1461146' | \
~/program/bedtools2/bin/bedtools merge -i stdin -d 5000 \
-c 5,5 \
-o sum,count \
-g <(cut -f 1-2 ~/analysis/data/genome/unsplit_merge.fasta.fai) > ~/analysis/data/fst/individual_pbs/VB_pbs_5kb.bed

#PB

cat ~/analysis/data/fst/allpbs5kb | awk '$9>0' | \
grep -v NA | \
awk '$6>0.1765786' | \
~/program/bedtools2/bin/bedtools merge -i stdin -d 5000 \
-c 6,6 \
-o sum,count \
-g <(cut -f 1-2 ~/analysis/data/genome/unsplit_merge.fasta.fai) > ~/analysis/data/fst/individual_pbs/PB_pbs_5kb.bed

#SJ

cat ~/analysis/data/fst/allpbs5kb | awk '$9>0' | \
grep -v NA | \
awk '$7>0.2164984' | \
~/program/bedtools2/bin/bedtools merge -i stdin -d 5000 \
-c 7,7 \
-o sum,count \
-g <(cut -f 1-2 ~/analysis/data/genome/unsplit_merge.fasta.fai) > ~/analysis/data/fst/individual_pbs/SJ_pbs_5kb.bed

#BNP

cat ~/analysis/data/fst/allpbs5kb | awk '$9>0' | \
grep -v NA | \
awk '$8>0.05471878' | \
~/program/bedtools2/bin/bedtools merge -i stdin -d 5000 \
-c 8,8 \
-o sum,count \
-g <(cut -f 1-2 ~/analysis/data/genome/unsplit_merge.fasta.fai) > ~/analysis/data/fst/individual_pbs/BNP_pbs_5kb.bed

#Put it back into R to look at it

bb<-read.table("~/analysis/data/fst/individual_pbs/BB_pbs_5kb.bed",stringsAsFactors = FALSE)
vb<-read.table("~/analysis/data/fst/individual_pbs/VB_pbs_5kb.bed",stringsAsFactors = FALSE)
pb<-read.table("~/analysis/data/fst/individual_pbs/PB_pbs_5kb.bed",stringsAsFactors = FALSE)
sj<-read.table("~/analysis/data/fst/individual_pbs/SJ_pbs_5kb.bed",stringsAsFactors = FALSE)
bnp<-read.table("~/analysis/data/fst/individual_pbs/BNP_pbs_5kb.bed",stringsAsFactors = FALSE)

colnames(bb)<- c("Scaf","start","end","sum", "count")
colnames(vb)<- c("Scaf","start","end","sum", "count")
colnames(pb)<- c("Scaf","start","end","sum", "count")
colnames(sj)<- c("Scaf","start","end","sum", "count")
colnames(bnp)<- c("Scaf","start","end","sum", "count")

bbdist<-bb[,3]-bb[,2]
vbdist<-vb[,3]-vb[,2]
pbdist<-pb[,3]-pb[,2]
sjdist<-sj[,3]-sj[,2]
bnpdist<-bnp[,3]-bnp[,2]

###plotting distribution of 1% blocks
bbd<-density(bbdist,na.rm=TRUE)
vbd<-density(vbdist,na.rm=TRUE)
pbd<-density(pbdist,na.rm=TRUE)
sjd<-density(sjdist,na.rm=TRUE)
bnpd<-density(bnpdist,na.rm=TRUE)

plot(bbd,xlim=c(0,30000),ylim=c(0,0.0004),
     col="black",bty="l",cex.lab=2,cex.lab=2,xlab="Size of region with high diversity",main='',lwd=3)
lines(vbd,col="grey",lwd=3)
lines(pbd,col="red",lwd=3)
lines(sjd,col="darkorange",lwd=3)
lines(bnpd,col="gold",lwd=3)
polygon(bbd,col="black",density=70)
polygon(vbd,col="grey",density=70)
polygon(pbd,col="red",density=70)
polygon(sjd,col="darkorange",density=70)
polygon(bnpd,col="gold",density=70)

###############################
#ordering controlled by percent of total drive for each population
BBtot<-sum(PBSout[,4])
VBtot<-sum(PBSout[,6])
PBtot<-sum(PBSout[,8])
SJtot<-sum(PBSout[,10])
BNPtot<-sum(PBSout[,12])

interest2<-c()
for (i in 1:2119){
  interest2<-(PBSout[,4]/BBtot)*100+(PBSout[,6]/VBtot)*100+(PBSout[,8]/PBtot)*100+(PBSout[,10]/SJtot)*100+(PBSout[,12]/BNPtot)*100
}

ord<-order(interest2,decreasing=TRUE)
ord2<-ord[1:100]

par(mar=c(4.2,5,4,4))
plot(PBSout[ord2,"BBsum"],col='black',pch=20,cex=3,ylim=c(0,500),ylab="Level of divergence",xlab="Region number",
     cex.lab=2,cex.axis=2)
points(PBSout[ord2,"VBsum"],col='grey',pch=20,cex=3)
points(PBSout[ord2,"PBsum"],col='red',pch=20,cex=3)
points(PBSout[ord2,"SJsum"],col='darkorange',pch=20,cex=3)
points(PBSout[ord2,"BNPsum"],col="gold",pch=20,cex=3)

legend('topright',legend=c("BB","VB","PB","SJ","BNP"),col=c("black","grey","red","darkorange2","gold"),
       pch=20,cex=2.3,bty="n",y.intersp=.5,x.intersp=.5)

size<-PBSout[,3]-PBSout[,2]
quantile(size)
sized<-density(size,na.rm=TRUE)
plot(sized)
polygon(sized,col='black',density=50)

head(PBSout[ord2,],n=30)
