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

write.table(bb[,1:3],"~/analysis/data/fst/individual_pbs/bb_pbs.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(vb[,1:3],"~/analysis/data/fst/individual_pbs/vb_pbs.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(pb[,1:3],"~/analysis/data/fst/individual_pbs/pb_pbs.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(sj[,1:3],"~/analysis/data/fst/individual_pbs/sj_pbs.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(bnp[,1:3],"~/analysis/data/fst/individual_pbs/bnp_pbs.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)
