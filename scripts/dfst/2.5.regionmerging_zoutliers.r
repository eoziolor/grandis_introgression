pbs<-read.table("~/analysis/data/dfst/zpbs_5kb.bed",header=FALSE,stringsAsFactors = FALSE)
pbsname<-c("Scaf","start","end","BBzpbs","VBzpbs","PBzpbs","SJzpbs","BNPzpbs")
colnames(pbs)<-pbsname

col<-c()
for (i in 1:5){
  col[i]<-quantile(pbs[,i+3],prob=.99,na.rm=TRUE)
}

1% thresholds:

BB=2.608447
VB=2.118201
PB=2.738914
SJ=2.473380
BNP=2.075014


cat ~/analysis/data/dfst/zpbs_5kb.bed | grep -v NA | \
awk '$4>2.608447 || $5>2.118201  || $6>2.738914 || $7>2.473380 || $8>2.075014' | \
~/program/bedtools2/bin/bedtools merge -i stdin -d 50000 \
-c 4,4,5,5,6,6,7,7,8,8 \
-o max,count,max,count,max,count,max,count,max,count \
-g <(cut -f 1-2 ~/analysis/data/genome/unsplit_merge.fasta.fai) > ~/analysis/data/dfst/zregions_max.bed
