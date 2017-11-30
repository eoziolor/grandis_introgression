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

#check NA
table(rowSums(is.na(pbs[subw,4:8])))

#use awk to exclude lines that are 0

cat FILE | awk '$9>0' | \
grep -v NA | \
awk '$4<q1 || $5<q2 || $6<q3 || $7<q4 || $8<q5'

cat ~/analysis/data/fst/allpbs5kb | awk '$9>0' | \
grep -v NA | \
awk '$4<0.2533448 || $5<0.1461146  || $6<0.1765786 || $7<0.2164984 || $8<0.05471878' | \
~/program/bedtools2/bin/bedtools merge -i stdin -d 5000 \
-c 4,4,5,5,6,6,7,7,8,8 \
-o sum,count,sum,count,sum,count,sum,count,sum,count \
-g <(cut -f 1-2 ~/analysis/data/genome/unsplit_merge.fasta.fai) > ~/analysis/data/fst/PBS_nonoutliers.bed


neut<-read.table("~/analysis/data/fst/PBS_nonoutliers.bed",stringsAsFactors = FALSE)
write.table(neut[,1:3],"~/analysis/data/fst/list_neutral",row.names = FALSE,col.names = FALSE,quote = FALSE)
