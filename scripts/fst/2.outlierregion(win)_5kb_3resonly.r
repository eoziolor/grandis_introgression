#Open fsttopbs first and load the PBS files
#Allpbs<-cbind(fst[[1]][,1:3],bbpbs=BBpbs,vbpbs=VBpbs,pbpbs=PBpbs,sjpbs=SJpbs,
#              bnppbs=BNPpbs,keep=as.numeric(subw)) ###forallpops###,


#write.table(Allpbs,file="~/share/fst/Allpbs50kb",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")

quantile(Allpbs[,4],probs=.99,na.rm=TRUE)
quantile(Allpbs[,5],probs=.99,na.rm=TRUE)
quantile(Allpbs[,6],probs=.99,na.rm=TRUE)
quantile(Allpbs[,7],probs=.99,na.rm=TRUE)
quantile(Allpbs[,8],probs=.99,na.rm=TRUE)


Quantile info
BBpbs - 0.2533448
VBpbs - 0.1461146 
PBpbs - 0.1765786 
SJpbs - 0.2164984
BNP - 0.05471878

#check NA
table(rowSums(is.na(Allpbs[subw,4:8])))

#use awk to exclude lines that are 0

cat FILE | awk '$9>0' | \
grep -v NA | \
awk '$4>q1 || $5>q2 || $6>q3 || $7>q4 || $8>q5'

cat ~/analysis/fst/allpbs5kb | awk '$9>0' | \
grep -v NA | \
awk '$4>0.2533448 || $5>0.1461146  || $6>0.1765786' | \
~/program/bedtools2/bin/bedtools merge -i stdin -d 5000 \
-c 4,4,5,5,6,6 \
-o sum,count,sum,count,sum,count \
-g <(cut -f 1-2 ~/analysis/genome/unsplit_merge.fasta.fai) > ~/analysis/fst/PBSoutliers_5kb_3resonly.bed


###pipe to bedtools to merge the windows

#Put it back into R to look at it

read.table("~/analysis/fst/PBSoutliers_5kb_3resonly.bed",stringsAsFactors=FALSE)->PBSout
colnames(PBSout)<- c("Scaf","start","end","BBsum", "BBcount","VBsum","VBcount","PBsum","PBcount")

###############################
#ordering controlled by percent of total drive for each population
BBtot<-sum(PBSout[,4])
VBtot<-sum(PBSout[,6])
PBtot<-sum(PBSout[,8])


interest2<-c()
for (i in 1:755){
  interest2<-(PBSout[,4]/BBtot)*100+(PBSout[,6]/VBtot)*100+(PBSout[,8]/PBtot)*100
}

ord2<-order(interest2,decreasing=TRUE)
ord3<-ord2[1:50]
plot(PBSout[ord3,"BBsum"],col='black',pch=20,cex=1,ylim=c(0,500))
points(PBSout[ord3,"VBsum"],col='red',pch=20,cex=1)
points(PBSout[ord3,"PBsum"],col='orange',pch=20,cex=1)


size_res<-PBSout[,3]-PBSout[,2]
quantile(size_res)

