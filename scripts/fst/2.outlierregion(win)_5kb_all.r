#Open fsttopbs first and load the PBS files
#Allpbs<-cbind(fst[[1]][,1:3],bbpbs=BBpbs,vbpbs=VBpbs,pbpbs=PBpbs,sjpbs=SJpbs,
#              bnppbs=BNPpbs,keep=as.numeric(subw)) ###forallpops###,


#write.table(Allpbs,file="~/share/data/fst/Allpbs50kb",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")

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
awk '$4>q1 || $5>q2 || $6>q3 || $7>q4 || $8>q5'

cat ~/analysis/data/fst/allpbs5kb | awk '$9>0' | \
grep -v NA | \
awk '$4>0.2533448 || $5>0.1461146  || $6>0.1765786 || $7>0.2164984 || $8>0.05471878' | \
~/program/bedtools2/bin/bedtools merge -i stdin -d 50000 \
-c 4,4,5,5,6,6,7,7,8,8 \
-o max,count,max,count,max,count,max,count,max,count \
-g <(cut -f 1-2 ~/analysis/data/genome/unsplit_merge.fasta.fai) > ~/analysis/data/fst/PBSoutliers_5kb_all_max.bed


###pipe to bedtools to merge the windows

#Put it back into R to look at it

read.table("~/analysis/data/fst/PBSoutliers_5kb_all.bed",stringsAsFactors=FALSE)->PBSout
colnames(PBSout)<- c("Scaf","start","end","BBsum", "BBcount","VBsum","VBcount","PBsum","PBcount","SJsum","SJcount","BNPsum","BNPcount")

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


###For region discovery

write.table(PBSout[,1:3],"~/analysis/data/fst/individual_pbs/all_pbs.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)
