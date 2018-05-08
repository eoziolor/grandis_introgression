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
BNPpbs - 0.05471878

#check NA
table(rowSums(is.na(Allpbs[subw,4:8])))

#use awk to exclude lines that are 0

cat FILE | awk '$9>0' | \
grep -v NA | \
awk '$7>q4 || $8>q5'

cat ~/analysis/fst/allpbs5kb | awk '$9>0' | \
grep -v NA | \
awk '$7>0.2164984 || $8>0.05471878' | \
~/program/bedtools2/bin/bedtools merge -i stdin -d 5000 \
-c 7,7,8,8 \
-o sum,count,sum,count \
-g <(cut -f 1-2 ~/analysis/genome/unsplit_merge.fasta.fai) > ~/analysis/fst/PBSoutliers_5kb_2interm.bed


###pipe to bedtools to merge the windows

#Put it back into R to look at it

read.table("~/analysis/fst/PBSoutliers_5kb_2interm.bed",stringsAsFactors=FALSE)->PBSout
colnames(PBSout)<- c("Scaf","start","end","SJSPsum",
"SJSPcount","BNPsum","BNPcount")

###############################
#ordering controlled by percent of total drive for each population

SJSPtot<-sum(PBSout[,4])
BNPtot<-sum(PBSout[,6])

interest2<-c()
for (i in 1:1745){
	interest2<-(PBSout[,4]/SJSPtot)*100+(PBSout[,6]/BNPtot)*100
	}
	
ord2<-order(interest2,decreasing=TRUE)

plot(PBSout[ord2,"BNPsum"],col='green',pch=20,cex=0.5)
points(PBSout[ord2,"SJSPsum"],col='blue',pch=20,cex=0.5)

size_int<-PBSout[,3]-PBSout[,2]
quantile(size_int)


#have to run size on res too to do this
hist(size_res,breaks=1000)
hist(size_int,breaks=1000)


####lim
ordlim<-ord2[1:50]

plot(PBSout[ordlim,"BBsum"],col='black',pch=20,cex=1,ylim=c(0,600),xlab="Ordered regions of interest by highest overall divergence from reference populations",
     ylab="Level of interest compared to reference populations")
points(PBSout[ordlim,"VBsum"],col='grey',pch=20,cex=1)
points(PBSout[ordlim,"PBsum"],col='red',pch=20,cex=1)
points(PBSout[ordlim,"BNPsum"],col='darkorange',pch=20,cex=1)
points(PBSout[ordlim,"SJSPsum"],col='gold',pch=20,cex=1)

legend(x=c(40,45),y=c(300,600),c("BB","VB","PB","SJ","BNP"),pch=20,cex=1,
       col=c("black","grey","red","darkorange","gold"), x.intersp=.5,y.intersp=.5)

write.table(head(PBSout[ordlim,1:3],n=30),"~/share/fst/Rscript/5kb/top30_fst",row.names=FALSE,quote=FALSE)
