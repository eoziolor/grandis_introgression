###################################
#reading in the zscore of fst in R
read.table("~/analysis/data/dfst/zmerge_win_fst_pi_5kb.bed",stringsAsFactors=FALSE)->zmergeout
colnames(zmergeout)<- c("Scaf","start","end","BBVBsum", "BBVBcount","BBPBsum","BBPBcount","BBSJsum","BBSJcount","BBBNPsum",
                        "BBBNPcount","BBSPsum","BBSPcount","BBGBsum","BBGBcount","VBPBsum","VBPBcount","VBSJsum","VBSJcount","VBBNPsum","VBBNPcount",
                        "VBSPsum","VBSPcount","VBGBsum","VBGBcount","PBSJsum","PBSJcount","PBBNPsum","PBBNPcount","PBSPsum","PBSPcount","PBGBsum","PBGBcount",
                        "SJBNPsum","SJBNPcount","SJSPsum","SJSPcount","SJGBsum","SJGBcount","BNPSPsum","BNPSPcount","BNPGBsum","BNPGBcount","SPGBsum","SPGBcount")

#ordering by populations in R
#head(zmergeout[order(zmergeout$BBGBsum,decreasing=TRUE),])


################tryout to merge the distances from SP and GB

distsp<-cbind(zmergeout[1:3],zmergeout[,12],zmergeout[,22],zmergeout[,30],zmergeout[,36],zmergeout[,40])
distgb<-cbind(zmergeout[1:3],zmergeout[,14],zmergeout[,24],zmergeout[,32],zmergeout[,38],zmergeout[,42])

distnames<-c("Scaf","start","end","bb","vb","pb","sj","bnp")
colnames(distsp)<-distnames
colnames(distgb)<-distnames

colnam<-names(distsp)[4:8]

t1<-cbind(seq=seq(1:7097),distsp[4:8])
t2<-cbind(seq=seq(1:7097),distsp[4:8])

t<-cbind(zmergeout[1:3],t1[colnam]+t2[match(t1$seq,t2$seq),colnam])

# total_dist<-cbind(zmergeout[1:3],distsp[colnam]+distgb[match(distsp$Scaf,distgb$Scaf),colnam])

pbs_dist<-matrix(nrow=7097,ncol=5)

for (i in 1:5)
{
  pbs_dist[,i]<-(t[,i+3]-zmergeout[,44])/2
}

pbs_dist<-cbind(zmergeout[,1:3],pbs_dist)
colnames(pbs_dist)<-distnames

write.table(pbs_dist,"~/analysis/data/dfst/pbs_true_fst_pi_5kb",sep="\t")

################plotting fake pbs_dist

pbs_dist<-read.table("~/analysis/data/dfst/pbs_true_fst_pi_5kb",sep="\t")
distnames<-c("Scaf","start","end","bb","vb","pb","sj","bnp")
colnames(pbs_dist)<-distnames

########ordering by percent interest

BBtot<-sum(pbs_dist[,4])
VBtot<-sum(pbs_dist[,5])
PBtot<-sum(pbs_dist[,6])
SJSPtot<-sum(pbs_dist[,7])
BNPtot<-sum(pbs_dist[,8])

interest<-c()
for (i in 1:7097){
  interest<-(pbs_dist[,4]/BBtot)*100+(pbs_dist[,5]/VBtot)*100+(pbs_dist[,6]/PBtot)*100+(pbs_dist[,7]/SJSPtot)*100+(pbs_dist[,8]/BNPtot)*100
}

ord<-order(interest,decreasing=TRUE)

par(mfrow=c(1,1))
plot(pbs_dist[ord,4],col='black',pch=20,cex=0.8)
points(pbs_dist[ord,5],col='grey',pch=20,cex=0.8)
points(pbs_dist[ord,6],col='red',pch=20,cex=0.8)
points(pbs_dist[ord,7],col='darkorange',pch=20,cex=0.8)
points(pbs_dist[ord,8],col='gold',pch=20,cex=0.8)

#############################
#superlimited

ordlim<-ord[c(1:100)]
plot(pbs_dist[ordlim,4],col='black',pch=20,cex=0.8)
points(pbs_dist[ordlim,5],col='grey',pch=20,cex=0.8)
points(pbs_dist[ordlim,6],col='red',pch=20,cex=0.8)
points(pbs_dist[ordlim,7],col='darkorange',pch=20,cex=0.8)
points(pbs_dist[ordlim,8],col='gold',pch=20,cex=0.8)


###################

#for a figure
par(mfrow=c(1,1),mar=c(6,6,2,2))
plot(pbs_dist[ordlim,4],col='black',pch=20,cex=2,
     xlab="Ordered regions of interest by highest overall divergence from reference populations",
     ylab="Level of interest compared to reference populations",bg="black",ylim=c(0,4300),cex.lab=1.5,cex.axis=1.5)
points(pbs_dist[ordlim,5],col='grey',pch=20,cex=2,bg="grey")
points(pbs_dist[ordlim,6],col='red',pch=20,cex=2,bg="white")
points(pbs_dist[ordlim,7],col='darkorange',pch=20,cex=2,bg="red")
points(pbs_dist[ordlim,8],col='gold',pch=20,cex=2,bg="pink")

legend(x=c(80,92),y=c(2800,4200),c("BB","VB","PB","SJSP","BNP"),pch=20,cex=1.7,
       col=c("black","grey","red","darkorange","gold"), x.intersp=.3,y.intersp=.6)


write.table(pbs_dist[ordlim,1:3],"~/analysis/data/dfst/top50.bed",row.names=FALSE,quote=FALSE,col.names=FALSE,sep="\t")




#forget the rest



########Order 3

BBtot<-sum(pbs_dist[,4])
VBtot<-sum(pbs_dist[,5])
PBtot<-sum(pbs_dist[,6])
SJSPtot<-sum(pbs_dist[,7])
BNPtot<-sum(pbs_dist[,8])

interest3<-c()
for (i in 1:2345){
  interest3<-(pbs_dist[,4]/BBtot)*100+(pbs_dist[,5]/VBtot)*100+(pbs_dist[,6]/PBtot)*100
}

ord3<-order(interest3,decreasing=TRUE)

plot(pbs_dist[ord3,4],col='black',pch=20,cex=1.2,
     xlab="Ordered regions of interest by highest overall divergence from reference populations",
     ylab="Level of interest compared to reference populations",bg="black",ylim=c(-4000,4000))
points(pbs_dist[ord3,5],col='grey',pch=20,cex=1.2,bg="grey")
points(pbs_dist[ord3,6],col='red',pch=20,cex=1.2,bg="white")
points(pbs_dist[ord3,7],col='darkorange',pch=20,cex=1.2,bg="red")
points(pbs_dist[ord3,8],col='gold',pch=20,cex=1.2,bg="pink")

legend(x=c(24,26.5),y=c(160,220),c("BB","VB","PB","SJSP","BNP"),pch=20,cex=1,
       col=c("black","grey","red","darkorange","gold"), x.intersp=.5,y.intersp=.7)

###FINAL FIGURE

par(mfrow=c(1,1))
ordsuperlim<-ord3[c(1:50)]
plot(pbs_dist[ordsuperlim,4],col='black',pch=20,cex=1.2,
     xlab="Ordered regions of interest by highest overall divergence from reference populations",
     ylab="Level of interest compared to reference populations",bg="black",ylim=c(0,4000))
points(pbs_dist[ordsuperlim,5],col='grey',pch=20,cex=1.2,bg="grey")
points(pbs_dist[ordsuperlim,6],col='red',pch=20,cex=1.2,bg="white")
points(pbs_dist[ordsuperlim,7],col='darkorange',pch=20,cex=1.2,bg="red")
points(pbs_dist[ordsuperlim,8],col='gold',pch=20,cex=1.2,bg="pink")

legend(x=c(35,40),y=c(2500,3500),c("BB","VB","PB","SJSP","BNP"),pch=20,cex=1,
       col=c("black","grey","red","darkorange","gold"), x.intersp=.5,y.intersp=.7)

write.table(head(pbs_dist[ordsuperlim,1:3],n=30),"~/share/d_fst_analysis/top30_5kb_fst_pi",row.names=FALSE,quote=FALSE)

#for phylogeny
write.table(pbs_dist[,1],"~/share/d_fst_analysis/5kb/outlier_justify/both_out",col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
