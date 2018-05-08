###read all fsts in
fst<-matrix(nrow=1026857,ncol=21)

fst2<-matrix(nrow=1026857,ncol=21)

pops<-c("bbvb","bbpb","bbsj","bbbnp","bbsp","bbgb",
        "vbpb","vbsj","vbbnp","vbsp","vbgb",
        "pbsj","pbbnp","pbsp","pbgb",
        "sjbnp","sjsp","sjgb",
        "bnpsp","bnpgb",
        "spgb")

colnames(fst)<-pops

###BB vs VB
bbvbfst<-read.table("~/analysis/fst/raw/BB.VB.fst.5kb1kb.bed",stringsAsFactors=FALSE)

fst[,1]<-as.numeric(bbvbfst[,4])
fst2[,1]<-as.numeric(bbvbfst[,5])

###BB vs PB
bbpbfst<-read.table("~/analysis/fst/raw/BB.PB.fst.5kb1kb.bed",stringsAsFactors=FALSE)

fst[,2]<-as.numeric(bbpbfst[,4])
fst2[,2]<-as.numeric(bbpbfst[,5])

###BB vs SJ
bbsjfst<-read.table("~/analysis/fst/raw/BB.SJ.fst.5kb1kb.bed",stringsAsFactors=FALSE)

fst[,3]<-as.numeric(bbsjfst[,4])
fst2[,3]<-as.numeric(bbsjfst[,5])

###BB vs BNP
bbbnpfst<-read.table("~/analysis/fst/raw/BB.BNP.fst.5kb1kb.bed",stringsAsFactors=FALSE)

fst[,4]<-as.numeric(bbbnpfst[,4])
fst2[,4]<-as.numeric(bbbnpfst[,5])

###BB vs SP
bbspfst<-read.table("~/analysis/fst/raw/BB.SP.fst.5kb1kb.bed",stringsAsFactors=FALSE)

fst[,5]<-as.numeric(bbspfst[,4])
fst2[,5]<-as.numeric(bbspfst[,5])

###BB vs GB
bbgbfst<-read.table("~/analysis/fst/raw/BB.GB.fst.5kb1kb.bed",stringsAsFactors=FALSE)

fst[,6]<-as.numeric(bbspfst[,4])
fst2[,6]<-as.numeric(bbspfst[,5])

###VB vs PB
vbpbfst<-read.table("~/analysis/fst/raw/VB.PB.fst.5kb1kb.bed",stringsAsFactors=FALSE)

fst[,7]<-as.numeric(vbpbfst[,4])
fst2[,7]<-as.numeric(vbpbfst[,5])

###VB vs SJ
vbsjfst<-read.table("~/analysis/fst/raw/VB.SJ.fst.5kb1kb.bed",stringsAsFactors=FALSE)

fst[,8]<-as.numeric(vbsjfst[,4])
fst2[,8]<-as.numeric(vbsjfst[,5])

###VB vs BNP
vbbnpfst<-read.table("~/analysis/fst/raw/VB.BNP.fst.5kb1kb.bed",stringsAsFactors=FALSE)

fst[,9]<-as.numeric(vbbnpfst[,4])
fst2[,9]<-as.numeric(vbbnpfst[,5])

###VB vs SP
vbspfst<-read.table("~/analysis/fst/raw/VB.SP.fst.5kb1kb.bed",stringsAsFactors=FALSE)

fst[,10]<-as.numeric(vbspfst[,4])
fst2[,10]<-as.numeric(vbspfst[,5])

###VB vs GB
vbgbfst<-read.table("~/analysis/fst/raw/VB.GB.fst.5kb1kb.bed",stringsAsFactors=FALSE)

fst[,11]<-as.numeric(vbgbfst[,4])
fst2[,11]<-as.numeric(vbgbfst[,5])

###PB vs SJ
pbsjfst<-read.table("~/analysis/fst/raw/PB.SJ.fst.5kb1kb.bed",stringsAsFactors=FALSE)

fst[,12]<-as.numeric(pbsjfst[,4])
fst2[,12]<-as.numeric(pbsjfst[,5])

###PB vs BNP
pbbnpfst<-read.table("~/analysis/fst/raw/PB.BNP.fst.5kb1kb.bed",stringsAsFactors=FALSE)

fst[,13]<-as.numeric(pbbnpfst[,4])
fst2[,13]<-as.numeric(pbbnpfst[,5])

###PB vs SP
pbspfst<-read.table("~/analysis/fst/raw/PB.SP.fst.5kb1kb.bed",stringsAsFactors=FALSE)

fst[,14]<-as.numeric(pbspfst[,4])
fst2[,14]<-as.numeric(pbspfst[,5])

###PB vs GB
pbgbfst<-read.table("~/analysis/fst/raw/PB.GB.fst.5kb1kb.bed",stringsAsFactors=FALSE)

fst[,15]<-as.numeric(pbgbfst[,4])
fst2[,15]<-as.numeric(pbgbfst[,5])

###SJ vs BNP
sjbnpfst<-read.table("~/analysis/fst/raw/SJ.BNP.fst.5kb1kb.bed",stringsAsFactors=FALSE)

fst[,16]<-as.numeric(sjbnpfst[,4])
fst2[,16]<-as.numeric(sjbnpfst[,5])

###SJ vs SP
sjspfst<-read.table("~/analysis/fst/raw/SJ.SP.fst.5kb1kb.bed",stringsAsFactors=FALSE)

fst[,17]<-as.numeric(sjspfst[,4])
fst2[,17]<-as.numeric(sjspfst[,5])

###SJ vs GB
sjgbfst<-read.table("~/analysis/fst/raw/SJ.GB.fst.5kb1kb.bed",stringsAsFactors=FALSE)

fst[,18]<-as.numeric(sjgbfst[,4])
fst2[,18]<-as.numeric(sjgbfst[,5])

###BNP vs SP
bnpspfst<-read.table("~/analysis/fst/raw/BNP.SP.fst.5kb1kb.bed",stringsAsFactors=FALSE)

fst[,19]<-as.numeric(bnpspfst[,4])
fst2[,19]<-as.numeric(bnpspfst[,5])

###BNP vs GB
bnpgbfst<-read.table("~/analysis/fst/raw/BNP.GB.fst.5kb1kb.bed",stringsAsFactors=FALSE)

fst[,20]<-as.numeric(bnpgbfst[,4])
fst2[,20]<-as.numeric(bnpgbfst[,5])

###SP vs GB
spgbfst<-read.table("~/analysis/fst/raw/GB.SP.fst.5kb1kb.bed",stringsAsFactors=FALSE)

fst[,21]<-as.numeric(spgbfst[,4])
fst2[,21]<-as.numeric(spgbfst[,5])

######################################3
####now z score all of those

fstmeans<-colMeans(fst,na.rm=TRUE)
fstmeans<-as.matrix(fstmeans)

library(matrixStats)
fststdev<-colSds(fst,na.rm=TRUE)
fststdev<-as.matrix(fststdev)

zfst<-matrix(nrow=1026857,ncol=21)
colnames(zfst)<-pops

for (i in 1:21)
{
  for(k in 1:1026857)
  {
    zfst[k,i]<-(fst[k,i]-fstmeans[i])/fststdev[i]
  }
}

#####################################
###making sure we only have sites that are on average represented by at least 20 reads
nsnps<-fst2[,1]
for (i in 1:21){
  nsnps <- nsnps + fst2[,i]
}
nsnps <- nsnps/21

subw <- nsnps > 20

zfst_keep<-cbind(vbpbfst[,1:3],
                 zfst,
                 keep=as.numeric(subw))

write.table(zfst_keep,"~/analysis/fst/zfst_keep_5kb",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")

plot(zfst_keep[subw,8],pch=20,cex=.5,col=factor(zfst_keep[,1]))
