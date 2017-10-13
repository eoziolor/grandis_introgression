mf<-read.table("~/analysis/data/fst/raw/M.F.fst.5kb1kb.bed", header=FALSE,stringsAsFactors = FALSE)
mf[,4]<-as.numeric(mf[,4])
subw<-mf[,5]>20
quantile(mf[subw,4],probs=.99,na.rm=TRUE)
plot(mf[subw,4],pch=20,cex=.5,col=ifelse(mf[,4]<0.01453576,"black","red"))
