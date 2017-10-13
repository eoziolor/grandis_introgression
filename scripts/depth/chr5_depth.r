names<-c("chr","start","end","mean","median","count")
bb<-read.table("~/analysis/data/depth/BB_5kb_chr5.cov",header=FALSE)
colnames(bb)<-names
vb<-read.table("~/analysis/data/depth/VB_5kb_chr5.cov",header=FALSE)
colnames(vb)<-names
pb<-read.table("~/analysis/data/depth/PB_5kb_chr5.cov",header=FALSE)
colnames(pb)<-names
sj<-read.table("~/analysis/data/depth/SJ_5kb_chr5.cov",header=FALSE)
colnames(sj)<-names
bnp<-read.table("~/analysis/data/depth/BNP_5kb_chr5.cov",header=FALSE)
colnames(bnp)<-names
sp<-read.table("~/analysis/data/depth/SP_5kb_chr5.cov",header=FALSE)
colnames(sp)<-names
gb<-read.table("~/analysis/data/depth/GB_5kb_chr5.cov",header=FALSE)
colnames(gb)<-names

plot(log(2*bb[14000:16000,6]),pch=20,cex=0.7,col="black",ylab="log raw read coverage over 5kb chunks",xlab="chunk between 600-1000kb inside chr5",
     cex.lab=1.5,cex.axis=1.5,ylim=c(3,10))
points(log(vb[14000:16000,6]),pch=20,cex=0.7,col="grey")
points(log(pb[14000:16000,6]),pch=20,cex=0.7,col="red")
points(log(2*sj[14000:16000,6]),pch=20,cex=0.7,col="darkorange")
points(log(bnp[14000:16000,6]),pch=20,cex=0.7,col="gold")
points(log(sp[14000:16000,6]),pch=20,cex=0.7,col="cyan")
points(log(gb[14000:16000,6]),pch=20,cex=0.7,col="blue")


plot(log(2*bb[,6]),pch=20,cex=0.7,col="black", ylab="log raw read coverage over 5kb chunks",xlab="chr5 position",
     cex.lab=1.5,cex.axis=1.5)
points(log(vb[,6]),pch=20,cex=0.7,col="grey")
points(log(pb[,6]),pch=20,cex=0.7,col="red")
points(log(2*sj[,6]),pch=20,cex=0.7,col="darkorange")
points(log(bnp[,6]),pch=20,cex=0.7,col="gold")
points(log(sp[,6]),pch=20,cex=0.7,col="cyan")
points(log(gb[,6]),pch=20,cex=0.7,col="blue")
