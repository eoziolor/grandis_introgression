ord_temp<-read.table("~/analysis/data/fastngs/order.txt",header=FALSE)
ord<-unlist(ord_temp)

k2<-read.table("~/analysis/data/fastngs/combined_all")
k2ord<-k2[ord,]
barplot(t(as.matrix(k2ord)),col=c("deepskyblue2","black"), ylab="Ancestry",border=NA, xaxt="n",space=0,
        cex.lab=1.5,cex.axis=1.5)

abline(b=0,v=c(0,48,96,144,168,215,264,288),col="grey",lwd=3,lty=5)
abline(a=0,b=0,col="black",lwd=3)

ax<-c(24,72,124,156,192,240,276)
axnames<-c("GB","SP","BNP", "SJ","PB","VB","BB")
axis(side=1,at= ax, labels = axnames,tck=-.03)

gb<-k2ord[0:48,]
gbord<-order(gb[,2],decreasing=FALSE)

sp<-k2ord[49:96,]
spord<-order(sp[,2],decreasing=FALSE)

bnp<-k2ord[97:144,]
bnpord<-order(bnp[,2],decreasing=FALSE)

sj<-k2ord[145:168,]
sjord<-order(sj[,2],decreasing=FALSE)

pb<-k2ord[168:215,]
pbord<-order(pb[,2],decreasing=FALSE)

vb<-k2ord[215:264,]
vbord<-order(vb[,2],decreasing=FALSE)

bb<-k2ord[265:288,]
bbord<-order(bb[,2],decreasing=FALSE)

kall<-rbind(gb[gbord,],sp[spord,],bnp[bnpord,],sj[sjord,],pb[pbord,],vb[vbord,],bb[bbord,])

kall<-na.omit(kall)

par(mfrow=c(1,1))
barplot(t(as.matrix(kall)),col=c("deepskyblue2","black"), ylab="Likelihood of belonging to reference or resistant genotype",
        border=NA, xaxt="n",space=0,
        cex.lab=1.5,cex.axis=1.5)

abline(b=0,v=c(0,48,96,144,168,215,264,288),col="grey50",lwd=3)
abline(a=0,b=0,col="black",lwd=3)
abline(h=c(0.25,.5,.75),col="purple",lwd=1.5,lty=5)

ax<-c(24,72,124,156,192,240,276)
axnames<-c("GB","SP","BNP", "SJ","PB","VB","BB")
axis(side=1,at= ax, labels = axnames,tck=-.03)

###histograms and densities

gbh<-density(gb[,1],na.rm=TRUE)
sph<-density(sp[,1],na.rm=TRUE)
bnph<-density(bnp[,1],na.rm=TRUE)
sjh<-density(sj[,1],na.rm=TRUE)
pbh<-density(pb[,1],na.rm=TRUE)
vbh<-density(vb[,1],na.rm=TRUE)
bbh<-density(bb[,1],na.rm=TRUE)

plot(bbh,col="black",ylim=c(0,48),xlim=c(-.2,1.2))
lines(vbh,col="grey50")
lines(pbh,col="red")
lines(sjh,col="darkorange")
lines(bnph,col="gold")
lines(sph,col="cyan")
lines(gbh,col="blue")

#separate

par(mfrow=c(2,1),mar=c(2,4,1,1))

plot(bbh,col="black",ylim=c(0,7),xlim=c(-.2,1.2),
     main="Bimodality suggests admixture through migration",xaxt='n',
     cex.axis=1.2,lwd=2)
lines(vbh,col="grey50",lwd=2)
lines(pbh,col="red",lwd=2)

legend('topright',legend=c("BB","VB","PB"),col=c("black","grey","red","darkorange2","gold"),
       pch=20,cex=1.2,bty="n",y.intersp=.5,x.intersp=.5)

plot(sjh,col="darkorange",ylim=c(0,48),xlim=c(-.2,1.2),main='K0 is reference genotype, K1 is resistant',
     cex.axis=1.2,lwd=2)
lines(bnph,col="gold",lwd=2)
lines(sph,col="cyan",lwd=2)
lines(gbh,col="blue",lwd=2)

legend('topright',legend=c("SJ","BNP","SP","GB"),col=c("darkorange2","gold","cyan","blue"),
       pch=20,cex=1.2,bty="n",y.intersp=.5,x.intersp=.5)


###
#install.packages("diptest")
library(diptest)

dip.test(bb[,1],simulate.p.value = TRUE,B=2000)
#p=0.92
dip.test(vb[,1],simulate.p.value = TRUE,B=2000)
#p=0.22
dip.test(pb[,1],simulate.p.value = TRUE,B=2000)
#p=0.0015
dip.test(sj[,1],simulate.p.value = TRUE,B=2000)
#p=0.99
dip.test(bnp[,1],simulate.p.value = TRUE,B=2000)
#p=0.99
dip.test(sp[,1],simulate.p.value = TRUE,B=2000)
#p=1
dip.test(gb[,1],simulate.p.value = TRUE,B=2000)
#p=1