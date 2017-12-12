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
gbord<-order(gb[,2],decreasing=TRUE)

sp<-k2ord[49:96,]
spord<-order(sp[,2],decreasing=TRUE)

bnp<-k2ord[97:144,]
bnpord<-order(bnp[,2],decreasing=TRUE)

sj<-k2ord[145:168,]
sjord<-order(sj[,2],decreasing=TRUE)

pb<-k2ord[168:215,]
pbord<-order(pb[,2],decreasing=TRUE)

vb<-k2ord[215:264,]
vbord<-order(vb[,2],decreasing=TRUE)

bb<-k2ord[265:288,]
bbord<-order(bb[,2],decreasing=TRUE)

kall<-rbind(bb[bbord,],vb[vbord,],pb[pbord,],sj[sjord,],bnp[bnpord,],gb[gbord,],sp[spord,])

kall<-na.omit(kall)

par(mfrow=c(1,1),mar=c(3,3,1,1))
barplot(t(as.matrix(kall)),col=c("cadetblue2","black"), border=NA, xaxt="n",space=0,
        cex.lab=2,cex.axis=1.9)

ax<-c(12,48,96,132,168,216,264)
ax2<-c(0,120,121,192,193,288)
vert<-c(24,74,121,144,192,240)
axnames<-c("BB","VB","PB", "SJ","BNP","SP","GB")
axcol<-c("black","firebrick2","cadetblue3")
axis(side=1,at= ax[1:3], labels = axnames[1:3],tck=-.02,lwd=0,
     col.axis=axcol[1],col=axcol[1],cex.axis=1.6)
axis(side=1,at= ax[4:5], labels = axnames[4:5],tck=-.02,lwd=0,
     col.axis=axcol[2],col=axcol[2],cex.axis=1.6)
axis(side=1,at= ax[6:7], labels = axnames[6:7],tck=-.02,lwd=0,
     col.axis=axcol[3],col=axcol[3],cex.axis=1.6)
#Just the lines
axis(side=1,at= ax2[1:2], labels=c("",""), tck=0,lwd=4,
     col.axis=axcol[1],col=axcol[1],cex.axis=1.4)
axis(side=1,at= ax2[3:4], labels = c("",""),tck=0,lwd=4,
     col.axis=axcol[2],col=axcol[2],cex.axis=1.4)
axis(side=1,at= ax2[5:6], labels = c("",""),tck=-0,lwd=4,
     col.axis=axcol[3],col=axcol[3],cex.axis=1.4)
abline(v=vert,col="khaki2",lty=1,lwd=2.5)

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
     main="",xaxt='n',
     cex.axis=1.2,lwd=2,bty="l")
lines(vbh,col="black",lwd=2)
lines(pbh,col="black",lwd=2)

legend('topright',legend=c("resistant","intermediate","reference"),col=c("black","red","blue"),
       pch=20,cex=1.5,bty="n",y.intersp=1,x.intersp=.5)

plot(sjh,col="red",ylim=c(0,48),xlim=c(-.2,1.2),main="",
     cex.axis=1.2,lwd=2,bty="l")
lines(bnph,col="red",lwd=2)
lines(sph,col="blue",lwd=2)
lines(gbh,col="blue",lwd=2)
# 
# legend('topleft',legend=c("intermediate","reference"),col=c("red","blue"),
#        pch=20,cex=1.5,bty="n",y.intersp=1,x.intersp=.5)


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