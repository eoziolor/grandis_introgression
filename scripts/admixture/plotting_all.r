ord_temp<-read.table("~/analysis/data/admixture/seqlist.txt",header=FALSE)
ord<-unlist(ord_temp)

k2<-read.table("~/analysis/data/admixture/all/allscaf.2.Q")
k2ord<-k2[ord,]
barplot(t(as.matrix(k2ord)),col=c("deepskyblue2","black"), ylab="Ancestry",border=NA, xaxt="n",space=0,
        cex.lab=1.5,cex.axis=1.5)

abline(b=0,v=c(0,48,96,144,168,215,264,288),col="grey",lwd=3,lty=5)
abline(a=0,b=0,col="black",lwd=3)

ax<-c(24,72,124,156,192,240,276)
axnames<-c("SP","GB","BNP", "SJ","PB","VB","BB")
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
ax2<-c(0,121,122,193,194,288)
vert<-c(24,74,122,146,194,242)
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

########################################################
k3<-read.table("~/analysis/data/admixture/all/allscaf.3.Q")
k3ord<-k3[ord,]
barplot(t(as.matrix(k3ord)),col=c("black","deepskyblue2","red"), ylab="Ancestry",border=NA, xaxt="n",space=0)

abline(b=0,v=c(0,48,96,144,168,215,264,288),col="grey",lwd=3,lty=5)

ax<-c(24,72,124,156,192,240,276)
axnames<-c("SP","GB","BNP", "SJ","PB","VB","BB")
axis(side=1,at= ax, labels = axnames,tck=-.03)

##########################################################
k4<-read.table("~/analysis/data/admixture/all/allscaf.4.Q")
k4ord<-k4[ord,]

gb<-k4ord[0:48,]
gbord<-order(gb[,2],decreasing=TRUE)

sp<-k4ord[49:96,]
spord<-order(sp[,2],decreasing=TRUE)

bnp<-k4ord[97:144,]
bnpord<-order(bnp[,3],decreasing=TRUE)

sj<-k4ord[145:168,]
sjord<-order(sj[,3],decreasing=TRUE)

pb<-k4ord[168:215,]
pbord<-order(pb[,3],decreasing=TRUE)

vb<-k4ord[215:264,]
vbord<-order(vb[,3],decreasing=TRUE)

bb<-k4ord[265:288,]
bbord<-order(bb[,3],decreasing=TRUE)

kall<-rbind(bb[bbord,],vb[vbord,],pb[pbord,],sj[sjord,],bnp[bnpord,],gb[gbord,],sp[spord,])

kall<-na.omit(kall)
barplot(t(as.matrix(kall)),col=c("red","deepskyblue2","black","darkorange"), ylab="Ancestry",border=NA, xaxt="n",space=0,
        cex.lab=1.5,cex.axis = 1.5)

ax<-c(12,48,96,132,168,216,264)
ax2<-c(0,121,122,193,194,288)
vert<-c(24,74,122,146,194,242)
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

