ord_temp<-read.table("~/analysis/admixture/seqlist.txt",header=FALSE)
ord<-unlist(ord_temp)

k2<-read.table("~/analysis/admixture/all/allscaf.2.Q")
k2ord<-k2[ord,]
barplot(t(as.matrix(k2ord)),col=c("deepskyblue2","black"), ylab="Ancestry",border=NA, xaxt="n",space=0,
        cex.lab=1.5,cex.axis=1.5)

abline(b=0,v=c(0,48,96,144,168,215,264,288),col="grey",lwd=3,lty=5)
abline(a=0,b=0,col="black",lwd=3)

ax<-c(24,72,124,156,192,240,276)
axnames<-c("SP","GB","BNP", "SJ","PB","VB","BB")
axis(side=1,at= ax, labels = axnames,tck=-.03)

########################################################
k3<-read.table("~/analysis/admixture/all/allscaf.3.Q")
k3ord<-k3[ord,]
barplot(t(as.matrix(k3ord)),col=c("black","deepskyblue2","red"), ylab="Ancestry",border=NA, xaxt="n",space=0)

abline(b=0,v=c(0,48,96,144,168,215,264,288),col="grey",lwd=3,lty=5)

ax<-c(24,72,124,156,192,240,276)
axnames<-c("SP","GB","BNP", "SJ","PB","VB","BB")
axis(side=1,at= ax, labels = axnames,tck=-.03)

##########################################################
k4<-read.table("~/analysis/admixture/all/allscaf.4.Q")
k4ord<-k4[ord,]
barplot(t(as.matrix(k4ord)),col=c("red","deepskyblue2","black","darkorange"), ylab="Ancestry",border=NA, xaxt="n",space=0,
        cex.lab=1.5,cex.axis = 1.5)

abline(b=0,v=c(0,48,96,144,168,215,264,288),col="grey",lwd=3,lty=5)
abline(a=0,b=0,col="black",lwd=3)

ax<-c(24,72,124,156,192,240,276)
axnames<-c("SP","GB","BNP", "SJ","PB","VB","BB")
axis(side=1,at= ax, labels = axnames,tck=-.03)
