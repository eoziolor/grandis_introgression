theta<-read.table("~/analysis/data/angsd/thetas_neut_5kb",header=TRUE, sep=',')
pi<-read.table("~/analysis/data/angsd/pi_neut_5kb", header=TRUE, sep=',')
taj<-read.table("~/analysis/data/angsd/taj",header=TRUE, sep=',')

#Theta no sliding
#t_noslid<-matrix(nrow=22460,ncol=7)
#colnoslid<-colnames(theta[,4:10])

#for (i in 1:7){
#  k=0
#  for (j in seq(from=1, to=112300,by=5)){
#    k=k+1
#    t_noslid[k,i]<-theta[j,i+3]
#  }
#}
#colnames(t_noslid)<-colnoslid

#dim(t_noslid)

#theta_full<-matrix(nrow=1,ncol=7)
#for (i in 1:7){
#  theta_full[,i]<-median(t_noslid[,i],na.rm=TRUE)
#}
#colnames(theta_full)<-colnoslid
#theta_full

# par(mfrow=c(3,3))
# hist(theta[,4],breaks=3000,xlim=c(0.002,.025),col="black",border="black")
# hist(theta[,5],breaks=3000,xlim=c(0.002,.025),col="grey40",border="black")
# hist(theta[,6],breaks=3000,xlim=c(0.002,.025),col="grey80",border="black")
# hist(theta[,7],breaks=3000,xlim=c(0.002,.025),col="red",border="firebrick2")
# hist(theta[,8],breaks=3000,xlim=c(0.002,.025),col="lightpink",border="firebrick2")
# hist(theta[,9],breaks=3000,xlim=c(0.002,.025),col="cadetblue3",border="cadetblue3")
# hist(theta[,10],breaks=3000,xlim=c(0.002,.025),col="cadetblue1",border="cadetblue3")

pimean<-c()
for(i in 1:7){
  pimean[i]<-mean(pi[,i+3],na.rm=TRUE)
}

pimedian<-c()
for(i in 1:7){
  pimedian[i]<-median(pi[,i+3],na.rm=TRUE)
}
pops<-c("BB","VB","PB","SJ","BNP","SP","GB")
names(pimean)<-pops
names(pimedian)<-pops

bbt<-density(theta[,4],na.rm=TRUE)
vbt<-density(theta[,5],na.rm=TRUE)
pbt<-density(theta[,6],na.rm=TRUE)
sjt<-density(theta[,7],na.rm=TRUE)
bnpt<-density(theta[,8],na.rm=TRUE)
spt<-density(theta[,9],na.rm=TRUE)
gbt<-density(theta[,10],na.rm=TRUE)

par(mfrow=c(3,1),mar=c(4,5,2,2),mgp=c(3,2,0))
plot(bbt,xlim=c(0,.025),col="black",bty="l",ylim=c(0,450),cex.lab=2,cex.lab=2,xlab="",main="",ylab="",lwd=3,cex.axis=3)
polygon(bbt,col="black",density=100,border=NA)
lines(vbt,xlim=c(0,.025),col="black",lwd=3)
polygon(vbt,col="grey40",density=100,border=NA)
lines(pbt,xlim=c(0,.025),col="black",lwd=3)
polygon(pbt,col="grey80",density=100,border=NA)
lines(sjt,xlim=c(0,.025),col="firebrick2",lwd=3)
polygon(bnpt,col="lightpink",density=100,border=NA)
lines(bnpt,xlim=c(0,.025),col="firebrick2",lwd=3)
polygon(spt,col="cadetblue3",density=100,border=NA)
lines(spt,xlim=c(0,.025),col="cadetblue3",lwd=3)
polygon(sjt,col="red",density=100,border=NA)
lines(gbt,xlim=c(0,.025),col="cadetblue3",lwd=3)
polygon(gbt,col="cadetblue1",density=100,border=NA)

# legend('topright',legend=c("BB","VB","PB","SJ","BNP","SP","GB"),col=c("black","grey","red","darkorange2","gold","cyan","blue"),
#        pch=20,cex=1.8,bty="n",y.intersp=.5,x.intersp=.5)

#Pi no sliding
#p_noslid<-matrix(nrow=22460,ncol=7)
#colnoslid<-colnames(pi[,4:10])

#for (i in 1:7){
#  k=0
#  for (j in seq(from=1, to=112300,by=5)){
#    k=k+1
#    p_noslid[k,i]<-pi[j,i+3]
# }
#}
#colnames(p_noslid)<-colnoslid

#dim(p_noslid)

#pi_full<-matrix(nrow=1,ncol=7)
#for (i in 1:7){
#  pi_full[,i]<-median(p_noslid[,i],na.rm=TRUE)
#}
#colnames(pi_full)<-colnoslid
#pi_full

# par(mfrow=c(3,3))
# hist(pi[,4],breaks=5000,xlim=c(0,.025),col="black",border="black")
# hist(pi[,5],breaks=5000,xlim=c(0,.025),col="grey",border="grey")
# hist(pi[,6],breaks=5000,xlim=c(0,.025),col="red",border="red")
# hist(pi[,7],breaks=5000,xlim=c(0,.025),col="orange",border="orange")
# hist(pi[,8],breaks=5000,xlim=c(0,.025),col="yellow",border="yellow")
# hist(pi[,9],breaks=5000,xlim=c(0,.025),col="green",border="green")
# hist(pi[,10],breaks=5000,xlim=c(0,.025),col="blue",border="blue")

bbp<-density(pi[,4],na.rm=TRUE)
vbp<-density(pi[,5],na.rm=TRUE)
pbp<-density(pi[,6],na.rm=TRUE)
sjp<-density(pi[,7],na.rm=TRUE)
bnpp<-density(pi[,8],na.rm=TRUE)
spp<-density(pi[,9],na.rm=TRUE)
gbp<-density(pi[,10],na.rm=TRUE)

plot(bbp,xlim=c(0.001,.025),col="black",bty="l",ylim=c(0,450),cex.lab=2,cex.lab=2,xlab="",ylab="",lwd=3,main="",cex.axis=3)
polygon(bbp,col="black",density=100,border=NA)
lines(vbp,xlim=c(0.001,.025),col="black",lwd=3)
polygon(vbp,col="grey40",density=100,border=NA)
lines(pbp,xlim=c(0.001,.025),col="black",lwd=3)
polygon(pbp,col="grey80",density=100,border=NA)
lines(sjp,xlim=c(0.001,.025),col="firebrick2",lwd=3)
polygon(sjp,col="red",density=100,border=NA)
lines(bnpp,xlim=c(0.001,.025),col="firebrick2",lwd=3)
polygon(bnpp,col="lightpink",density=100,border=NA)
lines(spp,xlim=c(0.001,.025),col="cadetblue3",lwd=3)
polygon(spp,col="cadetblue3",density=100,border=NA)
lines(gbp,xlim=c(0.001,.025),col="cadetblue3",lwd=3)
polygon(gbp,col="cadetblue1",density=100,border=NA)

###Tajima's D

# par(mfrow=c(1,1))
# hist(taj[,4],breaks=5000,xlim=c(-.3,.1),col="black",border="black")
# hist(taj[,5],breaks=5000,xlim=c(-.3,.1),col="grey",border="grey")
# hist(taj[,6],breaks=5000,xlim=c(-.3,.1),col="red",border="red")
# hist(taj[,7],breaks=5000,xlim=c(-.3,.1),col="orange",border="orange")
# hist(taj[,8],breaks=5000,xlim=c(-.3,.1),col="yellow",border="yellow")
# hist(taj[,9],breaks=5000,xlim=c(-.3,.1),col="green",border="green")
# hist(taj[,10],breaks=5000,xlim=c(-.3,.1),col="blue",border="blue")

bbtaj<-density(taj[,4],na.rm=TRUE)
vbtaj<-density(taj[,5],na.rm=TRUE)
pbtaj<-density(taj[,6],na.rm=TRUE)
sjtaj<-density(taj[,7],na.rm=TRUE)
bnptaj<-density(taj[,8],na.rm=TRUE)
sptaj<-density(taj[,9],na.rm=TRUE)
gbtaj<-density(taj[,10],na.rm=TRUE)

#par(mfrow=c(1,1),mgp=c(3,2,0))
plot(bbtaj,xlim=c(-.15,.3),col="black",bty="l",ylim=c(0,10.5),xlab="",ylab="",main="",lwd=3,cex.axis=3)
polygon(bbtaj,col="black",density=100,border=NA)
lines(vbtaj,xlim=c(-.15,.3),col="black",lwd=3)
polygon(vbtaj,col="grey40",density=100,border=NA)
lines(pbtaj,xlim=c(-.15,.3),col="black",lwd=3)
polygon(pbtaj,col="grey80",density=100,border=NA)
lines(sjtaj,xlim=c(-.15,.3),col="firebrick2",lwd=3)
polygon(sjtaj,col="red",density=100,border=NA)
lines(bnptaj,xlim=c(-.15,.3),col="firebrick2",lwd=3)
polygon(bnptaj,col="lightpink",density=100,border=NA)
lines(sptaj,xlim=c(-.15,.3),col="cadetblue3",lwd=3)
polygon(sptaj,col="cadetblue3",density=100,border=NA)
lines(gbtaj,xlim=c(-.15,.3),col="cadetblue3",lwd=3)
polygon(gbtaj,col="cadetblue1",density=100,border=NA)

#####Plotting theta vs pi

par(mfrow=c(3,3))
plot(theta[,4],pi[,4],pch=20,cex=.5,col="black")
abline(a=0,b=1,col="red")
plot(theta[,5],pi[,5],pch=20,cex=.5,col="grey")
abline(a=0,b=1,col="red")
plot(theta[,6],pi[,6],pch=20,cex=.5,col="red")
abline(a=0,b=1,col="black")
plot(theta[,7],pi[,7],pch=20,cex=.5,col="orange")
abline(a=0,b=1,col="red")
plot(theta[,8],pi[,8],pch=20,cex=.5,col="yellow")
abline(a=0,b=1,col="red")
plot(theta[,9],pi[,9],pch=20,cex=.5,col="green")
abline(a=0,b=1,col="red")
plot(theta[,10],pi[,10],pch=20,cex=.5,col="blue")
abline(a=0,b=1,col="red")

