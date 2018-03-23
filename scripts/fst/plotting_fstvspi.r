library('ggplot2')
library("RColorBrewer")
require('gridExtra')

#loading fst and finding representative windows----
fs <- list.files("~/analysis/data/fst/raw/", "*5kb1kb.bed",full.names=TRUE)

fst <- list()

for (i in 1:21){
  fst[[i]] <- read.table(fs[i],stringsAsFactors=FALSE)
  fst[[i]][,4] <- as.numeric(fst[[i]][,4])
}

nfs <- gsub(".*\\/","",fs)
nfs <- gsub(".fst.*","",nfs)
names(fst)<-nfs

nsnps <-fst[[1]][,5]

for (i in 2:21){
  
  nsnps <- nsnps + fst[[i]][,5]
}

nsnps <- nsnps/21

subx <- nsnps > 20

#loading pi and theta and finding their rep windows----
theta<-read.table("~/analysis/data/angsd/thetas_neut_5kb",header=TRUE, sep=',')
pi<-read.table("~/analysis/data/angsd/pi_neut_5kb", header=TRUE, sep=',')

par(mfrow=c(2,5))

subw<-pi[,"keep"]>0
suby<-as.numeric(subw)+as.numeric(subx)
subz<-suby>1

##tiff(filename="~/share/VCF/angsd_res/theta_graphs/Pi_fst.tiff",width=1200,height=700,units="px")

#plotting them all in ggplot
a<-ggplot(pi[subz,],
       aes(x=bb,y=gb)) +
  geom_point(size=fst[["BB.GB"]][subz,4],
    aes(colour=fst[["BB.GB"]][subz,4])) +
  scale_color_gradient(low="yellow",high="black") +
  theme_classic()+
  scale_x_continuous(limits=c(0,0.1))+
  scale_y_continuous(limits=c(0,0.1))+
  geom_abline(slope=1,lty=2)

b<-ggplot(pi[subz,],
       aes(x=vb,y=gb)) +
  geom_point(size=fst[["VB.GB"]][subz,4],
             aes(colour=fst[["VB.GB"]][subz,4])) +
  scale_color_gradient(low="yellow",high="black") +
  theme_classic()+
  scale_x_continuous(limits=c(0,0.1))+
  scale_y_continuous(limits=c(0,0.1))+
  geom_abline(slope=1,lty=2)

c<-ggplot(pi[subz,],
          aes(x=pb,y=gb)) +
  geom_point(size=fst[["PB.GB"]][subz,4],
             aes(colour=fst[["PB.GB"]][subz,4])) +
  scale_color_gradient(low="yellow",high="black") +
  theme_classic()+
  scale_x_continuous(limits=c(0,0.1))+
  scale_y_continuous(limits=c(0,0.1))+
  geom_abline(slope=1,lty=2)

d<-ggplot(pi[subz,],
          aes(x=sj,y=gb)) +
  geom_point(size=fst[["SJ.GB"]][subz,4],
             aes(colour=fst[["SJ.GB"]][subz,4])) +
  scale_color_gradient(low="yellow",high="black") +
  theme_classic()+
  scale_x_continuous(limits=c(0,0.1))+
  scale_y_continuous(limits=c(0,0.1))+
  geom_abline(slope=1,lty=2)

e<-ggplot(pi[subz,],
       aes(x=bnp,y=gb)) +
  geom_point(size=fst[["BNP.GB"]][subz,4],
             aes(colour=fst[["BNP.GB"]][subz,4])) +
  scale_color_gradient(low="yellow",high="black") +
  theme_classic()+
  scale_x_continuous(limits=c(0,0.1))+
  scale_y_continuous(limits=c(0,0.1))+
  geom_abline(slope=1,lty=2)

f<-ggplot(pi[subz,],
             aes(x=sp,y=gb)) +
  geom_point(size=fst[["GB.SP"]][subz,4],
             aes(colour=fst[["GB.SP"]][subz,4])) +
  scale_color_gradient(low="yellow",high="black") +
  theme_classic()+
  scale_x_continuous(limits=c(0,0.1))+
  scale_y_continuous(limits=c(0,0.1))+
  geom_abline(slope=1,lty=2)

a2<-ggplot(pi[subz,],
          aes(x=bb,y=sp)) +
  geom_point(size=fst[["BB.SP"]][subz,4],
             aes(colour=fst[["BB.SP"]][subz,4])) +
  scale_color_gradient(low="yellow",high="black") +
  theme_classic()+
  scale_x_continuous(limits=c(0,0.1))+
  scale_y_continuous(limits=c(0,0.1))+
  geom_abline(slope=1,lty=2)

b2<-ggplot(pi[subz,],
          aes(x=vb,y=sp)) +
  geom_point(size=fst[["VB.SP"]][subz,4],
             aes(colour=fst[["VB.SP"]][subz,4])) +
  scale_color_gradient(low="yellow",high="black") +
  theme_classic()+
  scale_x_continuous(limits=c(0,0.1))+
  scale_y_continuous(limits=c(0,0.1))+
  geom_abline(slope=1,lty=2)

c2<-ggplot(pi[subz,],
          aes(x=pb,y=sp)) +
  geom_point(size=fst[["PB.SP"]][subz,4],
             aes(colour=fst[["PB.SP"]][subz,4])) +
  scale_color_gradient(low="yellow",high="black") +
  theme_classic()+
  scale_x_continuous(limits=c(0,0.1))+
  scale_y_continuous(limits=c(0,0.1))+
  geom_abline(slope=1,lty=2)

d2<-ggplot(pi[subz,],
          aes(x=sj,y=sp)) +
  geom_point(size=fst[["SJ.SP"]][subz,4],
             aes(colour=fst[["SJ.SP"]][subz,4])) +
  scale_color_gradient(low="yellow",high="black") +
  theme_classic()+
  scale_x_continuous(limits=c(0,0.1))+
  scale_y_continuous(limits=c(0,0.1))+
  geom_abline(slope=1,lty=2)

e2<-ggplot(pi[subz,],
          aes(x=bnp,y=sp)) +
  geom_point(size=fst[["BNP.SP"]][subz,4],
             aes(colour=fst[["BNP.SP"]][subz,4])) +
  scale_color_gradient(low="yellow",high="black") +
  theme_classic()+
  scale_x_continuous(limits=c(0,0.1))+
  scale_y_continuous(limits=c(0,0.1))+
  geom_abline(slope=1,lty=2)

grid.arrange(grobs=list(a,b,c,d,e,f,a2,b2,c2,d2,e2),ncol=6)

#plotting the same with base r----
######plot GB vs VB for pi and fst
plot(pi[,5],pi[,10],pch=20,cex=0.8,xlab="Pi Vince Bayou",ylab="Pi Gangs Bayou",
     xlim=c(0,.15),ylim=c(0,.15),col=(ifelse(fst[["VB.GB"]][,4]<0.1,"black",
	ifelse(fst[["VB.GB"]][,4]>0.1 & fst[["VB.GB"]][,4]<0.2,"orange","red"))))
legend(c(500,1000),c(0,180), c("fst<0.1","0.1<fst<0.2","fst>0.2"),col=c("black","orange","red"),pch=20,cex=.8)

quantile(as.numeric(fst[["VB.GB"]][,4]),na.rm=TRUE)
#mean fst is 0.014

###plot GB vs PB for pi and fst
plot(pi[,6],pi[,10],pch=20,cex=0.8,xlab="Pi Patrick Bayou",ylab="Pi Gangs Bayou",
     xlim=c(0,.15),ylim=c(0,.15),col=(ifelse(fst[["PB.GB"]][,4]<0.1,"black",
	ifelse(fst[["PB.GB"]][,4]>0.1 & fst[["PB.GB"]][,4]<0.2,"orange","red"))))
legend(c(500,1000),c(0,180), c("fst<0.1","0.1<fst<0.2","fst>0.2"),col=c("black","orange","red"),pch=20,cex=.8)

quantile(as.numeric(fst[["PB.GB"]][,4]),na.rm=TRUE)
#mean fst is 0.011

###plot GB vs SJSP for pi and fst
plot(pi[,7],pi[,10],pch=20,cex=0.8,xlab="Pi San Jacinto State Park",ylab="Pi Gangs Bayou",
     xlim=c(0,.15),ylim=c(0,.15),col=(ifelse(fst[["SJ.GB"]][,4]<0.1,"black",
	ifelse(fst[["SJ.GB"]][,4]>0.1 & fst[["SJ.GB"]][,4]<0.2,"orange","red"))))
legend(c(500,1000),c(0,180), c("fst<0.1","0.1<fst<0.2","fst>0.2"),col=c("black","orange","red"),pch=20,cex=.8)

quantile(as.numeric(fst[["SJ.GB"]][,4]),na.rm=TRUE)
#mean fst is 0.05

###Plot GB vs BNP for pi and fst
plot(pi[,8],pi[,10],pch=20,cex=0.8,xlab="Pi Baytown Nature Park",ylab="Pi Gangs Bayou",
     xlim=c(0,.15),ylim=c(0,.15),col=(ifelse(fst[["BNP.GB"]][,4]<0.1,"black",
	ifelse(fst[["BNP.GB"]][,4]>0.1 & fst[["BNP.GB"]][,4]<0.2,"orange","red"))))
legend(c(500,1000),c(0,180), c("fst<0.1","0.1<fst<0.2","fst>0.2"),col=c("black","orange","red"),pch=20,cex=.8)

quantile(as.numeric(fst[["BNP.GB"]][,4]),na.rm=TRUE)
#mean fst is 0.005

#SP

#####plot of SP vs BB for pi and fst
plot(pi[,4],pi[,9],pch=20,cex=0.8,xlab="Pi Buffalo Bayou",ylab="Pi Smith Point",
     xlim=c(0,.15),ylim=c(0,.15),col=(ifelse(fst[["BB.SP"]][,4]<0.1,"black",
	ifelse(fst[["BB.SP"]][,4]>0.1 & fst[["BB.SP"]][,4]<0.2,"orange","red"))))
legend(c(500,1000),c(0,180), c("fst<0.1","0.1<fst<0.2","fst>0.2"),col=c("black","orange","red"),pch=20,cex=.8)

quantile(as.numeric(fst[["BB.SP"]][,4]),na.rm=TRUE)
#mean fst is 0.037

######plot SP vs VB for pi and fst
plot(pi[,5],pi[,9],pch=20,cex=0.8,xlab="Pi Vince Bayou",ylab="Pi Smith Point",
     xlim=c(0,.15),ylim=c(0,.15),col=(ifelse(fst[["VB.SP"]][,4]<0.1,"black",
	ifelse(fst[["VB.SP"]][,4]>0.1 & fst[["VB.SP"]][,4]<0.2,"orange","red"))))
legend(c(500,1000),c(0,180), c("fst<0.1","0.1<fst<0.2","fst>0.2"),col=c("black","orange","red"),pch=20,cex=.8)

quantile(as.numeric(fst[["VB.SP"]][,4]),na.rm=TRUE)
#mean fst is 0.010

###plot SP vs PB for pi and fst
plot(pi[,6],pi[,9],pch=20,cex=0.8,xlab="Pi Patrick Bayou",ylab="Pi Smith Point",
     xlim=c(0,.15),ylim=c(0,.15),col=(ifelse(fst[["PB.SP"]][,4]<0.1,"black",
	ifelse(fst[["PB.SP"]][,4]>0.1 & fst[["PB.SP"]][,4]<0.2,"orange","red"))))
legend(c(500,1000),c(0,180), c("fst<0.1","0.1<fst<0.2","fst>0.2"),col=c("black","orange","red"),pch=20,cex=.8)

quantile(as.numeric(fst[["PB.SP"]][,4]),na.rm=TRUE)
#mean fst is 0.011

###plot SP vs SJSP for pi and fst
plot(pi[,7],pi[,9],pch=20,cex=0.8,xlab="Pi San Jacinto State Park",ylab="Pi Smith Point",
     xlim=c(0,.15),ylim=c(0,.15),col=(ifelse(fst[["SJ.SP"]][,4]<0.1,"black",
	ifelse(fst[["SJ.SP"]][,4]>0.1 & fst[["SJ.SP"]][,4]<0.2,"orange","red"))))
legend(c(500,1000),c(0,180), c("fst<0.1","0.1<fst<0.2","fst>0.2"),col=c("black","orange","red"),pch=20,cex=.8)

quantile(as.numeric(fst[["SJ.SP"]][,4]),na.rm=TRUE)
#mean fst is 0.05

###Plot SP vs BNP for pi and fst
plot(pi[,8],pi[,9],pch=20,cex=0.8,xlab="Pi Baytown Nature Park",ylab="Pi Smith Point",
     xlim=c(0,.15),ylim=c(0,.15),col=(ifelse(fst[["BNP.SP"]][,4]<0.1,"black",
	ifelse(fst[["BNP.SP"]][,4]>0.1 & fst[["BNP.SP"]][,4]<0.2,"orange","red"))))
legend(c(500,1000),c(0,180), c("fst<0.1","0.1<fst<0.2","fst>0.2"),col=c("black","orange","red"),pch=20,cex=.8)

quantile(as.numeric(fst[["SP.BNP"]][,4]),na.rm=TRUE)
#mean fst is 0.003

###GB vs SP

plot(pi[,9],pi[,10],pch=20,cex=0.8,xlab="Pi Smith Point",ylab="Pi Gangs Bayou",
     xlim=c(0,.15),ylim=c(0,.15),col=(ifelse(fst[["GB.SP"]][,4]<0.1,"black",
	ifelse(fst[["GB.SP"]][,4]>0.1 & fst[["GB.SP"]][,4]<0.2,"orange","red"))))
legend(c(500,1000),c(0,180), c("fst<0.1","0.1<fst<0.2","fst>0.2"),col=c("black","orange","red"),pch=20,cex=.8)

quantile(as.numeric(fst[["GB.SP"]][,4]),na.rm=TRUE)
#mean fst is 0.002

