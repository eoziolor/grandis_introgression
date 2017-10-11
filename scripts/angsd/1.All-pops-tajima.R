#Alternative Script

fs <- list.files("~/analysis/data/angsd/raw/", "*5kb1kb.bed",full.names=TRUE)
neut <- list()

for (i in 1:7){
  neut[[i]] <- read.table(fs[i],stringsAsFactors=FALSE)
  neut[[i]][,4] <- as.numeric(neut[[i]][,4])
  neut[[i]][,5] <- as.numeric(neut[[i]][,5])
  neut[[i]][,6] <- as.numeric(neut[[i]][,6])
  neut[[i]][,7] <- as.numeric(neut[[i]][,7])
}

nfs <- gsub(".*\\/","",fs)
nfs <- gsub("_neut_*","",nfs)
nfs <- gsub("5kb1kb.bed*","",nfs)
names(neut)<-nfs

source("~/analysis/scripts/angsd/tajimas.r")

#LEFT OFF HERE

###Buffalo Bayou loding pi(column6), theta(column4), and counts(columns 5,7)


bb<-read.table("~/analysis/data/angsd/raw/BB_neut_5kb1kb.bed",stringsAsFactors=FALSE)
colnam<-c("scaf","Start","End","Theta","Tcount","Pi","Pcount")
colnames(bb)<-colnam

bb[,4]<-as.numeric(bb[,4])
bb[,5]<-as.numeric(bb[,5])
bb[,6]<-as.numeric(bb[,6])
bb[,7]<-as.numeric(bb[,7])

bbbase<-cbind(bb[,1:3],bb[,4]/bb[,5],bb[,6]/bb[,7])
names<-c("scaf", "start","end",'theta/b',"pi/b")
colnames(bbbase)<-names

source("~/analysis/angsd/scripts/tajimas.r")

bbtaj <- c ()

for (i in 1:1026857){
	bbtaj <- c(bbtaj,tajimas(bbbase[i,5],bbbase[i,4],24))
	}

plot(bbtaj,pch=20,cex=0.1)
quantile(bbtaj,na.rm=TRUE)

write.table(bbtaj,"~/analysis/angsd/bbtaj",
            row.names = FALSE,col.names = FALSE,quote = FALSE)

###Vince Bayou

vb<-read.table("~/analysis/data/angsd/raw/VB_neut_5kb1kb.bed",stringsAsFactors=FALSE)
colnames(vb)<-colnam

vb[,4]<-as.numeric(vb[,4])
vb[,5]<-as.numeric(vb[,5])
vb[,6]<-as.numeric(vb[,6])
vb[,7]<-as.numeric(vb[,7])

vbbase<-cbind(vb[,1:3],vb[,4]/vb[,5],vb[,6]/vb[,7])
names<-c("scaf", "start","end",'theta/b',"pi/b")
colnames(vbbase)<-names

source("~/analysis/angsd/scripts/tajimas.r")

vbtaj <- c ()

for (i in 1:1026857){
  vbtaj <- c(vbtaj,tajimas(vbbase[i,5],vbbase[i,4],49))
}

plot(vbtaj,pch=20,cex=0.1)

write.table(vbtaj,"~/analysis/angsd/vbtaj",
            row.names = FALSE,col.names = FALSE,quote = FALSE)

###Patrick Bayou

pb<-read.table("~/analysis/data/angsd/raw/PB_neut_5kb1kb.bed",stringsAsFactors=FALSE)
colnames(pb)<-colnam

pb[,4]<-as.numeric(pb[,4])
pb[,5]<-as.numeric(pb[,5])
pb[,6]<-as.numeric(pb[,6])
pb[,7]<-as.numeric(pb[,7])

pbbase<-cbind(pb[,1:3],pb[,4]/pb[,5],pb[,6]/pb[,7])
names<-c("scaf", "start","end",'theta/b',"pi/b")
colnames(pbbase)<-names

source("~/analysis/angsd/scripts/tajimas.r")

pbtaj <- c ()

for (i in 1:1026857){
  pbtaj <- c(pbtaj,tajimas(pbbase[i,5],pbbase[i,4],47))
}

plot(pbtaj,pch=20,cex=0.1)

write.table(pbtaj,"~/analysis/angsd/pbtaj",
            row.names = FALSE,col.names = FALSE,quote = FALSE)


###San Jacinto State Park

sj<-read.table("~/analysis/data/angsd/raw/SJ_neut_5kb1kb.bed",stringsAsFactors=FALSE)
colnames(sj)<-colnam

sj[,4]<-as.numeric(sj[,4])
sj[,5]<-as.numeric(sj[,5])
sj[,6]<-as.numeric(sj[,6])
sj[,7]<-as.numeric(sj[,7])

sjbase<-cbind(sj[,1:3],sj[,4]/sj[,5],sj[,6]/sj[,7])
names<-c("scaf", "start","end",'theta/b',"pi/b")
colnames(sjbase)<-names

source("~/analysis/angsd/scripts/tajimas.r")

sjtaj <- c ()

for (i in 1:1026857){
  sjtaj <- c(sjtaj,tajimas(sjbase[i,5],sjbase[i,4],24))
}

plot(sjtaj,pch=20,cex=0.1)

write.table(sjtaj,"~/analysis/angsd/sjtaj",row.names = FALSE,col.names = FALSE,quote = FALSE)


###Baytown Nature Park

bnp<-read.table("~/analysis/data/angsd/raw/BNP_neut_5kb1kb.bed",stringsAsFactors=FALSE)
colnames(bnp)<-colnam

bnp[,4]<-as.numeric(bnp[,4])
bnp[,5]<-as.numeric(bnp[,5])
bnp[,6]<-as.numeric(bnp[,6])
bnp[,7]<-as.numeric(bnp[,7])

bnpbase<-cbind(bnp[,1:3],bnp[,4]/bnp[,5],bnp[,6]/bnp[,7])
names<-c("scaf", "start","end",'theta/b',"pi/b")
colnames(bnpbase)<-names

source("~/analysis/angsd/scripts/tajimas.r")

bnptaj <- c ()

for (i in 1:1026857){
  bnptaj <- c(bnptaj,tajimas(bnpbase[i,5],bnpbase[i,4],48))
}

plot(bnptaj,pch=20,cex=0.1)

write.table(bnptaj,"~/analysis/angsd/bnptaj",
            row.names = FALSE,col.names = FALSE,quote = FALSE)


###Smith Point

sp<-read.table("~/analysis/data/angsd/raw/SP_neut_5kb1kb.bed",stringsAsFactors=FALSE)
colnames(sp)<-colnam

sp[,4]<-as.numeric(sp[,4])
sp[,5]<-as.numeric(sp[,5])
sp[,6]<-as.numeric(sp[,6])
sp[,7]<-as.numeric(sp[,7])

spbase<-cbind(sp[,1:3],sp[,4]/sp[,5],sp[,6]/sp[,7])
names<-c("scaf", "start","end",'theta/b',"pi/b")
colnames(spbase)<-names

source("~/analysis/angsd/scripts/tajimas.r")

sptaj <- c ()

for (i in 1:1026857){
  sptaj <- c(sptaj,tajimas(spbase[i,5],spbase[i,4],48))
}

plot(sptaj,pch=20,cex=0.1)

write.table(sptaj,"~/analysis/angsd/sptaj",
            row.names = FALSE,col.names = FALSE,quote = FALSE)


###Gangs Bayou

gb<-read.table("~/analysis/data/angsd/raw/GB_neut_5kb1kb.bed",stringsAsFactors=FALSE)
colnames(gb)<-colnam

gb[,4]<-as.numeric(gb[,4])
gb[,5]<-as.numeric(gb[,5])
gb[,6]<-as.numeric(gb[,6])
gb[,7]<-as.numeric(gb[,7])

gbbase<-cbind(gb[,1:3],gb[,4]/gb[,5],gb[,6]/gb[,7])
names<-c("scaf", "start","end",'theta/b',"pi/b")
colnames(gbbase)<-names

source("~/analysis/angsd/scripts/tajimas.r")

gbtaj <- c ()

for (i in 1:1026857){
  gbtaj <- c(gbtaj,tajimas(gbbase[i,5],gbbase[i,4],48))
}

plot(gbtaj,pch=20,cex=0.1)

write.table(gbtaj,"~/analysis/angsd/gbtaj",
            row.names = FALSE,col.names = FALSE,quote = FALSE)

############################ Calculating coverage

cov<-cbind(bb[1:3],bb[,7],vb[,7],pb[,7],sj[,7],bnp[,7],sp[,7],gb[,7])

nsnps<-cov[,4]
for (i in 5:10){
  nsnps <- nsnps + cov[,i]
}
nsnps <- nsnps/7

subw <- nsnps > 20

###

taj<-cbind(bbbase[1:3],bbtaj,vbtaj,pbtaj,sjtaj,bnptaj,sptaj,gbtaj)
tajname<-c("scaf","start","end","bb","vb","pb","sjsp","bnp",
            "sp","gb")
colnames(taj)<-tajname
write.csv(taj,file="~/analysis/angsd/taj",quote=FALSE,row.names=FALSE)

theta<-cbind(bbbase[1:3],bbbase[4],vbbase[4],pbbase[4],sjbase[4],bnpbase[4]
             ,spbase[4],gbbase[4])
thetname<-c("scaf","start","end","bb","vb","pb","sjsp","bnp",
            "sp","gb")
colnames(theta)<-thetname
write.csv(theta,file="~/analysis/angsd/thetas_neut_5kb",quote=FALSE,row.names=FALSE)


pi<-cbind(bb[1:3],bbbase[5],vbbase[5],pbbase[5],sjbase[5],bnpbase[5],
          spbase[5],gbbase[5],keep=as.numeric(subw))
piname<-c("scaf","start","end","bb","vb","pb","sj","bnp","sp","gb","keep")
colnames(pi)<-piname
write.csv(pi,file="~/analysis/data/angsd/pi_neut_5kb",quote=FALSE,row.names=FALSE)

########



