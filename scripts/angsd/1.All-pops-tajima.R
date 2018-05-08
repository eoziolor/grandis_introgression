###Buffalo Bayou loding pi(column6), theta(column4), and counts(columns 5,7)


bb<-read.table("~/analysis/data/angsd/subsample/BB_neut_1kb.bed",stringsAsFactors=FALSE)
vb<-read.table("~/analysis/data/angsd/subsample/VB_neut_1kb.bed",stringsAsFactors=FALSE)
pb<-read.table("~/analysis/data/angsd/subsample/PB_neut_1kb.bed",stringsAsFactors=FALSE)
sj<-read.table("~/analysis/data/angsd/subsample/SJ_neut_1kb.bed",stringsAsFactors=FALSE)
bnp<-read.table("~/analysis/data/angsd/subsample/BNP_neut_1kb.bed",stringsAsFactors=FALSE)
sp<-read.table("~/analysis/data/angsd/subsample/SP_neut_1kb.bed",stringsAsFactors=FALSE)
gb<-read.table("~/analysis/data/angsd/subsample/GB_neut_1kb.bed",stringsAsFactors=FALSE)

colnam<-c("scaf","Start","End","Theta","Tcount","Pi","Pcount")
colnames(bb)<-colnam
colnames(vb)<-colnam
colnames(pb)<-colnam
colnames(sj)<-colnam
colnames(bnp)<-colnam
colnames(sp)<-colnam
colnames(gb)<-colnam

pops<-list(bb,vb,pb,sj,bnp,sp,gb)
popnames<-c("bb","vb","pb","sj","bnp","sp","gb")
names(pops)<-popnames

for(i in popnames){
  for(j in 4:7){
    pops[[i]][,j]<-as.numeric(pops[[i]][,j])
  }
}

popbase<-list()

for(i in popnames){
  popbase[[i]]<-cbind(pops[[i]][,1:3],pops[[i]][,4]/pops[[i]][,5],pops[[i]][,6]/pops[[i]][,7])
}

names<-c("scaf", "start","end",'theta/b',"pi/b")

for(i in popnames){
  colnames(popbase[[i]])<-names
}

source("~/analysis/scripts/angsd/tajimas.r")

taj<-list()

for(i in popnames){
  for(j in 1:1027369){
    taj[[i]]<-c(taj[[i]],tajimas(popbase[[i]][j,5],popbase[[i]][j,4],24))
  }
}


#OLD WAY
# ###Vince Bayou
# 
# vb<-read.table("~/analysis/data/angsd/subsample/VB_neut_1kb.bed",stringsAsFactors=FALSE)
# colnames(vb)<-colnam
# 
# vb[,4]<-as.numeric(vb[,4])
# vb[,5]<-as.numeric(vb[,5])
# vb[,6]<-as.numeric(vb[,6])
# vb[,7]<-as.numeric(vb[,7])
# 
# vbbase<-cbind(vb[,1:3],vb[,4]/vb[,5],vb[,6]/vb[,7])
# names<-c("scaf", "start","end",'theta/b',"pi/b")
# colnames(vbbase)<-names
# 
# vbtaj <- c ()
# 
# for (i in 1:1027369){
#   vbtaj <- c(vbtaj,tajimas(vbbase[i,5],vbbase[i,4],24))
# }
# 
# # plot(vbtaj,pch=20,cex=0.1)
# # 
# # write.table(vbtaj,"~/analysis/angsd/vbtaj",
# #             row.names = FALSE,col.names = FALSE,quote = FALSE)



############################ Calculating coverage

cov<-cbind(bb[1:3],bb[,7],vb[,7],pb[,7],sj[,7],bnp[,7],sp[,7],gb[,7])

nsnps<-cov[,4]
for (i in 5:10){
  nsnps <- nsnps + cov[,i]
}
nsnps <- nsnps/7

subw <- nsnps > 20

###

taj<-cbind(bbbase[1:3],taj[["bb"]],taj[["vb"]],taj[["pb"]],taj[["sj"]],taj[["bnp"]],taj[["sp"]],taj[["gb"]])

taj<-cbind(taj,keep=as.numeric(subw))

tajname<-c("scaf","start","end","bb","vb","pb","sjsp","bnp",
            "sp","gb","keep")
colnames(taj)<-tajname

write.csv(taj,file="~/analysis/data/angsd/taj_sub",quote=FALSE,row.names=FALSE)

theta<-cbind(popbase[["bb"]][1:3],popbase[["bb"]][4],popbase[["vb"]][4],popbase[["pb"]][4],popbase[["sj"]][4],popbase[["bnp"]][4],
             popbase[["sp"]][4],popbase[["gb"]][4],keep=as.numeric(subw))

thetname<-c("scaf","start","end","bb","vb","pb","sjsp","bnp",
            "sp","gb","keep")
colnames(theta)<-thetname

write.csv(theta,file="~/analysis/data/angsd/thetas_sub",quote=FALSE,row.names=FALSE)


pi<-cbind(popbase[["bb"]][1:3],popbase[["bb"]][5],popbase[["vb"]][5],popbase[["pb"]][5],popbase[["sj"]][5],popbase[["bnp"]][5],
          popbase[["sp"]][5],popbase[["gb"]][5],keep=as.numeric(subw))

piname<-c("scaf","start","end","bb","vb","pb","sj","bnp","sp","gb","keep")
colnames(pi)<-piname
write.csv(pi,file="~/analysis/data/angsd/pi_sub",quote=FALSE,row.names=FALSE)

########



