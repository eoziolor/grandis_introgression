p<-read.table("~/analysis/data/fst/raw/BB.GB.fst.5kb1kb_p.bed", header=FALSE,stringsAsFactors = FALSE)
p[,4]<-as.numeric(p[,4])
names<-c("Scaf","start","end","p","count")
colnames(p)<-names

library(stringr)
library(dplyr)
library(gtools)
pc<-p %>% filter(str_detect(Scaf,"chr"))

subw<-pc[,5]>50

palette(c("black","grey","grey30","grey50"))
plot(-log(pc[subw,4]),pch=20,cex=.5,col=ifelse(p[subw,4]>0.05,factor(pc[subw,1]),"red"),
     ylab="-log p value",xlab="chr position",cex.lab=1.5,cex.axis=1.5)
pc1<-pc[subw,]
subs<-pc1[,4]<0.05
pcs<-pc1[subs,]

ord<-order(pcs[,4],decreasing=FALSE)
table(pcs[,1])

chr16<-pcs %>% filter(str_detect(Scaf,"chr16"))
