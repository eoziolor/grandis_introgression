library(stringr)
library(dplyr)
library(gtools)
library('naturalsort')
load("~/analysis/scripts/introgression/dxyscan.Rdata") #loading data into R


subw<-dxyscan[,"keep"]==TRUE
dxy<-dxyscan[subw,1:9]

ord<-mixedorder(dxy$chr)
dxyt<-dxy[ord,]

dxy2<-dxyt %>% filter(str_detect(chr,"chr"))

col<-c()
for(i in 1:3){
  col[i]<-quantile(dxy2[,i+6],probs=0.005)
}

chr16<-dxy2 %>% filter(str_detect(chr,"chr16"))#checking out chr16
ord<-order(chr16[,"BB"],decreasing=FALSE) #ordering them so I can remove the crappy scaffold
crap<-str_detect(chr16$scaf,"NW_012224891.1")
crap2<-str_detect(chr16$scaf,"NW_012224806.1")
ord<-order(chr16[!crap,"BB"],decreasing=FALSE) #ordering them so I can remove the crappy scaffold


crappy<-str_detect(dxy2$scaf,"NW_012224891.1")
dxy3<-dxy2[!crappy,]
crappy2<-str_detect(dxy3$scaf,"NW_012224806.1")
dxy4<-dxy3[!crappy2,]

dxy4$chr<-factor(dxy4$chr,levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                                     "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19",
                                     "chr20","chr21","chr22","chr23","chr24"))

#tiff("~/analysis/data/introgression/introgression.tiff",width=1200,height=600)
palette(c("grey40","grey80"))
par(mfrow=c(3,1),mar=c(2,4,2,0))
plot(dxy4[,"BB"],pch=20,cex=1,col=ifelse(dxy4[,"BB"]<col[2],"red",as.factor(dxy4[,"chr"])),bty='n',
     xlab="",ylab="",xaxt='n',ylim=c(-.0271,.021),cex.axis=2.5)
plot(dxy4[,"VB"],pch=20,cex=1,col=ifelse(dxy4[,"VB"]<col[1],"red",as.factor(dxy4[,"chr"])),bty='n',
     xlab="",ylab="",xaxt='n',ylim=c(-.0271,.021),cex.axis=2.5)
plot(dxy4[,"PB"],pch=20,cex=1,col=ifelse(dxy4[,"PB"]<col[3],"red",as.factor(dxy4[,"chr"])),bty='n',
     xlab="",ylab="",xaxt='n',ylim=c(-.0271,.021),cex.axis=2.5)
#dev.off()
     
