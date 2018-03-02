library(stringr)
library(dplyr)
library(gtools)
library('naturalsort')
load("~/analysis/scripts/introgression/dxyscan.Rdata") #loading data into R


subw<-dxyscan[,"keep"]==TRUE
dxy<-dxyscan[subw,1:9]

ord<-mixedorder(dxy$chr)
dxy2<-dxy[ord,]
dxy3<-dxy2 %>% filter(str_detect(chr,"chr"))

col<-c()
for(i in 1:3){
  col[i]<-quantile(dxy3[,i+6],probs=0.005)
}

chr16<-dxy3 %>% filter(str_detect(chr,"chr16"))#checking out chr16
ord<-order(chr16[,"BB"],decreasing=FALSE) #ordering them so I can remove the crappy scaffold
crap<-str_detect(chr16$scaf,"NW_012224891.1")
ord<-order(chr16[!crap,"BB"],decreasing=FALSE) #ordering them so I can remove the crappy scaffold


crappy<-str_detect(dxy3$scaf,"NW_012224891.1")
dxy4<-dxy3[!crappy,]
crappy2<-str_detect(dxy4$scaf,"NW_012224806.1")
dxy5<-dxy4[!crappy2,]

tiff("~/analysis/data/introgression/introgression.tiff",width=1200,height=600)
palette(c("grey50","grey70","black","grey30"))
par(mfrow=c(3,1),mar=c(2,4,2,0))
plot(dxy5[,"BB"],pch=20,cex=1,col=ifelse(dxy5[,"BB"]<col[2],"red",as.factor(dxy5[,"chr"])),bty='n',
     xlab="",ylab="",xaxt='n',ylim=c(-.0271,.021),cex.axis=2.5)
plot(dxy5[,"VB"],pch=20,cex=1,col=ifelse(dxy5[,"VB"]<col[1],"red",as.factor(dxy5[,"chr"])),bty='n',
     xlab="",ylab="",xaxt='n',ylim=c(-.0271,.021),cex.axis=2.5)
plot(dxy5[,"PB"],pch=20,cex=1,col=ifelse(dxy5[,"PB"]<col[3],"red",as.factor(dxy5[,"chr"])),bty='n',
     xlab="",ylab="",xaxt='n',ylim=c(-.0271,.021),cex.axis=2.5)
dev.off()
     
