library(stringr)
library(dplyr)
library(gtools)
load("~/analysis/scripts/introgression/dxyscan.Rdata")
#install.packages("naturalsort")
library('naturalsort')

subw<-dxyscan[,"keep"]==TRUE
dxy<-dxyscan[subw,1:9]

ord<-mixedorder(dxy$chr)
dxy2<-dxy[ord,]
dxy3<-dxy2 %>% filter(str_detect(chr,"chr"))

col<-c()
for(i in 1:3){
  col[i]<-quantile(dxy3[,i+6],probs=0.005)
}


palette(c("grey50","grey70","black","grey30"))
par(mfrow=c(3,1),mar=c(2,2,0,0))
plot(dxy3[,"BB"],pch=20,cex=.5,col=ifelse(dxy3[,"BB"]<col[1],"red",as.factor(dxy3[,"chr"])),bty='n',
     xlab="",ylab="",xaxt='n')
plot(dxy3[,"VB"],pch=20,cex=.5,col=ifelse(dxy3[,"VB"]<col[2],"red",as.factor(dxy3[,"chr"])),bty='n',
     xlab="",ylab="",xaxt='n')
plot(dxy3[,"PB"],pch=20,cex=.5,col=ifelse(dxy3[,"PB"]<col[3],"red",as.factor(dxy3[,"chr"])),bty='n',
     xlab="",ylab="",xaxt='n')

     