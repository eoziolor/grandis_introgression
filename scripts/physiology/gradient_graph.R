grad<-read.table("~/analysis/scripts/physiology/ec50_chem.csv",header=TRUE,sep=',')


###Plotting EC,ER,dist chem

plot(grad[,"Max_EROD"],log(grad[,"ec50"]),pch=16,cex=(log((grad[,'contam']/10)+1)*6),
     ylab="Level of resistance", xlab="Maximum CYP1A activity (as % of highest reference)",
     col=(ifelse(grad[,"ec50"]<.01,"cadetblue3",ifelse(grad[,"ec50"]<0.04,"gold2",
     ifelse(grad[,"ec50"]<1,"firebrick2","black")))),ylim=c(-7,9),bty="l",
     xlim=c(10,100),cex.lab=1.5)

legend(x=c(90,100),y=c(4,9),c("1","20","100"),col="grey",
       pch=16,bty="n",pt.cex=c(0.5718611,6.5916737,14.3873716), x.intersp=1.3,y.intersp=c(1.25,1,1.25))

legend(x=c(85.5,100),y=c(10,9),c("Relative"),
       bty="n")

legend(x=c(83.5,100),y=c(8.5,9.5),c("contamination"),
       bty="n")

legend(x=c(88,102),y=c(9,1.2),box.col="grey",bg=NA,c(""))

grad[order(grad[,"Max_EROD"],decreasing = TRUE),]

legend(x=c(6.4,7.2),y=c(10.2,11.3),c("BB"),col="black",bty ="n",cex=1.2)
legend(x=c(10.8,12.3),y=c(7,8.6),c("PB"),col="black",bty ="n",cex=1.2)
legend(x=c(26.7,28),y=c(7,8.6),c("VB"),col="black",bty ="n",cex=1.2)
legend(x=c(38.5,40),y=c(2,2.3),c("SJ"),text.col="firebrick2",bty ="n",cex=1.2)
legend(x=c(56.1,57.1),y=c(.2,.6),c("CB"),text.col="firebrick2",bty ="n",cex=1.2)
legend(x=c(61.6,64.5),y=c(0.7,1.7),c("BNP"),text.col="firebrick2",bty ="n",cex=1.2)
legend(x=c(66.2,68),y=c(-1.2,-.7),c("PG"),text.col="gold2",bty ="n",cex=1.2)
legend(x=c(75.7,77),y=c(-4,-3),c("HP"),text.col="gold2",bty ="n",cex=1.2)
legend(x=c(76.2,78),y=c(-1.3,-.4),c("FP"),text.col="gold2",bty ="n",cex=1.2)
legend(x=c(81.5,83.5),y=c(-4.5,-3.5),c("GB"),text.col="cadetblue3",bty ="n",cex=1.2)
legend(x=c(85.5,88.5),y=c(-3,-2),c("FB"),text.col="cadetblue3",bty ="n",cex=1.2)
legend(x=c(90.5,93),y=c(-4,-3),c("SP"),text.col="cadetblue3",bty ="n",cex=1.2)

#legend(x=c(90,100),y=c(8,-1),c("SP","GB","FB","HP","CB","PG","FP","BNP","SJSP","PB","VB","BB"),
#       col=c("deepskyblue","deepskyblue","deepskyblue","orange","red","orange",
#             "orange","red","red","black","black","black"),pch=20,pt.cex=1.5,
#       y.intersp=.6,bty="n")

#legend(x=c(8500,11000),y=c(9,8),c("Populations"),bty="n")


# ###OLDER GRAPH
# plot(-grad[,"dist"],log(grad[,"ec50"]),pch=16,cex=(log(grad[,'contam'])*2),
#      ylab="Log EC50 (mg/L)", xlab="distance from most inland population (km)",
#      col=(ifelse(grad[,"ec50"]<.01,"deepskyblue",ifelse(grad[,"ec50"]<0.04,"orange",
#                                                         ifelse(grad[,"ec50"]<1,"red","black")))),ylim=c(-7,8),bty="l")
# 
# abline(lty=5,col="grey",h=-4.61)
# 
# legend(x=c(-155,-130),y=c(7,-2),c("SP","GB","FB","HP","CB","PG","FP","BNP","SJSP","PB","VB","BB"),
#        col=c("deepskyblue","deepskyblue","deepskyblue","orange","red","orange",
#              "orange","red","red","black","black","black"),pch=20,pt.cex=1.5,
#        y.intersp=.4,bty="n")
# 
# legend(x=c(-125,-100),y=c(7.5,0),c("1","20","100"),col="grey",
#        pch=16,bty="n",pt.cex=c(.19,6,9.21), x.intersp=1.5,y.intersp=1.5)
# 
# legend(x=c(-160,-130),y=c(9,8),c("Populations"),bty="n")
# legend(x=c(-145,-110),y=c(9,8),c("Natural log of relative contamination"),
#        bty="n")
# legend(x=c(-155,-95),y=c(8.5,-3),box.col="grey",bg=NA,c(""))
