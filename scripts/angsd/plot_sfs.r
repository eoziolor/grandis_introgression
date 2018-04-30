sf<-list.files("~/analysis/data/angsd/raw/","*.sfs",full.names=TRUE)
cols<-c("black","grey40","grey80","firebrick2","lightpink","cadetblue1","cadetblue3")
pop<-list("bb","vb","pb","sj","bnp","sp","gb")

for(i in 1:7){
  pop[[i]]<-scan(sf[[i]])
}

par(mfrow=c(3,3),mar=c(2,2,2,2))
for(i in 1:7){
  plot(log(pop[[i]]),col=cols[i],pch=20,lwd=3)
}

# points(log(pop[[1]]),col="black",xlim=c(0,49),pch=20,cex=2,bty='n')
# cols<-c("black","grey40","grey80","firebrick2","lightpink","cadetblue1","cadetblue3")
# 
# for(i in 1:7){
#   lines(log(pop[[i]]),col=cols[i],pch=20,lwd=3)
# }