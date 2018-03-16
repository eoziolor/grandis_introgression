qual<-scan("~/Downloads/quality.txt")
hist(qual,breaks=1000)
subw<-qual>100
hist(qual[subw],breaks=1000)
length(qual)-length(qual[subw])
length(qual)

