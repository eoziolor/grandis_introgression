setwd("~/analysis/data/admixture/")
cov<-read.table("~/analysis/data/angsd/cov_10Mbrand.txt.gz",header=F)
names<-c("chrom","pos","cov")

hist(cov$cov,breaks=1000)

subw<-cov$cov<300
hist(cov[subw,"cov"],breaks=1000)
summary(cov$cov)
summary(cov[subw,"cov"])
