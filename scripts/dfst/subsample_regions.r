source("http://bioconductor.org/biocLite.R")
biocLite()
library("rtracklayer")

bed1=import("~/analysis/data/dfst/PBS_keep_5kb.bed")
bedall=import("~/analysis/data/dfst/pbs_regions_sharedall.bed")
bed1overlall=bed1[bed1 %over% bedall]
hitsall<-findOverlaps(bedall,bed1)
allhit<-subjectHits(hitsall)
regions<-queryHits((hitsall))

orig<-read.table("~/analysis/data/dfst/zregions_split5kb.bed",header=F)

# p<-c()
# for(i in 1:69){
#   p[i]<-sum(regions==i)/length(regions)
# }
#
# pp<-c()
# for(i in 1:69){
#   pp<-c(pp,rep(p[i],sum(regions==i)))
# }
# pp<-unlist(pp)

v<-c()
for(i in 1:69){
  mini<-sample(x=which(regions==i),size=7)
  v<-c(v,mini)
}
v<-sort(v)

write.table(orig[v,],"~/analysis/data/dfst/zshared_10perregion.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)
