orig<-read.table("~/analysis/data/dfst/zsharedregions_split5kb.bed",header=F)
library("rtracklayer")
bed1=import("~/analysis/data/dfst/PBS_keep_5kb.bed")
bedall=import("~/analysis/data/dfst/outlier_regions/pbs_regions_sharedall.bed")
bed1overlall=bed1[bed1 %over% bedall]
hitsall<-findOverlaps(bedall,bed1)
allhit<-subjectHits(hitsall)
regions<-queryHits((hitsall))

v<-c()
for(i in 1:69){
  mini<-sample(x=which(regions==i),size=7)
  v<-c(v,mini)
}
v<-sort(v)

write.table(orig[v,],"~/analysis/data/dfst/zshared_10perregion.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)


###grabbing all regions except chr1

no1<-grep("chr1\\b|chr10\\b",orig[,1],invert=TRUE)
orig2<-orig[no1,]
regions2<-regions[-(grep("chr1\\b",orig[,1]))]

v<-c()
for(i in 7:69){
  mini<-sample(x=which(regions2==i),size=7)
  v<-c(v,mini)
}
v<-sort(v)

write.table(orig2[v,],"~/analysis/data/dfst/outlier_regions/zshared_10perregion_no1or10.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)

###Doing the same for random regions
orig<-read.table("~/analysis/data/dfst/random/5kb1kb.bed",header=F)
v<-sample(x=length(orig[,1]),size=483)
v<-sort(v)
write.table(orig[v,],"~/analysis/data/dfst/random/random_483win.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)


#tried it other ways
# source("http://bioconductor.org/biocLite.R")
# biocLite()


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
