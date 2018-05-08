#r script to keep 50Mb into angsd files so that we can do SFS on SAF files with only 50Mb and not overload memory

p<-(sites[,3]-sites[2])/(sum(sites[,3]-sites[,2]))

p<-unlist(p)

v<-sample(x=1:253790,size=2500,prob=p)

v<-sort(v)

sum(sites[v,3]-sites[v,2])

write.table(sites[v,], file=(filename),row.names=FALSE,quotes=FALSE,col.names=FALSE)
