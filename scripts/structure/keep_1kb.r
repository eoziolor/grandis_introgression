#!/usr/bin/env Rscript

orig<-read.table("~/analysis/data/structure/top50.snps", header=F)

sub1k<-orig[sample(nrow(orig),1000),]

write.table(sub1k, file="~/analysis/data/structure/top50.1ksnps.bed",row.names=FALSE,col.names=FALSE,quote=FALSE)
