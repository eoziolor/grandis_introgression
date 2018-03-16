#R script to start new scaffold analysis in R
read.table("~/analysis/data/plink/full.genome/full.genome.mds",stringsAsFactors=FALSE)->full
read.table("~/analysis/data/plink/metadata",stringsAsFactors=FALSE,sep="\t")->meta
row.names(meta)<-meta[,1]
read.table("~/analysis/data/plink/popcolors.txt",stringsAsFactors=FALSE,sep="\t")->pops
row.names(pops)<-pops[,1]
plot(full[,4:5],pch=20,cex=.8,col=pops[meta[full[,1],2],2])
library(scatterplot3d)
scatterplot3d(full[,4:6],pch=20,color=pops[meta[full[,1],2],2])

#Separating by sex
read.table("~/analysis/data/plink/popcolors_sex.txt",stringsAsFactors=FALSE,sep="\t")->pops_sex
row.names(pops_sex)<-pops_sex[,1]
plot(full[,4:5],pch=20,cex=.8,col=pops_sex[meta[full[,1],3],2])
scatterplot3d(full[,4:6],pch=20,color=pops_sex[meta[full[,1],3],2])


#Separating by collection year
read.table("~/share/metadata_2.txt",stringsAsFactors=FALSE,sep="\t")->meta
row.names(meta)<-meta[,1]
read.table("~/share/popcolors_2.txt",stringsAsFactors=FALSE,sep="\t")->pops
row.names(pops)<-pops[,1]
plot(full[,4:5],pch=20,cex=.8,col=pops[meta[full[,1],3],2])
library(scatterplot3d)
scatterplot3d(full[,4:6],pch=20,color=pops[meta[full[,1],3],2])
