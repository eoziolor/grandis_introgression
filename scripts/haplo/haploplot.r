library(viridis)
library(dplyr)
library(magrittr)
library(Rphylip)
library(ape)
library(stringr)
vcf<-read.table("~/analysis/data/dfst/outliers/zshared_haplo.vcf.bgz",stringsAsFactors = FALSE) #vcf that has been filtered out to only present one allele call per site per individual
sexscore<-read.table("~/analysis/scripts/depth/sexscore",header=TRUE)
cname<-c(seq(1:9),as.character(sexscore[,1])) # colnames for the vcf
colnames(vcf)<-cname

#toss sites that happen to have more than one allele
keep<-!grepl(",",vcf[,5]) #none because of the way I've called haplotypes table to check that

gt<-as.matrix(vcf[,10:297])
class(gt)<-"numeric"

#toss individuals with greater than 90% missing data
keep<-colMeans(is.na(gt))<0.9
gt<-gt[,keep]

#population identifiers
pop<-ifelse(grepl("BB",sexscore[,1]),"BB",
            ifelse(grepl("VB",sexscore[,1]),"VB",
                   ifelse(grepl("PB",sexscore[,1]),"PB",
                          ifelse(grepl("SJ",sexscore[,1]),"SJ",
                                 ifelse(grepl("BNP",sexscore[,1]),"BNP",
                                        ifelse(grepl("SP",sexscore[,1]),"SP",
                                               ifelse(grepl("GB",sexscore[,1]),"GB","WRONG")))))))
popcol<-ifelse(grepl("BB",sexscore[,1]),"black",
            ifelse(grepl("VB",sexscore[,1]),"black",
                   ifelse(grepl("PB",sexscore[,1]),"black",
                          ifelse(grepl("SJ",sexscore[,1]),"firebrick1",
                                 ifelse(grepl("BNP",sexscore[,1]),"firebrick1",
                                        ifelse(grepl("SP",sexscore[,1]),"cadetblue3",
                                               ifelse(grepl("GB",sexscore[,1]),"cadetblue3","WRONG")))))))


pop2<-pop[keep]
popcol2<-popcol[keep]

popname<-c("BB","VB","PB","SJ","BNP","SP","GB")
popnamec<-c("black","black","black","firebrick1","firebrick1","cadetblue3","cadetblue3")
#distance matrix and nj tree
d <- t(gt) %>% dist()
mds <- cmdscale(d)
tr <- nj(d)

#plotting tree
plot(tr,"unrooted",show.tip.label=FALSE)
tiplabels(pch=20,col=popcol2,bg="white",cex=1.2)
legend("topright",pch=20,cex=1.2,legend=popname,col=popnamec)

#plotting mds
plot(mds,pch=20,col=popcol2)

#Plotting both
par(mfrow=c(1,2))
plot(tr,"unrooted",show.tip.label=FALSE)
tiplabels(pch=20,col=popcol2,bg="white",cex=1.2)
plot(mds,pch=20,col=popcol2)
legend("bottomright",pch=20,cex=.8,legend=popname,col=popnamec,y.intersp=.7)

