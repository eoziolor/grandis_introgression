# install.packages('viridis')
# install.packages('dplyr')
# install.packages('magrittr')
# install.packages('Rphylip')
# install.packages('ape')
# install.packages('stringr')

library(viridis)
library(dplyr)
library(magrittr)
library(Rphylip)
library(ape)
library(stringr)

source("~/analysis/scripts/depth/pullgenos.r")
source("~/analysis/scripts/depth/pullgenos.phased.r")

gts<-pullgenos("chr1:700000-900000")
hap<-pullgenos.phased("chr1:700000-900000")
colnames(gts)[1] <- "pos"
dep<-read.table("~/analysis/data/depth/coverage_ahr.txt.gz",stringsAsFactors=FALSE)
sexscore<-read.table("~/analysis/scripts/depth/sexscore",header=TRUE)
mat<-read.table("~/analysis/data/plink/popcolors.txt",stringsAsFactors=FALSE,sep="\t")
pops <- gsub("-.*","",sexscore[,1])
popord <- c(grep("BB",sexscore$sample),grep("VB",sexscore$sample),grep("PB",sexscore$sample),grep("SJ",sexscore$sample),
            grep("BNP",sexscore$sample),grep("SP",sexscore$sample),grep("GB",sexscore$sample))

##deletion boundaries:
  #approximate break points:740000-800000
  #will use a conservative window that doesn't include a little spike in coverage: 760000-800000

#no deletion
nodel<-dep[,2]>810000&dep[,2]<890000
del<-dep[,2]>760000&dep[,2]<800000

scalevec <- colSums(dep[nodel,-c(1,2)])/(sum(nodel))

copies_per_ind <- data.frame(cbind(n_copies=rep(2,288),whichdel=rep(0,288),pop=gsub("-.*","",colnames(gts)[-1])),stringsAsFactors=FALSE)
rownames(copies_per_ind) <- colnames(gts)[-1]

