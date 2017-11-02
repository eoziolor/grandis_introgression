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
bbo<-grep("BB",sexscore$sample)
vbo<-grep("VB",sexscore$sample)
pbo<-grep("PB",sexscore$sample)
sjo<-grep("SJ",sexscore$sample)
bnpo<-grep("BNP",sexscore$sample)
spo<-grep("SP",sexscore$sample)
gbo<-grep("GB",sexscore$sample)


##deletion boundaries:
  #approximate break points:740000-800000
  #will use a conservative window that doesn't include a little spike in coverage: 760000-800000

#no deletion
nodel<-dep[,2]>810000&dep[,2]<890000
del<-dep[,2]>760000&dep[,2]<800000

scalevec <- colSums(dep[nodel,-c(1,2)])/(sum(nodel))
scalevec2<- colSums(dep[del,-c(1,2)])/sum(del)

copies_per_ind <- data.frame(cbind(n_copies=rep(2,288),whichdel=rep(0,288),pop=gsub("-.*","",sexscore[,1]),stringsAsFactors=FALSE))
rownames(copies_per_ind) <- colnames(gts)[-1]

bbcopy<-t(t(dep[del,bbo+2])/scalevec[bbo])%>%
  (function(x){
    s <- colSums(x)
    o <- order(s,decreasing=TRUE)
    o}) %>%
   (gts[,bbo+1])[,.] %>%
  colnames() %>%
  cbind(
    .,
    c(rep(2,28),rep(1,18),rep(0,2)),
    c(rep(0,28),rep(2,7),rep(1,13))
  )

gbcopy<-t(t(dep[del,gbo+2])/scalevec[gbo])%>%
  (function(x){
    s <- colSums(x)
    o <- order(s,decreasing=TRUE)
    o}) %>%
  (gts[,gbo+1])[,.] %>%
  colnames() %>%
  cbind(
    .,
    c(rep(2,28),rep(1,18),rep(0,2))
  )


#my alternative

bbcopy<-t(t(dep[del,bbo+2])/scalevec[bbo])%>%
  (function(x){
    s <- colSums(x)
    o <- order(s,decreasing=TRUE)
    s[o]})

vbcopy<-t(t(dep[del,vbo+2])/scalevec[vbo])%>%
  (function(x){
    s <- colSums(x)
    o <- order(s,decreasing=TRUE)
    s[o]})

pbcopy<-t(t(dep[del,pbo+2])/scalevec[pbo])%>%
  (function(x){
    s <- colSums(x)
    o <- order(s,decreasing=TRUE)
    s[o]})

sjcopy<-t(t(dep[del,sjo+2])/scalevec[sjo])%>%
  (function(x){
    s <- colSums(x)
    o <- order(s,decreasing=TRUE)
    s[o]})

bnpcopy<-t(t(dep[del,bnpo+2])/scalevec[bnpo])%>%
  (function(x){
    s <- colSums(x)
    o <- order(s,decreasing=TRUE)
    s[o]})

spcopy<-t(t(dep[del,spo+2])/scalevec[spo])%>%
  (function(x){
    s <- colSums(x)
    o <- order(s,decreasing=TRUE)
    s[o]})

gbcopy<-t(t(dep[del,gbo+2])/scalevec[gbo])%>%
  (function(x){
    s <- colSums(x)
    o <- order(s,decreasing=TRUE)
    s[o]})

plot(bbcopy,col="black",ylim=c(0,50000),xlim=c(1,50),pch=20,cex=1)
points(vbcopy,col="grey",pch=20,cex=1)
points(pbcopy,col="red",pch=20,cex=1)
points(sjcopy,col="darkorange",pch=20,cex=1)
points(bnpcopy,col="gold",pch=20,cex=1)
points(spcopy,col="cyan",pch=20,cex=1)
points(gbcopy,col="blue",pch=20,cex=1)

abline(h=10000,lty=2,lwd=3,col="red")
abline(h=20000,lty=2,lwd=3,col="grey")
abline(h=40000,lty=2,lwd=3,col="black")
