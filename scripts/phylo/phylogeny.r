#this function pulls genotypes from a tabix-indexed vcf containing biallelic snps
#it tosses most of the information and produces a table of positions, and genotypes coded as 0,1,2
#install.packages("dplyr")
library(dplyr)
#install.packages("mscr")
library(mscr)
#install.packages("Rphylip")
library(Rphylip)

sexscore<-read.table("~/analysis/scripts/depth/sexscore",header=TRUE)
sexscore<-as.data.frame(sexscore)
rownames(sexscore)<-sexscore[,1]


mat <- read.table("~/analysis/scripts/phylo/popcolors.txt",stringsAsFactors=FALSE,row.names=1)


pullgenos_ALL<-function(){
  
  header<-scan("~/analysis/data/phylo/header.txt",what="character")
  gts<-read.table(pipe("gunzip -c ~/analysis/data/phylo/25Mb.vcf.bgz | sed -n '1~2p' "),stringsAsFactors=FALSE) #deleted g from gsed
  
  colnames(gts)<-header
  gts2<-as.matrix(gts[,10:297])
  gts2<-gsub(":.*","",gts2)
  gts2[gts2=="."]<-NA
  gts2[gts2=="0/0"]<-0
  gts2[gts2=="0/1"]<-1
  gts2[gts2=="1/1"]<-2
  
  gts2<-cbind(gts[,2],gts2)
  class(gts2)<-"numeric"
  return(gts2)
}

freq <- function(df,pop,sex="M|F"){
  sexscore <- sexscore[sexscore[,1]%in%colnames(df),]
  vec <- c(FALSE,grepl(sex,sexscore$sex))&grepl(pop,colnames(df))
  apply(df[,vec],
        MAR=1,
        FUN=function(x){
          sum(x,na.rm=TRUE)/sum(!is.na(x))/2
        }
  )
}

pco <- function(df,pop,sex="M|F"){
  sexscore <- sexscore[sexscore[,1]%in%colnames(df),]
  vec <- c(FALSE,grepl(sex,sexscore$sex))&grepl(pop,colnames(df))
  dfn <- apply(df[,vec],
               MAR=1,
               FUN=function(x){
                 sum(x,na.rm=TRUE)
               }
  )
  
  dfd <- apply(df[,vec],
               MAR=1,
               FUN=function(x){
                 sum(!is.na(x))*2
               }
  )
  out <- list(co=dfn,sa=dfd)
  return(out)
}


gts <- pullgenos_ALL()
pops <- colnames(gts[,-1]) %>% gsub("-.*","",.)

popfreqs <- cbind(
  freq(gts,"BB",sex="M|F"),
  freq(gts,"VB",sex="M|F"),
  freq(gts,"PB",sex="M|F"),
  freq(gts,"SJ",sex="M|F"),
  freq(gts,"BNP",sex="M|F"),
  freq(gts,"SP",sex="M|F"),
  freq(gts,"GB",sex="M|F")
)
colnames(popfreqs) <- c("BB","VB","PB","SJ","BNP","SP","GB")
popfreqs <- popfreqs[which(rowSums(is.na(popfreqs))==0),]

#write.table(gts,"~/analysis/data/phylo/gts",col.names=TRUE)
#write.table(popfreqs,"~/share/phylogeny/popfreqs",col.names =TRUE)

#gts2<-read.table("~/analysis/data/phylo/gts",header=TRUE)
#popfreqs<-read.table("~/share/phylogeny/popfreqs",header=TRUE)

#using contml
subl <- gts[,1]>0
tr <- Rcontml(X=t(popfreqs),path="~/phylip-3.696/exe")

subi <- which(colMeans(is.na(gts[,-1]))<0.8)
names(subi) <- NULL

d <- dist(t(gts[,subi+1]))

# par(mfrow=c(1,3))
# fit <- cmdscale(d,eig=TRUE, k=5)
# plot(fit$points[,1],fit$points[,2],pch=20,col=mat[subi,1],cex=1.5,
#      xlab="First Dimension",ylab="Second Dimension",bty="l")
# popcol <- c("black","grey","red","darkorange","gold","cyan","blue")
# popnames <- c("BB","VB","PB","SJ","BNP","SP","GB")
# legend("topright",legend=popnames,col=popcol,pch=20,bty="n",cex=1.7,
#        x.intersp=.5,y.intersp = .5)

par(mfrow=c(1,1),mar=c(0,0,0,0))
plot(tr,type="unrooted",show.tip.label=FALSE)
popord <- c(3, 4, 5, 6, 7, 2, 1)
popcol<-c("black","black","black","firebrick2","firebrick2","cadetblue3","cadetblue3")
tiplabels(pch=20,col=popcol[popord],cex=4)
legend("topright",legend=popnames,col=popcol,pch=20,bty="n",cex=1.7,
       x.intersp=.5,y.intersp = .5)

# tr2 <- nj(d)
# plot(tr2,show.tip.label=FALSE,type="unrooted")
# #tr2$tip.label %>% gsub("-.*","",.) %>% mat[.,1] %>%
# tiplabels(pch=20,col=mat[subi,1],cex=1.5)