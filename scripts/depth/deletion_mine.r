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

gts<-pullgenos("chr1:650000-900000")
hap<-pullgenos.phased("chr1:650000-900000")
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

copies_per_ind <- data.frame(cbind(n_copies=rep(2,288),pop=gsub("-.*","",sexscore[,1])),stringsAsFactors=FALSE)
rownames(copies_per_ind) <- colnames(gts)[-1]

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
plot(vbcopy,col="grey",pch=20,cex=1)
plot(pbcopy,col="red",pch=20,cex=1)
plot(sjcopy,col="darkorange",pch=20,cex=1)
plot(bnpcopy,col="gold",pch=20,cex=1)
plot(spcopy,col="cyan",pch=20,cex=1)
plot(gbcopy,col="blue",pch=20,cex=1)

msp<-mean(spcopy[1:47])
ssp<-sd(spcopy[1:47])

abline(h=msp-2*ssp,lty=2,lwd=3,col="black")
abline(h=10000,lty=2,lwd=3,col="grey")


###Counting the copies

bbnum<-t(t(dep[del,bbo+2])/scalevec[bbo])%>%
  (function(x){
    s <- colSums(x)
    o <- order(s,decreasing=TRUE)
    o}) %>%
  (gts[,bbo+1])[,.] %>%
  colnames() %>%
  cbind(
    .,
    c(rep(2,18),rep(1,5),rep(0,1))
  )

copies_per_ind[bbnum[,1],1]<-bbnum[,2]

vbnum<-t(t(dep[del,vbo+2])/scalevec[vbo])%>%
  (function(x){
    s <- colSums(x)
    o <- order(s,decreasing=TRUE)
    o}) %>%
  (gts[,vbo+1])[,.] %>%
  colnames() %>%
  cbind(
    .,
    c(rep(2,38),rep(1,10),rep(0,1))
  )

copies_per_ind[vbnum[,1],1]<-vbnum[,2]

pbnum<-t(t(dep[del,pbo+2])/scalevec[pbo])%>%
  (function(x){
    s <- colSums(x)
    o <- order(s,decreasing=TRUE)
    o}) %>%
  (gts[,pbo+1])[,.] %>%
  colnames() %>%
  cbind(
    .,
    c(rep(2,21),rep(1,22),rep(0,4))
  )

copies_per_ind[pbnum[,1],1]<-pbnum[,2]

sjnum<-t(t(dep[del,sjo+2])/scalevec[sjo])%>%
  (function(x){
    s <- colSums(x)
    o <- order(s,decreasing=TRUE)
    o}) %>%
  (gts[,sjo+1])[,.] %>%
  colnames() %>%
  cbind(
    .,
    c(rep(2,1),rep(1,4),rep(0,19))
  )

copies_per_ind[sjnum[,1],1]<-sjnum[,2]

bnpnum<-t(t(dep[del,bnpo+2])/scalevec[bnpo])%>%
  (function(x){
    s <- colSums(x)
    o <- order(s,decreasing=TRUE)
    o}) %>%
  (gts[,bnpo+1])[,.] %>%
  colnames() %>%
  cbind(
    .,
    c(rep(2,0),rep(1,5),rep(0,43))
  )

copies_per_ind[bnpnum[,1],1]<-bnpnum[,2]

spnum<-t(t(dep[del,spo+2])/scalevec[spo])%>%
  (function(x){
    s <- colSums(x)
    o <- order(s,decreasing=TRUE)
    o}) %>%
  (gts[,spo+1])[,.] %>%
  colnames() %>%
  cbind(
    .,
    c(rep(2,0),rep(1,0),rep(0,48))
  )

copies_per_ind[spnum[,1],1]<-spnum[,2]

gbnum<-t(t(dep[del,gbo+2])/scalevec[gbo])%>%
  (function(x){
    s <- colSums(x)
    o <- order(s,decreasing=TRUE)
    o}) %>%
  (gts[,gbo+1])[,.] %>%
  colnames() %>%
  cbind(
    .,
    c(rep(2,0),rep(1,1),rep(0,47))
  )

copies_per_ind[gbnum[,1],1]<-gbnum[,2]

copies_per_ind[,1] <- as.numeric(copies_per_ind[,1])


# table(copies_per_ind[,c(1,2)]) %>% 
#   (function(x){as.array(x[,c(3,6,2,5,4,7,1)])}) %>% 
#   (function(x){barplot(x,beside=TRUE,col=c("darkolivegreen1","chartreuse3","darkgreen"),
#                        ylab='',cex.axis = 2.5,cex.names = 1.5)})
par(mfrow=c(1,1),mar=c(2,4,2,4))
copies<-as.array(table(copies_per_ind[,c(1,2)]) %>% 
                (function(x){x[,c(3,6,2,5,4,7,1)]}) )
copies.prop<-prop.table(copies,2)*100 
barplot(copies.prop,beside=TRUE,col=c("darkolivegreen1","chartreuse3","darkgreen"),
       ylab='',cex.axis = 2.2,cex.names = 1.5)

# legend("topright",legend=c("wt","delHet","delHom"),pch=20,cex=1.2,col=c("blue","red","pink","green","lightgrey","grey50","black"),
#        x.intersp = 0.4,y.intersp = .7)

###smoothing vector function

subsmooth <- function(vec,by=10,width=1000){
  
  len <- length(vec)
  subl <- seq(from=by,to=len,by=by)
  submax <- length(subl)
  width <- width/2
  test <- vec[subl]
  
  for(i in 1:submax){
    
    j <- i - width
    k <- i + width
    if(j < 1) {j <- 1}
    if(k > submax) {k <- submax}
    test[i] <- mean(test[j:k],na.rm=TRUE)
  }
  
  return(test)
  
}


###smoothing pop coverage

depsub <- as.data.frame(apply(dep[,-c(1,2)],MAR=2,FUN=subsmooth,by=20,width=500))
jump<-dep[seq(from=20,to=dim(dep)[1],by=20),2]
depsub <- cbind("Chr1sub",jump,depsub)

par(mfrow=c(7,1),mar=c(2,4,1.5,0.5))


bb <- grep("BB",sexscore$sample)
t(t(depsub[,bbo+2])/scalevec[bbo]) %>%
  (function(x){
    s <- colSums(x)
    o <- order(s,decreasing=TRUE)
    x[,o]
  }) %>%
  
  image(z=.,zlim=c(0,4),col=viridis(30),x=depsub[,2],y=(1:length(bbo)),ylab="BB")

vb <- grep("vb",sexscore$sample)
t(t(depsub[,vbo+2])/scalevec[vbo]) %>%
  (function(x){
    s <- colSums(x)
    o <- order(s,decreasing=TRUE)
    x[,o]
  }) %>%
  
  image(z=.,zlim=c(0,4),col=viridis(30),x=depsub[,2],y=(1:length(vbo)),ylab="VB")

pb <- grep("pb",sexscore$sample)
t(t(depsub[,pbo+2])/scalevec[pbo]) %>%
  (function(x){
    s <- colSums(x)
    o <- order(s,decreasing=TRUE)
    x[,o]
  }) %>%
  
  image(z=.,zlim=c(0,4),col=viridis(30),x=depsub[,2],y=(1:length(pbo)),ylab="PB")

sj <- grep("sj",sexscore$sample)
t(t(depsub[,sjo+2])/scalevec[sjo]) %>%
  (function(x){
    s <- colSums(x)
    o <- order(s,decreasing=TRUE)
    x[,o]
  }) %>%
  
  image(z=.,zlim=c(0,4),col=viridis(30),x=depsub[,2],y=(1:length(sjo)),ylab="SJ")

bnp <- grep("bnp",sexscore$sample)
t(t(depsub[,bnpo+2])/scalevec[bnpo]) %>%
  (function(x){
    s <- colSums(x)
    o <- order(s,decreasing=TRUE)
    x[,o]
  }) %>%
  
  image(z=.,zlim=c(0,4),col=viridis(30),x=depsub[,2],y=(1:length(bnpo)),ylab="BNP")

sp <- grep("sp",sexscore$sample)
t(t(depsub[,spo+2])/scalevec[spo]) %>%
  (function(x){
    s <- colSums(x)
    o <- order(s,decreasing=TRUE)
    x[,o]
  }) %>%
  
  image(z=.,zlim=c(0,4),col=viridis(30),x=depsub[,2],y=(1:length(spo)),ylab="SP")

gb <- grep("gb",sexscore$sample)
t(t(depsub[,gbo+2])/scalevec[gbo]) %>%
  (function(x){
    s <- colSums(x)
    o <- order(s,decreasing=TRUE)
    x[,o]
  }) %>%
  
  image(z=.,zlim=c(0,4),col=viridis(30),x=depsub[,2],y=(1:length(gbo)),ylab="GB")

par(mfrow=c(2,1),mar=c(2,4,1.5,0.5))

###Sex

f <- grep("F",sexscore$sex)
t(t(depsub[,f+2])/scalevec[f]) %>%
  (function(x){
    s <- colSums(x)
    o <- order(s,decreasing=TRUE)
    x[,o]
  }) %>%
  
  image(z=.,zlim=c(0,4),col=viridis(30),x=depsub[,2],y=(1:length(f)),ylab="F")


m <- grep("M",sexscore$sex)
t(t(depsub[,m+2])/scalevec[m]) %>%
  (function(x){
    s <- colSums(x)
    o <- order(s,decreasing=TRUE)
    x[,o]
  }) %>%
  
  image(z=.,zlim=c(0,4),col=viridis(30),x=depsub[,2],y=(1:length(m)),ylab="M")


####representative graphs
#make representative coverage figures for three individuals,yes/het/no duplication
#264 - homozygous BB
#277 - heterozygous SJ
#19 - no deletion SP
par(mfrow=c(3,1),mar=c(2,2,1.5,0.5))

ind <- 264
subl <- seq(from=10,to=177830,by=10)
plot(dep[subl,2],dep[subl,ind]/scalevec[ind-2],pch=20,cex=.2,ylim=c(0,3))
abline(h=seq(from=0.5,to=2.5,by=0.5),col="darkgray",lwd=2)
points(dep[subl,2],subsmooth(dep[,ind])/scalevec[ind-2],pch=20,cex=.2,col="red")

ind <- 277
subl <- seq(from=10,to=177830,by=10)
plot(dep[subl,2],dep[subl,ind]/scalevec[ind-2],pch=20,cex=.2,ylim=c(0,3))
abline(h=seq(from=0.5,to=2.5,by=0.5),col="darkgray",lwd=2)
points(dep[subl,2],subsmooth(dep[,ind])/scalevec[ind-2],pch=20,cex=.2,col="red")

ind <- 19
subl <- seq(from=10,to=177830,by=10)
plot(dep[subl,2],dep[subl,ind]/scalevec[ind-2],pch=20,cex=.2,ylim=c(0,3))
abline(h=seq(from=0.5,to=2.5,by=0.5),col="darkgray",lwd=2)
points(dep[subl,2],subsmooth(dep[,ind])/scalevec[ind-2],pch=20,cex=.2,col="red")


##Alternative plot
ind1 <- 264
ind2 <- 277
ind3 <- 19
par(mfrow=c(1,1),mar=c(4,4,1.5,0.5))
subl <- seq(from=10,to=177830,by=10)
plot(dep[subl,2],subsmooth(dep[,ind1])/scalevec[ind1-2],pch=20,cex=.3,col="black",ylim=c(0,2.2),
     ylab="Smoothed region coverage",xlab="Location on chromosome1")
abline(h=seq(from=0.5,to=2.5,by=0.5),col="darkgray",lwd=2)
points(dep[subl,2],subsmooth(dep[,ind2])/scalevec[ind2-2],pch=20,cex=.3,col="firebrick2")
points(dep[subl,2],subsmooth(dep[,ind3])/scalevec[ind3-2],pch=20,cex=.3,col="cadetblue3")

legend("bottomright",lty=2,lwd=2,legend=c("wt/wt S2","del/wt IH1","del/del R1"),col=c("cadetblue3","firebrick2","black"),cex=1.4,
       x.intersp = .7,y.intersp = 1)
abline(v=c(729500,806000),lty=2,lwd=2,col="purple")
abline(v=c(728500,730500,765000,767000,805000,807000),lty=2,lwd=2,col="purple")


###returns population frequency given pop,sex,data

freq <- function(df,pop){
  vec <-grepl(pop,colnames(df))
  apply(df[,vec],
        MAR=1,
        FUN=function(x){
          sum(x,na.rm=TRUE)/sum(!is.na(x))/2
        }
  )
}

###make a tree!

Fpopfreqs <- cbind(
  freq(gts,"BB"),
  freq(gts,"VB"),
  freq(gts,"PB"),
  freq(gts,"SJ"),
  freq(gts,"BNP"),
  freq(gts,"SP"),
  freq(gts,"GB")
)
colnames(Fpopfreqs) <- c("BB","VB","PB","SJ","BNP","SP","GB")
Fpopfreqs <- Fpopfreqs[which(rowSums(is.na(Fpopfreqs))==0),]

#using contml
tr <- Rcontml(X=t(Fpopfreqs),path="~/phylip-3.696/exe")

#already ran bootstraps, don't re-run unless necessary. 
#trboot <- boot.phylo(tr,t(Fpopfreqs),FUN=function(x){Rcontml(X=x,path="~/phylip-3.696/exe")},trees=TRUE)

#for some reason boot.phylo calculates BP wrong. 
plot(tr,type="unrooted",show.tip.label=FALSE)
# nodelabels(trboot$BP)
popord <- c(3, 4, 5, 6, 7, 2, 1)
popcol <- c("black","grey","red","darkorange","gold","cyan","blue")
popnames <- c("BB","VB","PB","SJ","BNP","SP","GB")
tiplabels(pch=20,col=popcol[popord],cex=6)
legend("topright",legend=popnames,col=popcol,pch=20,cex=1)

###


trN <- Rgendist(t(Fpopfreqs),path="~/phylip-3.696/exe/") %>% 
  Rneighbor(.,path="~/phylip-3.696/exe/")


trbootN <- boot.phylo(trN,t(Fpopfreqs),FUN=function(x){
  Rgendist(t(Fpopfreqs),path="~/phylip-3.696/exe/") %>% 
    Rneighbor(.,path="~/phylip-3.696/exe/")
},trees=TRUE)

plot(trN,type="unrooted")
nodelabels(trbootN$BP)
