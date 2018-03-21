install.packages("dplyr")
library(viridis)
library(dplyr)
library(magrittr)
library(stringr)
library(gtools)
library(rtracklayer)
library(ggplot2)
library(ggsci)

###Finding overlaps for the bigger regions if we want to----
pbs_dist<-read.table("~/analysis/data/dfst/zpbs_5kb.bed",header=FALSE)
distnames<-c("Scaf","start","end","bb","vb","pb","sj","bnp")
colnames(pbs_dist)<-distnames
col<-c()
for (i in 1:5){
  col[i]<-quantile(pbs_dist[,i+3],prob=.99,na.rm=TRUE)
}

###Overall outliers
z<-read.table("~/analysis/data/expectations/zregions_max.bed",stringsAsFactors = FALSE) #loads a pbs vector with windows merged within 50kb of each other and with max and windows count statistics
names<-c("Scaf","start","end","BBmax","BBcount","VBmax","VBcount","PBmax","PBcount","SJmax","SJcount","BNPmax","BNPcount")
colnames(z)<-names

all<-z[,4]>col[1] & z[,6]>col[2] & z[,8]>col[3] & z[,10]>col[4] & z[,12]>col[5]

write.table(na.omit(z[all,1:3]),"~/analysis/data/expectations/zall_max.bed",row.names=FALSE,col.names = FALSE,quote = FALSE)
###Take this into kodiak to make a vcf file of the outlier regions----

#Plug that VCF in the following

###Reading in VCF, formatting ----
vcf<-read.table("~/analysis/data/expectations/zall_max_haplo.vcf.bgz",stringsAsFactors = FALSE) #vcf that has been filtered out to only present one allele call per site per individual
sexscore<-read.table("~/analysis/scripts/depth/sexscore",header=TRUE)
#vcf<-cbind(vcf[,1],vcf[,2]-1,vcf[,2:297]) #Converting the first three columns into bed type format
cname<-c(seq(1:9),as.character(sexscore[,1])) # colnames for the vcf
colnames(vcf)<-cname

subw<-nchar(vcf[,5])<2 #using only biallelic SNPs

vcf<-vcf[subw,] #restricting to biallelic SNPs

table(nchar(vcf[,4])>2) #Confirming biallelic SNPs with reference allele after restriction

gt<-as.matrix(vcf[,10:297]) #convert to numeric
class(gt)<-"numeric"

keep<-colMeans(is.na(gt))<0.9 #toss sites with low coverage
gt<-gt[,keep]

#Summing SNP frequencies----

popfr<-cbind(vcf[,1:2])
nam<-c("chr","pos","BB","VB","PB","SJ","BNP","SP","GB")

for(j in 1:length(nam[-1:-2])){
  b<-gt[,grepl(nam[j+2],colnames(gt))]
  f<-as.numeric(rowSums(gt[,grepl(nam[j+2],colnames(gt))],na.rm=TRUE)/ncol(b))
  popfr<-cbind(popfr,f)
  print(nam[j+2])
  print(ncol(b))
}
#summed frequencies over all individuals 
colnames(popfr)<-nam

res<-(rowSums(popfr[,3:5])/3) #summing SNP frequencies over resistant sites and below for all other
int<-(rowSums(popfr[,6:7])/2)
ref<-(rowSums(popfr[,8:9])/2)
rr<-cbind(popfr[,1:2],res,ref,int) #merging these into a single DF

quantile(with(rr,res-ref),probs=0.97)
plot(rr[,"res"]-rr[,"ref"],pch=20,cex=.4,col=ifelse(res-ref>.13,"red","black")) #Plotting difference between res and ref
subw<-with(rr,res-ref)>.13

plot(rr[subw,"res"]-rr[subw,"int"],pch=20,cex=.4,col="firebrick2")
points(rr[subw,"res"]-rr[subw,"ref"],pch=20,cex=.4,col="cadetblue3")

#just plotting where intermediates would lie
# plot(with(rr,(res-int)/(res-ref)),pch=20,cex=.5)
# 
# ggplot(data=rr,
#        aes((res-int)/(res-ref)))+
#   geom_histogram()

###Frequency of haplotypes in individuals based on those SNPs----

vcfo<-vcf[subw,] #using only the biallelic SNPs that are in the 97th percentile of allele frequency difference
zn<-na.omit(z[all,1:3]) #using the rows of merged divergent regions to create regions from which I will merge SNPs
nam<-c("chr","start","end")
colnames(zn)<-nam

vcfo<-cbind(vcfo[,1],vcfo[,2]-1,vcfo[,2:297]) #making vcfo into a bed file
write.table(vcfo[,1:3],"~/analysis/data/expectations/vcf_out.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(zn[,1:3],"~/analysis/data/expectations/zshared_all.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)

bed1=import("~/analysis/data/expectations/vcf_out.bed")
bedall=import("~/analysis/data/expectations/zshared_all.bed")
bed1overlall=bed1[bed1 %over% bedall]
hitsall<-countOverlaps(bedall,bed1) #crScienceeating a vector of overlap counts - # SNPs present in that overlap

gt<-as.matrix(vcfo[,10:297]) #convert to numeric
class(gt)<-"numeric"
keep<-colMeans(is.na(gt))<0.9 #toss sites with low coverage
gt<-gt[,keep]

fr<-c() #this is used to create a matrix of frequencies of the haplotype in each outlier window by counting the number of SNPs in it
j=0
for(i in 1:length(hitsall)){
  if(hitsall[i]>1){
  f<-as.numeric(colSums(gt[(j+1):(hitsall[i]+j),],na.rm=TRUE)/hitsall[i])
  fr<-rbind(fr,f)
  j=j+hitsall[i]
  print(j==sum(hitsall[1:i]))
  print(i)
  } else if (hitsall[i]==1){
    j=j+1
  }
}

colnames(fr)<-colnames(gt)

#Calling individuals as 0,0.5,1 frequency of haplotype ----
#I'll make a judgment call that anything below 0.1 is absence, anything above 0.6 is presence in full, anything in between is half
frr<-fr
for(i in 1:dim(fr)[[2]]){
  for(j in 1:dim(fr)[[1]]){
    if(fr[j,i]<0.1){
      frr[j,i]<-as.numeric(0)
    } else if (fr[j,i]>0.6){
      frr[j,i]<-as.numeric(1)
    } else{
      frr[j,i]<-as.numeric(0.5)
    }
  }
}

###Summing allele frequencies over populations----

nam<-c("BB","VB","PB","SJ","BNP","SP","GB")
subw<-hitsall>1

pfr<-zn[subw,1:3]
for(i in 1:length(nam)){
  f<-as.numeric(rowSums(frr[,grepl(nam[i],colnames(frr))],na.rm=TRUE))/ncol(frr[,grepl(nam[i],colnames(frr))])
  pfr<-cbind(pfr,f)
}

colnames(pfr)<-c(colnames(pfr[1:3]),nam)

pfr<-t(pfr[4:10])

# plot(pfr[,1],ylim=c(0,.5))
# for(i in 1:dim(pfr)[[2]]){
#   lines(pfr[,i],ylim=c(0,.5))
# }

#Seeing if those frequencies meet expectations from admixture----
a<-read.table("~/analysis/data/fastngs/alpha",header=TRUE)
a<-as.numeric(a)

e<-c() #creating a vector e of expected values for each allele
r<-c()
for(i in 1:dim(pfr)[[2]]){
  for(j in 1:dim(pfr)[[1]]){
    b<-a[j]*pfr[1,i]+(1-a[j])*pfr[7,i]
    r[j]<-b
  }
  e<-cbind(e,r)
}

rownames(e)<-rownames(pfr)
colnames(e)<-colnames(pfr)

efr<-pfr-e
colnames(efr)<-(zn[subw,1])
efr1<-efr[,grepl("chr",colnames(efr))]

par(mar=c(4,4,2,2))
plot(efr[,1],ylim=c(-.2,.4),type='n',xaxt='n',
     ylab="Deviation from expected allele frequencies due to admixture",bty='l')
for(i in 1:dim(efr)[[2]]){
  lines(efr[,i],ylim=c(0,.5),lwd=ifelse(zn[i,1]=="chr1",2.5,0.8)
        ,col=ifelse(zn[i,1]=="chr1","red","black"))
}
axis(side=1,at=1:7,labels=c(rownames(efr)))
abline(h=1,lwd=3,lty=2,col="black")

#Color by chromosome ----
chr<-unique(colnames(efr1))
# color = pal_igv(palette=c("default"),alpha=1)(17)
color=pal_ucscgb(palette=c("default"),alpha=1)(17)
col=color[1:17]
chr<-cbind(chr,col)
rownames(chr)<-chr[,1]

par(mar=c(4,4,2,2))
plot(efr1[,1],ylim=c(-.2,.4),type='n',xaxt='n',
     ylab="Deviation from expected allele frequencies due to admixture",bty='l')
for(i in 1:dim(efr1)[[2]]){
  lines(efr1[,i],ylim=c(0,.4),lwd=1.5,col=chr[colnames(efr1)[i],2])
}
axis(side=1,at=1:7,labels=c(rownames(efr1)))
abline(h=1,lwd=3,lty=2)
legend("topright",legend=c(chr[,1]),col=c(chr[,2]),pch=20,cex=1)
