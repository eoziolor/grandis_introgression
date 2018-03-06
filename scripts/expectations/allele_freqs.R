library(viridis)
library(dplyr)
library(magrittr)
library(Rphylip)
library(ape)
library(stringr)
library(gtools)
library(rtracklayer)
library(ggplot2)

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

###Ranking individuals based on those SNPs----

vcfo<-vcf[subw,]
zn<-na.omit(z[all,1:3])
zn2<-matrix(nrow = length(zn[,1]),ncol=288)
zn<-cbind(zn,zn2)
nam<-c("chr","start","end",colnames(vcf[,10:297]))
colnames(zn)<-nam

vcfo<-cbind(vcfo[,1],vcfo[,2]-1,vcfo[,2:297])
write.table(vcfo[,1:3],"~/analysis/data/expectations/vcf_out.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(zn[,1:3],"~/analysis/data/expectations/zshared_all.bed",row.names = FALSE,col.names = FALSE,quote = FALSE)

bed1=import("~/analysis/data/expectations/vcf_out.bed")
bedall=import("~/analysis/data/expectations/zshared_all.bed")
bed1overlall=bed1[bed1 %over% bedall]
hitsall<-countOverlaps(bedall,bed1)
subset<-subsetByOverlaps(bedall,bed1)

gt<-as.matrix(vcfo[,10:297]) #convert to numeric
class(gt)<-"numeric"
keep<-colMeans(is.na(gt))<0.9 #toss sites with low coverage
gt<-gt[,keep]

fr<-c()
j=0
for(i in 1:length(hitsall)){
  if(hitsall[i]>1){
  f<-as.numeric(colSums(gt[(j+1):(hitsall[i]+j),],na.rm=TRUE)/hitsall[i])
  fr<-rbind(fr,f)
  j=j+hitsall[i]
  print(j==sum(hitsall[1:i]))
  print(i)
}
}

# for(i in 1:length(zn[,1])){
#   j=1
#   while(zn[i,1]==vcfo[j,1]){
#     if(zn[i,2]<=vcfo[j,2] && zn[i,3]>=vcfo[j,2])
#     zn[i,4:291]<-zn[i,4:291]+vcfo[j,10:297]
#     j=j+1
#   }
# }

as.numeric(colSums(gt[sum(hitsall[1:3]):sum(hitsall[1:4]),],na.rm=TRUE)/hitsall[4])