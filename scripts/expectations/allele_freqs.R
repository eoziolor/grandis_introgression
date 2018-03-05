library(viridis)
library(dplyr)
library(magrittr)
library(Rphylip)
library(ape)
library(stringr)
library(gtools)
library(rtracklayer)
library(ggplot2)
library(cluster)
###Finding overlaps for the bigger regions if we want to----

pbs_out<-read.table("~/analysis/data/dfst/zregions_max.bed",stringsAsFactors = FALSE) #loads a pbs vector with windows merged within 50kb of each other and with max and windows count statistics
names<-c("Scaf","start","end","BBmax","BBcount","VBmax","VBcount","PBmax","PBcount","SJmax","SJcount","BNPmax","BNPcount")
colnames(pbs_out)<-names

all<-pbs_out[,4]>col[1] & pbs_out[,6]>col[2] & pbs_out[,8]>col[3] & pbs_out[,10]>col[4] & pbs_out[,12]>col[5]

write.table(pbs_out[all,1:3],"~/analysis/data/expectations/zall_max.bed",row.names=FALSE,col.names = FALSE,quote = FALSE)
###Take this into kodiak to make a vcf file of the outlier regions----

#Plug that VCF in the following

###Reading in VCF, formatting ----
vcf<-read.table("~/analysis/data/dfst/outliers/zshared_haplo.vcf.bgz",stringsAsFactors = FALSE) #vcf that has been filtered out to only present one allele call per site per individual
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

quantile(with(rr,res-ref),probs=0.9675)
plot(rr[,"res"]-rr[,"ref"],pch=20,cex=.4,col=ifelse(res-ref>.184,"red","black")) #Plotting difference between res and ref
subw<-with(rr,res-ref)>.2

plot(rr[subw,"res"]-rr[subw,"int"],pch=20,cex=.4,col="firebrick2")
points(rr[subw,"res"]-rr[subw,"ref"],pch=20,cex=.4,col="cadetblue3")

#just plotting where intermediates would lie
# plot(with(rr,(res-int)/(res-ref)),pch=20,cex=.5)
# 
# ggplot(data=rr,
#        aes((res-int)/(res-ref)))+
#   geom_histogram()

                 