###Read both keep tables in
zfst<-read.table("~/analysis/data/fst/zfst_keep_5kb",sep='\t')

zpi<-read.table("~/analysis/data/angsd/zpi_keep_5kb",sep="\t")

zfst<-as.data.frame(zfst)
zpi<-as.data.frame(zpi)


zmerge<-matrix(nrow=1026857,ncol=25)

names<- c("Scaf","start","end","bbvb","bbpb","bbsjsp","bbbnp","bbsp","bbgb",
          "vbpb","vbsjsp","vbbnp","vbsp","vbgb",
          "pbsjsp","pbbnp","pbsp","pbgb",
          "sjspbnp","sjspsp","sjspgb",
          "bnpsp","bnpgb",
          "spgb","keep")


colnames(zfst)<-names
colnames(zpi)<-names
colnames(zmerge)<-names


zfst2<-cbind(seq=seq(1:1026857),zfst)
zpi2<-cbind(seq=seq(1:1026857),zpi)

colnam<-names(zfst)[4:24]
zmerge_temp<-cbind(zfst2[1:3],zfst2[colnam]-zpi2[match(zfst2$seq,zpi2$seq),colnam])

# zmerge_temp<-zfst[1:24]
# 
# for (i in 4:5){
#   for (j in 1:1026857){
#     zmerge_temp[j,i]<-zfst[j,i]-zpi[j,i]
#   }h
# }

zkeep<-zfst[,25]+zpi[,25]

zmerge<-cbind(zfst[,1:3],zmerge_temp[4:24],zkeep)

write.table(zmerge,"~/analysis/data/dfst/zmerge_fst_pi_5kb",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
####### read in the zscore merged table and merge them into windows of interest based on outliers

zmerge<-read.table("~/analysis/data/dfst/zmerge_fst_pi_5kb",sep="\t")

quant<-matrix(nrow=21,ncol=1)

for (i in 1:21)
{
  quant[i,]<-quantile(zmerge[,(i+3)],probs = 0.99,na.rm=TRUE)
} 

############
#now remove all rows with 2 and everything that is below the quantiles

cat ~/analysis/data/dfst/zmerge_fst_pi_5kb | awk '$25>1' | \
grep -v NA | \
awk '$4>3.201965 || $5>3.145603  || $6>3.767753  || $7>3.614671 || $8>3.659491 || $9>3.633674 || $10>3.299948 || $11>3.992378 || $12>3.059818 || $13>3.054864 || $14>3.183091 || $15>3.948423 || $16>3.768981 || $17>3.683639 || $18>3.791544 || $19>3.218964 || $20>3.323305 || $21>3.326156 || $22>3.221562 || $23>3.293004 || $24>3.091664' > ~/analysis/data/dfst/zmerge_outliers_fst_pi_5kb

~/program/bedtools2/bin/bedtools merge -i ~/analysis/data/dfst/zmerge_outliers_fst_pi_5kb -d 5000 \
-c 4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13,14,14,15,15,16,16,17,17,18,18,19,19,20,20,21,21,22,22,23,23,24,24 \
-o sum,count,sum,count,sum,count,sum,count,sum,count,sum,count,sum,count,sum,count,sum,count,sum,count,sum,count,sum,count,sum,count,sum,count,sum,count,sum,count,sum,count,sum,count,sum,count,sum,count,sum,count \
-g <(cut -f 1-2 ~/analysis/data/genome/unsplit_merge.fasta.fai) > ~/analysis/data/dfst/zmerge_win_fst_pi_5kb.bed

###continue to reading_merge_outliers.r
