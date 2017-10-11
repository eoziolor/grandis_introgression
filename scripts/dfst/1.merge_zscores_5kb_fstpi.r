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

#colnam<-names(zfst)[4:24]
# zmerge_temp<-cbind(zfst[1:3],zfst[colnam]-zpi[match(zfst$Scaf,zpi$Scaf),colnam])

zmerge_temp<-zfst[1:24]

for (i in 4:24){
  for (j in 1:1026857){
    zmerge_temp[j,i]<-zfst[j,i]-zpi[j,i]
  }
}

zkeep<-zfst[,25]+zpi[,25]

zmerge<-cbind(zmerge_temp[1:24],zkeep)

write.table(zmerge,"~/analysis/data/dfst/zmerge_fst_pi_5kb",row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")

####### read in the zscore merged table and merge them into windows of interest based on outliers

zmerge<-read.table("~/analysis/data/dfst/zmerge_fst_pi_5kb",sep="\t")

quant<-matrix(nrow=21,ncol=1)

for (i in 1:21)
{
  quant[i,]<-quantile(zmerge[,(i+3)],probs = 0.99,na.rm=TRUE)
} 

############
#now remove all rows with 2 and everything that is below the quantiles

cat ~/analysis/d_fst_analysis/zmerge_fst_pi_5kb | awk '$25>1' | \
grep -v NA | \
awk '$4>4.802706 || $5>4.372692  || $6>5.987987  || $7>5.287479 || $8>4.859941 || $9>5.847102 || $10>5.335611 || $11>6.951480 || $12>4.360933 || $13>4.572003 || $14>4.540147 || $15>5.608441 || $16>4.954220 || $17>5.171635 || $18>4.943881 || $19>3.858078 || $20>3.901198 || $21>3.942247 || $22>4.316083 || $23>4.240024 || $24>4.439515' > ~/analysis/d_fst_analysis/zmerge_outliers_fst_pi_5kb

~/bedtools2/bin/bedtools merge -i ~/analysis/d_fst_analysis/zmerge_outliers_fst_pi_5kb -d 5000 \
-c 4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13,14,14,15,15,16,16,17,17,18,18,19,19,20,20,21,21,22,22,23,23,24,24 \
-o sum,count,sum,count,sum,count,sum,count,sum,count,sum,count,sum,count,sum,count,sum,count,sum,count,sum,count,sum,count,sum,count,sum,count,sum,count,sum,count,sum,count,sum,count,sum,count,sum,count,sum,count \
-g <(cut -f 1-2 ~/analysis/genome/reference_unmasked.fasta.fai) > ~/analysis/d_fst_analysis/zmerge_win_fst_pi_5kb.bed

###continue to reading_merge_outliers.r
