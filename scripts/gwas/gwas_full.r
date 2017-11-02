gwas<-read.table(gzfile("~/analysis/data/gwas/gwas.essentials.chr.gz"), header=TRUE,stringsAsFactors = FALSE)
gwas[,3]<-as.numeric(gwas[,3])

palette(c("black","grey","grey30","grey50"))
plot(-log(gwas[,3]),pch=20,cex=.5,col=ifelse(gwas[,3]<4e-9,"red",factor(gwas[,1])),
     ylab="-log p for gwas",xlab="chromosomes + unplaced scaffolds",cex.lab=1.5,cex.axis=1.5)

sign<-gwas[,3]<4e-9

gsig<-na.omit(gwas[sign,])

write.table(gsig,"~/analysis/data/gwas/significant_regions.txt",row.names = FALSE,col.names = FALSE,quote=FALSE)

####In bash, merging into 5kb windows for easier find

cat significant_regions.txt | awk '{OFS="\t"}{f=chr}{s=$2-1}{print f$1,s,$2,$3}' > sign_gwas.bed

~/program/bedtools2/bin/bedtools merge -i ~/analysis/data/gwas/sign_gwas.bed -d 50000 \
-c 4,4 \
-o mean,count \
-g <(cut -f 1-2 ~/analysis/data/genome/unsplit_merge.fasta.fai) > ~/analysis/data/gwas/sign_merged_50k.bed

####loading windows
gwin<-read.table("~/analysis/data/gwas/sign_merged_50k.bed",header=FALSE)
gsig<-read.table("~/analysis/data/gwas/sign_gwas.bed",header=FALSE)
colnames(gwin)<-c("chr","start","end","mean p","SNPcount")
colnames(gsig)<-c("chr","start","end","p")

###Ranges
source("http://bioconductor.org/biocLite.R")
biocLite()
# biocLite("GenomicFeatures")
library(GenomicRanges)
biocLite("GenomicFeatures")

species<-"F.heteroclitus"
fhet_txdb<-makeTranscriptDbFromGFF("~/analysis/data/genome/unsplit_ncbi.gff",species)

gr<-GRanges(seqname=gwin[,1], ranges=IRanges(start=gwin[,2],end=gwin[,3]))
gs<-GRanges(seqname=gsig[,1], ranges=IRanges(start=gsig[,3],end=gsig[,3]))
