gwas<-read.table(gzfile("~/analysis/data/gwas/gwas.essentials.gz"), header=TRUE,stringsAsFactors = FALSE)
gwas[,3]<-as.numeric(gwas[,3])

palette(c("black","grey","grey30","grey50"))
plot(-log(gwas[,3]),pch=20,cex=.5,col=ifelse(gwas[,3]<0.05,"red",factor(gwas[,1])),
     ylab="-log p for gwas",xlab="chromosomes + unplaced scaffolds",cex.lab=1.5,cex.axis=1.5)
