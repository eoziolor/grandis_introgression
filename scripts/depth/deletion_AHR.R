#scaffold10074 analysis
library(mscr)
library(viridis)
library(dplyr)
library(magrittr)
library(Rphylip)
library(ape)
library(stringr)


data(sexscore)


#gts <- pullgenos("Scaffold10074:500000-900000")
#hap <- pullgenos.phased("Scaffold10074:500000-900000")
#colnames(gts)[1] <- "pos"
#dep <- read.table(pipe("samtools depth -f ~/popgen/alignments/bowtie/ALL/meta/all.list -Q 40 -r Scaffold10074:500000-900000"),stringsAsFactors=FALSE)
#colnames(dep) <- c("scaffold",colnames(gts))

#rna <- read.table(pipe("samtools depth -f ~/tolerance_rnaseq/bowtie.bams.list -Q 40 -r Scaffold10074:500000-900000"),stringsAsFactors=FALSE)


# #remove sites with very high coverage
# subbi <- rowSums(dep[,grep("BI",colnames(dep))])<400
# plot(rowSums(dep[subbi,grep("BI",colnames(dep))]),pch=20,cex=.2)

# subf <- rowSums(dep[,grep("F",colnames(dep))])<55
# plot(rowSums(dep[subf,grep("F",colnames(dep))]),pch=20,cex=.2)

# subsh <- rowSums(dep[,grep("SH",colnames(dep))])<50
# plot(rowSums(dep[subsh,grep("SH",colnames(dep))]),pch=20,cex=.2)

# subkc <- rowSums(dep[,grep("KC",colnames(dep))])<50
# plot(rowSums(dep[subkc,grep("KC",colnames(dep))]),pch=20,cex=.2)

# dep <- dep[subbi&subf&subsh&subkc,]

load("Scaffold10074.Rdata")

mat <- read.table("symbols_colors.txt",row.names=1,stringsAsFactors=FALSE)
pops <- gsub("-.*","",sexscore[,1])
popord <- c(grep("NBH",sexscore$sample),grep("BI",sexscore$sample),grep("BP",sexscore$sample),grep("F",sexscore$sample),grep("NYC",sexscore$sample),grep("SH",sexscore$sample),grep("ER",sexscore$sample),grep("KC",sexscore$sample))
# #remove sites with very low coverage from rna-seq
# rnanames <- scan("~/tolerance_rnaseq/bowtie.bams.list",what="character",sep="\n")
# rnanames <- gsub(".out.*","",rnanames)
# rnanames <- gsub(".*Sample_","",rnanames)
# colnames(rna) <- c("scaffold","pos",rnanames)
# rna <- rna[rowSums(rna[,-c(1:2)])>5,]

# subl <- (rna[,2]>703897&rna[,2]<705556)|(rna[,2]>614173&rna[,2]<615712)
# plot(rna[subl,50+2],pch=20,cex=.2,ylim=c(0,45))
# points(rna[subl,51+2],pch=20,cex=.2,col="red")
# abline(v=c(641929,694091),col="red")
# plot(rna[subl,2],log((rna[subl,51+2]+1)/(rna[subl,53+2]+1),base=2),pch=20,cex=.2)


# plot(rowSums(dep[,grep("BI|F|SH|KC",colnames(dep))]),pch=20,cex=.2)
# par(mfrow=c(2,1),mar=c(2,2,1.5,0.5))
# plot(rowSums(dep[,grep("ER",colnames(dep))]),pch=20,cex=.2)
# plot(rowSums(dep[,grep("KC",colnames(dep))]),pch=20,cex=.2)

##deletion boundaries:
	#approximate break points:
	# NBH1 640000-710000
	# NBH2 688000-733000
	# ER1 620000-703000

	#no deletions:
	reg1 <- dep[,2]>500000&dep[,2]<615000
	reg2 <- dep[,2]>735000&dep[,2]<900000
	#er1
	er1 <- dep[,2]>620000&dep[,2]<703000
	#nhb1
	nbh1 <- dep[,2]>640000&dep[,2]<710000
	#nbh2
	nbh2 <- dep[,2]>688000&dep[,2]<733000

	#all three overlap
	alldel <- dep[,2]>688000&dep[,2]<703000
	#only nbh2
	nbh2only <- dep[,2]>711000&dep[,2]<734000
	#only nbh1+er1
	nbh1er1 <- dep[,2]>641000&dep[,2]<687000
	#only er1
	er1only <- dep[,2]>620000&dep[,2]<639000


	scalevec <- colSums(dep[reg1|reg2,-c(1,2)])/(sum(reg2)+sum(reg1))


#the approach I used on Scaffold9924 isn't working so well here. 
	#call them by eyeballing coverage plots below. 
copies_per_ind <- data.frame(cbind(n_copies=rep(2,384),whichdel=rep(0,384),pop=gsub("-.*","",colnames(gts)[-1])),stringsAsFactors=FALSE)
rownames(copies_per_ind) <- colnames(gts)[-1]

nbhcopy <- t(t(depsub[,nbh+2])/scalevec[nbh]) %>%
	(function(x){
		s <- colSums(x)
		o <- order(s,decreasing=TRUE)
		o
		}) %>%
	(gts[,nbh+1])[,.] %>%
	colnames() %>%
	cbind(
		.,
		c(rep(2,28),rep(1,18),rep(0,2)),
		c(rep(0,28),rep(2,7),rep(1,13))
		)

copies_per_ind[nbhcopy[,1],1] <- nbhcopy[,2]
copies_per_ind[nbhcopy[,1],2] <- nbhcopy[,3]

nyccopy <- t(t(depsub[,nyc+2])/scalevec[nyc]) %>%
	(function(x){
		s <- colSums(x)
		o <- order(s,decreasing=TRUE)
		o
		}) %>%
	(gts[,nyc+1])[,.] %>%
	colnames() %>%
	cbind(
		.,
		c(rep(2,24),rep(1,18),rep(0,1)),
		c(rep(0,24),rep(1,19))
		)

copies_per_ind[nyccopy[,1],1] <- nyccopy[,2]
copies_per_ind[nyccopy[,1],2] <- nyccopy[,3]

shcopy <- t(t(depsub[,sh+2])/scalevec[sh]) %>%
	(function(x){
		s <- colSums(x)
		o <- order(s,decreasing=TRUE)
		o
		}) %>%
	(gts[,sh+1])[,.] %>%
	colnames() %>%
	cbind(
		.,
		c(rep(2,47),rep(1,3)),
		c(rep(0,47),rep(1,3))
		)

copies_per_ind[shcopy[,1],1] <- shcopy[,2]
copies_per_ind[shcopy[,1],2] <- shcopy[,3]

ercopy <- t(t(depsub[,er+2])/scalevec[er]) %>%
	(function(x){
		s <- colSums(x)
		o <- order(s,decreasing=TRUE)
		o
		}) %>%
	(gts[,er+1])[,.] %>%
	colnames() %>%
	cbind(
		.,
		c(rep(2,3),rep(1,13),rep(0,33)),
		c(rep(0,3),rep(3,46))
		)

copies_per_ind[ercopy[,1],1] <- ercopy[,2]
copies_per_ind[ercopy[,1],2] <- ercopy[,3]

copies_per_ind[,1] <- as.numeric(copies_per_ind[,1])
copies_per_ind[,2] <- as.numeric(copies_per_ind[,2])

table(copies_per_ind[,c(1,3)]) %>% 
	(function(x){x[,c(6,1,2,4,7,8,3,5)]}) %>% 
	(function(x){barplot(x,beside=TRUE)})


# ###assign NBH individuals to haplogroups
# subl <- hap[,1]>850000&hap[,1]<900000
# subi <- gsub("_.","",colnames(hap))%in%rownames(copies_per_ind[copies_per_ind[,2]>=0&grepl("NBH|BI",copies_per_ind[,3]),])
# dist(t(hap[subl,subi]),method="manhattan") %>% hclust() %>% plot(.,cex=.75)
# hap[subl,subi] %>%
# 	t() %>%
# 	dist(.,method="manhattan") %>%
# 	(function(x){ hclust(x)$order}) %>%
# 	(hap[,subi])[,.] %>%
# 	cbind(0,.) %>%
# 	pimage(.)
# 
# hapgrps <- hap[subl,subi] %>% t() %>% dist(.,method="manhattan") %>%
# 	(function(x){ 
# 		x <- hclust(x) 
# 		x <- cutree(x,5)
# 		x <- cbind(sam=gsub("_.","",names(x)),grp=x)
# 		x
# 		}) 
# hapgrps2 <- hapgrps %>% data.frame() %>% group_by(.,sam) %>% summarize(.,haps=paste(sort(grp),sep="",collapse="_")) %>% data.frame()
# rownames(hapgrps2) <- hapgrps2[,1]
# 


# return allele frequency given population, sex, data frame containing genotypes
freq <- function(df,pop,sex="F"){
	vec <- c(FALSE,grepl(sex,sexscore$sex))&grepl(pop,colnames(df))
	apply(df[,vec],
		MAR=1,
		FUN=function(x){
			sum(x,na.rm=TRUE)/sum(!is.na(x))/2
			}
		)
	}



#function to subsample and smooth a vector. in this case depth
	#can't use loess or anything complicated because it takes annoyingly long
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

#smooth and subsample depth
depsub <- as.data.frame(apply(dep[,-c(1,2)],MAR=2,FUN=subsmooth,by=20,width=500))
depsub <- cbind("Chr1sub",dep[seq(from=20,to=dim(dep)[1],by=20),2],depsub)
colnames(depsub)[1:2] <- c("scaffold","pos")

par(mfrow=c(2,1),mar=c(2,2,1.5,0.5))

bi <- grep("BI",sexscore$sample)
t(t(depsub[,bi+2])/scalevec[bi]) %>%
	(function(x){
		s <- colSums(x)
		o <- order(s,decreasing=TRUE)
		x[,o]
		}) %>%

	image(z=.,zlim=c(0,4),col=viridis(30),x=depsub[,2],y=(1:length(bi)))

nbh <- grep("NBH",sexscore$sample)
t(t(depsub[,nbh+2])/scalevec[nbh]) %>%
	(function(x){
		s <- colSums(x)
		o <- order(s,decreasing=TRUE)
		x[,o]
		}) %>%

	image(.,zlim=c(0,4),col=viridis(30),x=depsub[,2],y=(1:length(nbh)))

bi <- grep("BI",sexscore$sample)
t(t(depsub[,bi+2])/scalevec[bi]) %>%
	(function(x){
		s <- colSums(x)
		o <- order(s,decreasing=TRUE)
		o
		}) %>%
	(gts[,bi+1])[,.] %>%

	pimage(.)

nbh <- grep("NBH",sexscore$sample)
t(t(depsub[,nbh+2])/scalevec[nbh]) %>%
	(function(x){
		s <- colSums(x)
		o <- order(s,decreasing=TRUE)
		o
		}) %>%
	(gts[,nbh+1])[,.] %>%

	pimage(.)


bi <- grep("BI",sexscore$sample)
gts[,bi+1] %>%
	cor(.,use="pairwise") %>%
	(function(x){(1 / x)}) %>%
	as.dist() %>%
	(function(x){ hclust(x)$order}) %>%
	(gts[,bi+1])[,.] %>%

	pimage(.)

nbh <- grep("NBH",sexscore$sample)
gts[,nbh+1] %>%
	cor(.,use="pairwise") %>%
	(function(x){(1 / x)}) %>%
	as.dist() %>%
	(function(x){ hclust(x)$order}) %>%
	(gts[,nbh+1])[,.] %>%

	pimage(.)

bi <- grep("BI",colnames(hap))
hap[,bi] %>%
	t() %>%
	dist(.,method="manhattan") %>%
	(function(x){ hclust(x)$order}) %>%
	(hap[,bi])[,.] %>%

	pimage(.)

nbh <- grep("NBH",colnames(hap))
hap[subl,nbh] %>%
	t() %>%
	dist(.,method="manhattan") %>%
	(function(x){ hclust(x)$order}) %>%
	(hap[,nbh])[,.] %>%

	pimage(.)


f <- grep("F",sexscore$sample)
t(t(depsub[,f+2])/scalevec[f]) %>%
	(function(x){
		s <- colSums(x)
		o <- order(s,decreasing=TRUE)
		x[,o]
		}) %>%

	image(z=.,zlim=c(0,4),col=viridis(30),x=depsub[,2],y=(1:length(f)))

bp <- grep("BP",sexscore$sample)
t(t(depsub[,bp+2])/scalevec[bp]) %>%
	(function(x){
		s <- colSums(x)
		o <- order(s,decreasing=TRUE)
		x[,o]
		}) %>%

	image(z=.,zlim=c(0,4),col=viridis(30),x=depsub[,2],y=(1:length(bp)))

sh <- grep("SH",sexscore$sample)
t(t(depsub[,sh+2])/scalevec[sh]) %>%
	(function(x){
		s <- colSums(x)
		o <- order(s,decreasing=TRUE)
		x[,o]
		}) %>%

	image(z=.,zlim=c(0,4),col=viridis(30),x=depsub[,2],y=(1:length(sh)))

nyc <- grep("NYC",sexscore$sample)
t(t(depsub[,nyc+2])/scalevec[nyc]) %>%
	(function(x){
		s <- colSums(x)
		o <- order(s,decreasing=TRUE)
		x[,o]
		}) %>%

	image(z=.,zlim=c(0,4),col=viridis(30),x=depsub[,2],y=(1:length(nyc)))



nbh <- grep("NBH",sexscore$sample)
t(t(depsub[,nbh+2])/scalevec[nbh]) %>%
	(function(x){
		s <- colSums(x)
		o <- order(s,decreasing=TRUE)
		x[,o]
		}) %>%

	image(.,zlim=c(0,4),col=viridis(30),x=depsub[,2],y=(1:length(nbh)))

nbh <- grep("NBH",sexscore$sample)
t(t(depsub[,nbh+2])/scalevec[nbh]) %>%
	(function(x){
		s <- colSums(x)
		o <- order(s,decreasing=TRUE)
		x[,o]
		}) %>%

	image(.,zlim=c(0,4),col=viridis(30),x=depsub[,2],y=(1:length(nbh)))

kc <- grep("KC",sexscore$sample)
t(t(depsub[,kc+2])/scalevec[kc]) %>%
	(function(x){
		s <- colSums(x)
		o <- order(s,decreasing=TRUE)
		x[,o]
		}) %>%

	image(z=.,zlim=c(0,4),col=viridis(30),x=depsub[,2],y=(1:length(kc)))

er <- grep("ER",sexscore$sample)
t(t(depsub[,er+2])/scalevec[er]) %>%
	(function(x){
		s <- colSums(x)
		o <- order(s,decreasing=TRUE)
		x[,o]
		}) %>%

	image(.,zlim=c(0,4),col=viridis(30),x=depsub[,2],y=(1:length(er)))


kc <- grep("KC",sexscore$sample)
t(t(depsub[,kc+2])/scalevec[kc]) %>%
	(function(x){
		s <- colSums(x)
		o <- order(s,decreasing=TRUE)
		o
		}) %>%
	(gts[,kc+1])[,.] %>%

	pimage(.)

er <- grep("ER",sexscore$sample)
t(t(depsub[,er+2])/scalevec[er]) %>%
	(function(x){
		s <- colSums(x)
		o <- order(s,decreasing=TRUE)
		o
		}) %>%
	(gts[,er+1])[,.] %>%

	pimage(.)

#make representative coverage figures for three individuals, two duplications (large and smallest) and no duplication
	#281 - no deletion
	#253 - homozygous NBH1
	#284 - heterozygous NBH1
	#292,265 - heterozygous NBH2
	#148 - homozygous ER1
	par(mfrow=c(3,1),mar=c(2,2,1.5,0.5))
	ind <- 281
	subl <- seq(from=10,to=319910,by=10)
	plot(dep[subl,2],dep[subl,ind]/scalevec[ind-2],pch=20,cex=.2,ylim=c(0,3))
	abline(h=seq(from=0.5,to=2.5,by=0.5),col="darkgray",lwd=2)
	points(dep[subl,2],subsmooth(dep[,ind])/scalevec[ind-2],pch=20,cex=.2,col="red")
	ind <- 253
	subl <- seq(from=10,to=319910,by=10)
	plot(dep[subl,2],dep[subl,ind]/scalevec[ind-2],pch=20,cex=.2,ylim=c(0,3))
	abline(h=seq(from=0.5,to=2.5,by=0.5),col="gray",lwd=3)
	points(dep[subl,2],subsmooth(dep[,ind])/scalevec[ind-2],pch=20,cex=.2,col="red")
	ind <- 284
	subl <- seq(from=10,to=319910,by=10)
	plot(dep[subl,2],dep[subl,ind]/scalevec[ind-2],pch=20,cex=.2,ylim=c(0,3))
	abline(h=seq(from=0.5,to=2.5,by=0.5),col="gray",lwd=3)
	points(dep[subl,2],subsmooth(dep[,ind])/scalevec[ind-2],pch=20,cex=.2,col="red")
	ind <- 287
	subl <- seq(from=10,to=319910,by=10)
	plot(dep[subl,2],dep[subl,ind]/scalevec[ind-2],pch=20,cex=.2,ylim=c(0,3))
	abline(h=seq(from=0.5,to=2.5,by=0.5),col="gray",lwd=3)
	points(dep[subl,2],subsmooth(dep[,ind])/scalevec[ind-2],pch=20,cex=.2,col="red")
	ind <- 148
	subl <- seq(from=10,to=319910,by=10)
	plot(dep[subl,2],dep[subl,ind]/scalevec[ind-2],pch=20,cex=.2,ylim=c(0,3))
	abline(h=seq(from=0.5,to=2.5,by=0.5),col="gray",lwd=3)
	points(dep[subl,2],subsmooth(dep[,ind])/scalevec[ind-2],pch=20,cex=.2,col="red")


plot(
	(1.0*colSums(dep[nbh1,-c(1,2)]/sum(nbh1)) / 
	colSums(dep[reg1|reg2,-c(1,2)]/(sum(reg2)+sum(reg1))))[popord],
	pch=(sexscore[,3]=="M")[popord]+15,col=mat[pops[popord],5]
	)
abline(h=seq(from=0,to=6,by=0.5))


plot(
	1.0*colSums(dep[nbh3,-c(1,2)]/sum(nbh3)) / 
	colSums(dep[reg1|reg2,-c(1,2)]/(sum(reg2)+sum(reg1))),
	pch=(sexscore[,3]=="M")+15,col=mat[pops,5]
	)
abline(h=seq(from=0,to=6,by=0.5))

plot(
	1.0*colSums(dep[nbh2,-c(1,2)]/sum(nbh2)) / 
	colSums(dep[reg1|reg2,-c(1,2)]/(sum(reg2)+sum(reg1))),
	pch=(sexscore[,3]=="M")+15,col=mat[pops,5]
	)
abline(h=seq(from=0,to=6,by=0.5))

plot(
	1.0*colSums(dep[er1,-c(1,2)]/sum(er1)) / 
	colSums(dep[reg1|reg2,-c(1,2)]/(sum(reg2)+sum(reg1))),
	pch=(sexscore[,3]=="M")+15,col=mat[pops,5]
	)
abline(h=seq(from=0,to=6,by=0.5))

plot(
	1.0*colSums(dep[nbh1er1,-c(1,2)]/sum(nbh1er1)) / 
	colSums(dep[reg1|reg2,-c(1,2)]/(sum(reg2)+sum(reg1))),
	pch=(sexscore[,3]=="M")+15,col=mat[pops,5]
	)
abline(h=seq(from=0,to=6,by=0.5))





#####below is some unix code to generate MDS plots of this scaffold

tabix -h hallsnps.vcf.gz Scaffold10074:400000-610000 | \
~/bin/vcftools/bin/vcftools --vcf - --recode --stdout --keep er_kc.keep | \
gzip -c >Scaffold10074.erkc.vcf.gz
#gzip Scaffold10074.binbh.vcf

~/bin/plink_linux_x86_64/plink \
-vcf ~/popgen/variants/bowfree/ALL1/Scaffold10074.erkc.vcf.gz \
--allow-extra-chr \
--cluster \
--mds-plot 5 \
--out scaffold10074.erkc

mds <- read.table("scaffold10074.mds",stringsAsFactors=FALSE,header=TRUE)
mdsmost <- read.table("scaffold10074.most.mds",stringsAsFactors=FALSE,header=TRUE)
mdsnorth <- read.table("scaffold10074.north.mds",stringsAsFactors=FALSE,header=TRUE)
mdsbinbh <- read.table("scaffold10074.binbh.mds",stringsAsFactors=FALSE,header=TRUE)
mdsbpf <- read.table("scaffold10074.bpf.mds",stringsAsFactors=FALSE,header=TRUE)
mdsshnyc <- read.table("scaffold10074.shnyc.mds",stringsAsFactors=FALSE,header=TRUE)
mdserkc <- read.table("scaffold10074.erkc.mds",stringsAsFactors=FALSE,header=TRUE)

par(mfrow=c(1,3),mar=c(2,2,1.5,0.5))
tab <- mdsbinbh
i <- 6 ; j <- 7
plot(tab[,i],tab[,j],pch=(sexscore$sex=="M")+15,col=mat[gsub("-.*","",tab[,1]),5],cex=1.5)
text(x=tab[,i],y=tab[,j],labels=hapgrps2[tab[,1],2],adj=0)
plot(tab[,i],tab[,j],pch=(sexscore$sex=="M")+15,col=viridis(3)[copies_per_ind[tab[,1],1]+1],cex=1.5)
plot(tab[,i],tab[,j],pch=(sexscore$sex=="M")+15,col=copies_per_ind[tab[,1],2]+1,cex=1.5)

par(mfrow=c(1,3),mar=c(2,2,1.5,0.5))
tab <- mdsmost
i <- 4 ; j <- 5
plot(tab[,i],tab[,j],pch=(sexscore$sex=="M")+15,col=mat[gsub("-.*","",tab[,1]),5],cex=1.5)
#text(x=tab[,i],y=tab[,j],labels=hapgrps2[tab[,1],2],adj=0)
plot(tab[,i],tab[,j],pch=(sexscore$sex=="M")+15,col=viridis(3)[copies_per_ind[tab[,1],1]+1],cex=1.5)
plot(tab[,i],tab[,j],pch=(sexscore$sex=="M")+15,col=copies_per_ind[tab[,1],2]+1,cex=1.5)


subi <- rownames(copies_per_ind)[copies_per_ind[,1]<2&copies_per_ind[,3]=="NYC"]
subi2 <- rownames(copies_per_ind)[copies_per_ind[,1]==2&copies_per_ind[,3]=="NYC"]
del1 <- rowSums(gts[,subi],na.rm=TRUE)/rowSums(!is.na(gts[,subi]))/2
nond <- rowSums(gts[,subi2],na.rm=TRUE)/rowSums(!is.na(gts[,subi2]))/2



###make a tree!
Fpopfreqs <- cbind(
	freq(gts,"BI",sex="M|F"),
	freq(gts,"NBH",sex="M|F"),
	freq(gts,"F",sex="M|F"),
	freq(gts,"BP",sex="M|F"),
	freq(gts,"SH",sex="M|F"),
	freq(gts,"NYC",sex="M|F"),
	freq(gts,"KC",sex="M|F"),
	freq(gts,"ER",sex="M|F")
	)
colnames(Fpopfreqs) <- c("BI","NBH","F","BP","SH","NYC","KC","ER")
Fpopfreqs <- Fpopfreqs[which(rowSums(is.na(Fpopfreqs))==0),]

#using contml
tr <- Rcontml_quieter(X=t(Fpopfreqs),path="~/bin/phylip-3.696/exe/")

#already ran bootstraps, don't re-run unless necessary. 
trboot <- boot.phylo(tr,t(Fpopfreqs),FUN=function(x){Rcontml_quieter(X=x,path="~/bin/phylip-3.696/exe/")},trees=TRUE)

#for some reason boot.phylo calculates BP wrong. 
	#table of BP
	as.prop.part(bitsplits(trboot$trees))
	#replace with proper BP:
	trboot$BP <- c(100,100,100,100,100,100)

plot(tr,type="unrooted")
nodelabels(trboot$BP)

trN <- Rgendist(t(Fpopfreqs),path="~/bin/phylip-3.696/exe/") %>% 
	Rneighbor(.,path="~/bin/phylip-3.696/exe/")


trbootN <- boot.phylo(trN,t(Fpopfreqs),FUN=function(x){
	Rgendist(t(Fpopfreqs),path="~/bin/phylip-3.696/exe/") %>% 
		Rneighbor(.,path="~/bin/phylip-3.696/exe/")
	},trees=TRUE)

#for some reason boot.phylo calculates BP wrong. 
	#table of BP
	as.prop.part(bitsplits(trbootN$trees))
	#replace with proper BP:
	trboot$BP <- c(100,100,100,100,100,100)








