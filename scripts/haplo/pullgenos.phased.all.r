pullgenos.phased <- function(vcf){

	cline<-paste("tabix ", vcf, " ",sep="")
	header<-scan(pipe(paste("zcat",vcf, "| grep -v '##' | head -n 1",sep=" ")),what="character")
	gts<-read.table(pipe(cline),stringsAsFactors=FALSE)

	colnames(gts)<-header
	gts2<-as.matrix(gts[,10:297])
	gts2<-gsub(":.*","",gts2)
	gts3 <- c()
	for(i in 1:dim(gts2)[2]){
		temp <- strsplit(gts2[,i],"\\/")
		temp <- do.call(rbind,temp)
		colnames(temp) <- paste(colnames(gts2)[i],c("_a","_b"),sep="")
		gts3 <- cbind(gts3,temp)
		}
	gts3<-cbind(gts[,2],gts3)
	class(gts3)<-"numeric"
	return(gts3)

}
