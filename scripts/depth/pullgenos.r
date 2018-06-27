pullgenos<-function(){

	cline<-paste("tabix ", vcf, " ", scaf,sep="")
	header<-scan(pipe(paste("zcat",vcf, "| grep -v '##' | head -n 1",sep=" ")),what="character")
	gts<-read.table(pipe(cline),stringsAsFactors=FALSE)

	colnames(gts)<-header
	gts2<-as.matrix(gts[,10:585])
	gts2<-gsub(":.*","",gts2)
	gts2[gts2=="."]<-NA
	gts2[gts2=="0/0"]<-0
	gts2[gts2=="0/1"]<-1
	gts2[gts2=="1/1"]<-2

	gts2<-cbind(gts[,2],gts2)
	class(gts2)<-"numeric"
	return(gts2)
	}

