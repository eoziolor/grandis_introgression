#!/usr/bin/env Rscript

library(stringr)
library(magrittr)


	nal <- function(x){
		x <- str_split(x,",") %>% unlist() %>% length()
		x + 1
		}

	sam<-function(x,len){
		if(x=="."){return(".")}

		x <- str_split(x,":")[[1]][3] %>% 
			str_split(.,",") %>% 
			unlist() %>% 
			as.numeric() # pull allele depth TAG field
		x <- x!=0 # alleles present/absent
		if(sum(x)==0){return(".")} # no genotype if no alleles
		x <- x/sum(x) # probability vector
		x <- sample(x=len,size=1,prob=x) - 1
		return(x)
		}


#f <- file("stdin")
f<-file("~/analysis/data/dfst/outliers/zshared.vcf.bgz")
open(f)
while(length(line <- readLines(f,n=1)) > 0) {
  if(grepl("^#",line)){write(line,stdout());next()}

  line <- str_split(line,"\\t") %>% unlist()
  nl <- nal(line[5])
  line[10:297] <- sapply(line[10:297],sam,len=nl)
  line <- paste(line,collapse="\t")
  write(line,stdout())
  # write(line, stderr())
  # process line
}



