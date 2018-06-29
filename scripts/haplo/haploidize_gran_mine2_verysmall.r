#!/usr/bin/env Rscript

library(stringr)
library(magrittr)

sam<-function(x){
      y <- str_split(x,":")[[1]][1] %>% 
    			 unlist() %>% 
    			 as.numeric() 
        return(y)
      }

f <- file("stdin")
#f<-file("~/analysis/data/dfst/outliers/zshared.vcf.bgz")
open(f)
while(length(line <- readLines(f,n=1)) > 0) {
  if(grepl("^#",line)){write(line,stdout());next()}

  line <- str_split(line,"\\t") %>% unlist()
  line[10:106] <- sapply(line[10:106],sam)
  line <- paste(line,collapse="\t")
  write(line,stdout())
  # write(line, stderr())
  # process line
}



