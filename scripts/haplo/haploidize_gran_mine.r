#!/usr/bin/env Rscript

library(stringr)
library(magrittr)

sam<-function(x){
      y <- str_split(x,":")[[1]][1] %>% 
    			 str_split(.,"\\/") %>% 
    			 unlist() %>% 
    			 as.numeric() # pull genotypes
      if(is.na(y)){
        return(".")
      } else {
          z<-str_split(x,":")[[1]][3] %>% 
             str_split(.,",") %>% 
             unlist() %>% 
             as.numeric() # pull allele coverage
          if(y[1]==y[2]){
            h<-y[1]
          } else {
            if(z[1]>z[2]){
              h<-y[1]
            } else {
              h<-sample(y,size=1)
            }
          }
        return(h)
      }
}

f <- file("stdin")
#f<-file("~/analysis/data/dfst/outliers/zshared.vcf.bgz")
open(f)
while(length(line <- readLines(f,n=1)) > 0) {
  if(grepl("^#",line)){write(line,stdout());next()}

  line <- str_split(line,"\\t") %>% unlist()
  line[10:297] <- sapply(line[10:297],sam)
  line <- paste(line,collapse="\t")
  write(line,stdout())
  # write(line, stderr())
  # process line
}



