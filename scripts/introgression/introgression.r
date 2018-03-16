library(viridis)
library(dplyr)
library(magrittr)
library(Rphylip)
library(ape)
library(stringr)

load("~/analysis/scripts/introgression/AHR_elias_plotdata.Rdata")
names<-row.names(mds)
popcol<-ifelse(grepl("-",names),"plum",
               ifelse(grepl("BB",names),"black",
                      ifelse(grepl("VB",names),"black",
                             ifelse(grepl("PB",names),"black",
                                    ifelse(grepl("SJ",names),"firebrick2",
                                           ifelse(grepl("BNP",names),"firebrick2",
                                                  ifelse(grepl("SP",names),"cadetblue3",
                                                         ifelse(grepl("GB",names),"cadetblue3","WRONG"))))))))
popcol2<-ifelse(grepl("-",names),"plum",
               ifelse(grepl("BB",names),"black",
                      ifelse(grepl("VB",names),"grey40",
                             ifelse(grepl("PB",names),"grey80",
                                    ifelse(grepl("SJ",names),"firebrick2",
                                           ifelse(grepl("BNP",names),"lightpink",
                                                  ifelse(grepl("SP",names),"cadetblue3",
                                                         ifelse(grepl("GB",names),"cadetblue1","WRONG"))))))))               
popsym<-ifelse(grepl("d2",pop2),"23",
               ifelse(grepl("d1",pop2),"24",
                      ifelse(grepl("d0",pop2),"21","22")))

#plotting NJ tree
plot(tr,type="unrooted",show.tip.label=FALSE)
tiplabels(pch=as.numeric(popsym),col=popcol,bg=popcol2,cex=1.5)

#plotting MDS
plot(mds,col=popcol,cex=2,cex.axis=2,pch=as.numeric(popsym),
     bg=popcol2,bty='n')

#Plotting with ER standing out
#popcol<-ifelse(grepl("ER",names),"seagreen1",
#                ifelse(grepl("BB",names),"black",
#                       ifelse(grepl("VB",names),"black",
#                              ifelse(grepl("PB",names),"black",
#                                     ifelse(grepl("SJ",names),"firebrick2",
#                                            ifelse(grepl("BNP",names),"firebrick2",
#                                                   ifelse(grepl("SP",names),"cadetblue3",
#                                                          ifelse(grepl("GB",names),"cadetblue3","plum"))))))))
# popcol2<-ifelse(grepl("ER",names),"seagreen1",
#                 ifelse(grepl("BB",names),"black",
#                        ifelse(grepl("VB",names),"grey40",
#                               ifelse(grepl("PB",names),"grey80",
#                                      ifelse(grepl("SJ",names),"firebrick2",
#                                             ifelse(grepl("BNP",names),"lightpink",
#                                                    ifelse(grepl("SP",names),"cadetblue3",
#                                                           ifelse(grepl("GB",names),"cadetblue1","plum"))))))))     
