out<-read.table("~/analysis/data/fst/individual_pbs/overlaps",stringsAsFactors = FALSE)

bb<-out[,6]>0 & out[,7]<1 & out[,8]<1 & out[,9]<1 & out[,10]<1 
vb<-out[,7]>0 & out[,6]<1 & out[,8]<1 & out[,9]<1 & out[,10]<1 
pb<-out[,8]>0 & out[,7]<1 & out[,6]<1 & out[,9]<1 & out[,10]<1 
sj<-out[,9]>0 & out[,7]<1 & out[,8]<1 & out[,6]<1 & out[,10]<1 
bnp<-out[,10]>0 & out[,7]<1 & out[,8]<1 & out[,9]<1 & out[,6]<1 

res<-out[,6]>0 & out[,7]>0 & out[,8]>0 & out[,9]<1 & out[,10]<1
interm<-out[,6]<1 & out[,7]<1 & out[,8]<1 & out[,9]>0 & out[,10]>0
