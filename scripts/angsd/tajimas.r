#tajima's d from pi, tw, and n samples.
#pi and tw assumed to be summed, NOT per site

tajimas<-function(pi,tw,n){

	if(any(is.na(c(pi,tw)))){
		return(NA)
		}

	nm<-n-1
	a1<-sum(1/1:nm)
	a2<-sum(1/(1:nm)^2)
	b1<-(n+1)/(3*(n-1))
	b2<-(2*((n^2)+n+3))/(9*n*(n-1))
	c1<-b1-(1/a1)
	c2<-b2-((n+2)/(a1*n))+a2/(a1^2)
	e1<-c1/a1
	e2<-c2/((a1^2)+a2)



	S<-tw*a1
	num<-pi-tw
	denom<-sqrt((e1*S)+(e2*S*(S-1)))

	tajd<-num/denom

	return(tajd)

	}
