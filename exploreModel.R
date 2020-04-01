source("abmEpi.R")

poptest=generatePopulation(500,recovery=2500) 
neutral=replicate(50,abmSIR(poptest,2500,speed=.6,p=c(1,1),di=2,i0=1,visu=F)$timeseries[,2])
twiceless=replicate(50,abmSIR(poptest,2500,speed=.6,p=c(1,.5),di=2,i0=1,visu=F,inf=.5)$timeseries[,2])
tenless=replicate(50,abmSIR(poptest,2500,speed=.6,p=c(1,.1),di=2,i0=1,visu=F,inf=.5)$timeseries[,2])

neutral=apply(replicate(50,abmSIR(poptest,1500,speed=.6,p=c(1,1),di=2,i0=1,visu=F)$timeseries[,2]),1,mean)
twiceless=apply(replicate(50,abmSIR(poptest,1500,speed=.6,p=c(1,.5),di=2,i0=1,visu=F,inf=.5)$timeseries[,2]),1,mean)
tenless=apply(replicate(50,abmSIR(poptest,1500,speed=.6,p=c(1,.1),di=2,i0=1,visu=F,inf=.5)$timeseries[,2]),1,mean)

probas=seq(0.1,1,.1)
inpoint=c(0.1,.5,.9)
library(parallel)
cl <- makeForkCluster(40,outfile="")
all=lapply(probas,function(p)lapply(inpoint,function(i)parLapply(cl,1:300,function(j){print(j);abmSIR(poptest,1500,speed=.6,p=c(1,p),di=2,i0=1,visu=F,inf=i,sat=20)}$timeseries[,2])))
tenless=lapply(inpoint,function(i)apply(parSapply(cl,1:300,function(j){print(j);abmSIR(poptest,1500,speed=.6,p=c(1,.1),di=2,i0=1,visu=F,inf=i,sat=50)}$timeseries[,2]),1,quantile))
stopCluster(cl)

cols=heat.colors(length(inpoint))
names(cols)=as.character(inpoint)
plot(neutral,type="l",xlim=c(0,1500))
sapply(1:length(inpoint),function(i)lines(twiceless[,i],col=cols[as.character(inpoint[i])]))
sapply(1:length(inpoint),function(i)lines(tenless[,i],col=cols[as.character(inpoint[i])]))
apply(tenless,2,lines,col="red")

simuAndVisuFriday2703<-function(){
	poptest=generatePopulation(500)
	neutral=apply(replicate(50,abmSIR(poptest,1500,speed=.6,p=c(1,1),di=2,i0=1,visu=F)$timeseries[,2]),1,mean)
	inpoint=seq(0,1,.1)
	twiceless=sapply(inpoint,function(i)apply(sapply(1:30,function(j){print(j);abmSIR(poptest,1500,speed=.6,p=c(1,.5),di=2,i0=1,visu=F,inf=i,sat=50)}$timeseries[,2]),1,mean))
	tenless=sapply(inpoint,function(i)apply(sapply(1:30,function(j){print(j);abmSIR(poptest,1500,speed=.6,p=c(1,.1),di=2,i0=1,visu=F,inf=i,sat=50)}$timeseries[,2]),1,mean))

	cols=colorRampPalette(c("blue","red"))(length(inpoint))
	names(cols)=as.character(inpoint)
	plot(neutral,type="l",xlim=c(0,1500),main=expression(P*i[B] == 2*P*i[G]),ylab="#infected",xlab="time")
	sapply(1:length(inpoint),function(i)lines(twiceless[,i],col=cols[as.character(inpoint[i])],lwd=2))
	legend("bottomright",legend=c("neutral",paste("inflection=",inpoint[c(1,5,11)])),lty=c(1,rep(1,3)),col=c(1,cols[c(1,5,11)]),lwd=2)
	plot(neutral,type="l",xlim=c(0,1500),main=expression(P*i[B] == 10*P*i[G]),ylab="#infected",xlab="time")
	sapply(1:length(inpoint),function(i)lines(tenless[,i],col=cols[as.character(inpoint[i])],lwd=2,lty=2))
	legend("bottomright",legend=c("neutral",paste("inflection=",inpoint[c(1,5,11)])),lty=c(1,rep(2,3)),col=c(1,cols[c(1,5,11)]),lwd=2)
	#apply(tenless,2,lines,col="red")
}

