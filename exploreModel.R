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


simuDimanche<-function(){
poptest=generatePopulation(500,recovery=2500)

	probas=seq(0.1,1,.1)
	inpoint=c(0.1,.5,.9)

## generate all simu
	allres=lapply(probas,function(p)lapply(inpoint,function(i)parLapply(cl,1:300,function(j){print(j);abmSIR(poptest,1500,speed=.6,p=c(1,p),di=2,i0=1,visu=F,inf=i,sat=20)}$timeseries[,2])))

	clrsprobs=colorRampPalette(c("blue","red"))(length(probas))
	clrsinp=colorRampPalette(c("green","yellow"))(length(inpoint))


#apply summary
	hein=lapply(allres,lapply,function(i)do.call("cbind",i))
	means=lapply(hein,lapply,apply,1,mean)
	neutral=do.call("cbind",neutral)
	neutral=apply(neutral,1,mean)

neutral=lapply(1:300,function(j){print(j);abmSIR(poptest,1500,speed=.6,p=c(1,1),di=2,i0=1,visu=F,inf=1,sat=20)}$timeseries[,2])

#Visualize

pdf("twofirstdim.pdf",width=14,height=5)
	par(mfrow=c(2,3))
	par(mar=c(4,4,1,1))
	for(i in 1:length(inpoint)){
		plot(1:1500,neutral,type="l",xlim=c(1,1200),ylim=c(0,500),main=paste("inpoint=",inpoint[i]),ylab="#infected",xlab="",lwd=2)
		for(u in 1:length(probas)){
			lines(1:1500,means[[u]][[i]],lty=1,col=clrsprobs[u])
		}
		legend("bottomright",legend=c("neutral",paste0("PiB=",1/probas[c(1,2,4,10)],"PiG")),lty=c(1,rep(1,4)),col=c(1,clrsprobs[c(1,2,4,10)]),lwd=2)
	}
	par(mar=c(4,4,1,1))
	for(i in (1:length(probas))[c(1,5,9)]){
		plot(neutral,type="l",xlim=c(1,1200),ylim=c(0,500),main=bquote(P*i[B]==.(1/probas[i])*P*i[G]),xlab="#time",ylab="#infected",lwd=2)
		for(u in 1:length(inpoint)){
			lines(1:1500,means[[i]][[u]],lty=1,col=clrsinp[u])
		}
		legend("bottomright",legend=c("neutral",paste0("inpoininpoint=",inpoint)),lty=c(1,rep(1,3)),col=c(1,clrsinp),lwd=2)
	}
dev.off()

	p=.95
	n=500
	countmin=lapply(allres,lapply,lapply,function(u)min(which(u>(p*n))))
	countmin=lapply(countmin,lapply,unlist)
pdf("heamtapPiG.pdf")
	image(x=probas,y=inpoint,z=t(sapply(countmin,sapply,mean)),xlab=expression(P*i[G]),main="time to infect 90% of the population")
dev.off()
}



simuWithRecoverTime <- function(i){
    xsize=ysize=100
    poptest=generatePopulation(500,recovery=c(8,14)*25,speed=c(1,.2),xsize=xsize,ysize=ysize) 
    a=abmSIR(poptest,1500,p=c(1,1),di=2,i0=1,inf=1,sat=20,xsize=xsize,ysize=ysize,visu=T)
    neutral=lapply(1:100,function(j){print(j);abmSIR(poptest,1000,p=c(1,1),di=2,i0=1,visu=F,inf=1,sat=20,xsize=xsize,ysize=ysize)}$timeseries[,2])
    baseline=mean(lapply(neutral,max))

}
