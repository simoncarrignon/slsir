source("abmEpi.R")


poptest=generatePopulation(500) 
neutral=replicate(50,abmSIR(poptest,2500,speed=.6,p=c(1,1),di=2,i0=1,visu=F)$timeseries[,2])
twiceless=replicate(50,abmSIR(poptest,2500,speed=.6,p=c(1,.5),di=2,i0=1,visu=F,inf=.5)$timeseries[,2])
tenless=replicate(50,abmSIR(poptest,2500,speed=.6,p=c(1,.1),di=2,i0=1,visu=F,inf=.5)$timeseries[,2])

neutral=apply(replicate(50,abmSIR(poptest,1500,speed=.6,p=c(1,1),di=2,i0=1,visu=F)$timeseries[,2]),1,mean)
twiceless=apply(replicate(50,abmSIR(poptest,1500,speed=.6,p=c(1,.5),di=2,i0=1,visu=F,inf=.5)$timeseries[,2]),1,mean)
tenless=apply(replicate(50,abmSIR(poptest,1500,speed=.6,p=c(1,.1),di=2,i0=1,visu=F,inf=.5)$timeseries[,2]),1,mean)

inpoint=seq(0,1,.1)
library(parallel)
cl <- makeForkCluster(40,outfile="")
twiceless=sapply(inpoint,function(i)apply(parSapply(cl,1:300,function(j){print(j);abmSIR(poptest,1500,speed=.6,p=c(1,.5),di=2,i0=1,visu=F,inf=i,sat=50)}$timeseries[,2]),1,mean))
tenless=sapply(inpoint,function(i)apply(parSapply(cl,1:300,function(j){print(j);abmSIR(poptest,1500,speed=.6,p=c(1,.1),di=2,i0=1,visu=F,inf=i,sat=50)}$timeseries[,2]),1,mean))
stopCluster(cl)

cols=heat.colors(length(inpoint))
names(cols)=as.character(inpoint)
plot(neutral,type="l",xlim=c(0,1500))
sapply(1:length(inpoint),function(i)lines(twiceless[,i],col=cols[as.character(inpoint[i])]))
sapply(1:length(inpoint),function(i)lines(tenless[,i],col=cols[as.character(inpoint[i])]))
apply(tenless,2,lines,col="red")

