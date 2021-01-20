
png("5e8779c56517890001536101/figures/fullTraj.png",width=800,height=400,pointsize=17)
par(mfrow=c(1,2))
plot(1,1,type="n",ylim=c(0,500),xlim=c(1,1500),main="No learning")
legend("topright",legend=c("All bad","All good"),col=c(myred,mygreen),lwd=2)
apply(neutralGood,2,function(i)lines(i,col=alpha("blue",.1),lwd=2))
lines(apply(neutralGood,1,median),col=alpha(myred,.9),lwd=2)
apply(neutralBad,2,function(i)lines(i,col=alpha(mygreen,.1),lwd=2))
lines(apply(neutralBad,1,median),col=alpha(mygreen,.9),lwd=2)
plot(1,1,type="n",ylim=c(0,500),xlim=c(1,1500),main="Learning")
apply(repetBest,2,function(i)lines(i,col=alpha(mygreen,.1),lwd=2))
lines(apply(repetBest,1,median),col=alpha(mygreen,.9),lwd=2,lty=2)
apply(repetWorst,2,function(i)lines(i,col=alpha(myred,.1),lwd=2))
lines(apply(repetWorst,1,median),col=alpha(myred,.9),lwd=2,lty=2)
legend("topright",legend=c("Revert","No revert"),col=c(myred,mygreen),lwd=2)
dev.off()


repetBest=sapply(1:100,function(i){print(i);slsirSimu(poptest,1500,p=c(1,.2),di=2,i0=1,inf=.2,sat=1,xsize=xsize,ysize=ysize,visu=F,ap=F,ts=T,file=F,p_i=.5,log=F)$timeseries[,2]})
repetBestDos=parSapply(cl,1:100,function(i){print(i);slsirSimu(poptest,1500,p=c(1,.2),di=2,i0=1,inf=.9,sat=1000,xsize=xsize,ysize=ysize,visu=F,ap=F,ts=T,file=F,p_i=.5,log=F)$timeseries[,2]})
neutralbad=sapply(1:10,function(i){print(i);slsirSimu(poptest,1500,p=c(1,1),di=2,i0=1,inf=9,sat=1000,xsize=xsize,ysize=ysize,visu=F,ap=F,ts=T,file=F,p_i=1,log=F)$timeseries[,2]})
library(parallel)

load("best.bin")
cl <- makeForkCluster(3,outfile="")
neutralGood=parSapply(cl,1:100,function(i){print(i);slsirSimu(poptest,1500,p=c(.1,.1),di=2,i0=1,inf=9,sat=1000,xsize=xsize,ysize=ysize,visu=F,ap=F,ts=T,file=F,p_i=1,log=F)$timeseries[,2]})
neutralBad=parSapply(cl,1:100,function(i){n=sample(nrow(best),1);print(paste(i,n));slsirSimu(poptest,1500,p=c(1,1),di=2,i0=1,inf=best$inf[n],pind=best$pind[n],inf_r=best$inf_r[n],sat_r=best$sat_r[n],sat=best$sat[n],xsize=xsize,ysize=ysize,visu=F,ap=F,ts=T,file=F,p_i=1,log=F)$timeseries[,2]})
repetBest=sapply(1:100,function(i){print(i);slsirSimu(poptest,1500,p=c(1,.2),di=2,i0=1,inf=.2,sat=1,xsize=xsize,ysize=ysize,visu=F,ap=F,ts=T,file=F,p_i=.5,log=F)$timeseries[,2]})
stopCluster(cl)

png("5e8779c56517890001536101/figures/threeDimensionDistance.png",width=1000,height=400,pointsize=20)
par(mfrow=c(1,3),cex=.5)
test=allresults
plot(test$time_max,test$max_infect,col=alpha(color.gradient(log(test$sat),c(mygreen,myred)),.2),main=expression(kappa),pch=20,xlab=expression(tau),ylab=expression(I[max]))
legend("topright",legend=paste("steepness=10^",round(log10(unique(test$sat)[c(1,5,10)]))),col=clrssat[c(1,5,10)],lty=1)
plot(test$time_max,test$max_infect,col=alpha(color.gradient(test$inf,c(mygreen,myred)),.2),main=expression(nu),pch=20,xlab=expression(tau),ylab=expression(I[max]))
legend("topright",legend=paste("point=",round(unique(test$inf)[c(1,5,10)],digit=1)),col=clrssat[c(1,5,10)],lty=1)
plot(test$time_max,test$max_infect,col=alpha(color.gradient(test$pind,c(mygreen,myred)),.2),main="indiv learning",pch=20,xlab=expression(tau),ylab=expression(I[max]))
legend("topright",legend=paste("P indiv learning=",unique(test$pind)),col=clrssat[c(1,5,10)],lty=1)
dev.off()


### Figure4a
png("5e8779c56517890001536101/figures/distribSimulationGreyZonedS30.png",pointsize=17, res = 100, width = 7, height = 7, units = "in")
plot(allresults$time_max,allresults$max_infect,col=alpha(color.gradient(allresults$distances,c(mygreen,"yellow", myred)),.2),main="distance",pch=20,xlab=expression(tau),ylab=expression(I[max]))
sset=sort(allresults$distances)
cols=alpha(color.gradient(sset,c(mygreen,"yellow", myred)),1)
lsset=seq(1,length(sset),length.out=4)
v=m(t=t,Rm=max(allresults$max_infect)-min(allresults$max_infect),mm=min(allresults$max_infect),Rt=max(allresults$time_max)-min(allresults$time_max),mt=min(allresults$time_max),s=.3)
t=100:1500
polygon(c(200,t,1510),c(0,v,0),col=alpha(colorbest,.5),lwd=2,angle=60,density=5)
legend("topright",legend=c(sapply(round(sset[lsset],digit=2),function(d)as.expression(bquote(delta==.(d)))),expression(delta<.25)),col=c(cols[lsset],NA),fill=c(rep(0,4),alpha(colorbest,.5)),border=c(rep(NA,4),alpha(colorbest,.5)),angle=c(rep(NA,4),60),density=c(rep(NA,4),25),pch=20)
#points(best1000$time_max,best1000$max_infect,col="green",main="distance",pch=20,xlab=expression(tau),ylab=expression(I[max]))
#points(best500$time_max,best500$max_infect,col="red",main="distance",pch=20,xlab=expression(tau),ylab=expression(I[max]))
#points(best250$time_max,best250$max_infect,col=colorbest,main="distance",pch=20,xlab=expression(tau),ylab=expression(I[max]))
dev.off()


### Figure4b
png("5e8779c56517890001536101/figures/distribSimulationGreyZonedS60.png",pointsize=17, res = 100, width = 7, height = 7, units = "in")
plot(allresults$time_max,allresults$max_infect,col=alpha(color.gradient(allresults$distances,c(mygreen,"yellow", myred)),.2),main="distance",pch=20,xlab=expression(tau),ylab=expression(I[max]))
sset=sort(allresults$distances)
cols=alpha(color.gradient(sset,c(mygreen,"yellow", myred)),1)
lsset=seq(1,length(sset),length.out=4)
v=m(t=t,Rm=max(allresults$max_infect)-min(allresults$max_infect),mm=min(allresults$max_infect),Rt=max(allresults$time_max)-min(allresults$time_max),mt=min(allresults$time_max),s=.6)
polygon(c(100,t,1010),c(0,v,0),col=alpha(colorbest,.5),lwd=2,angle=60,density=5)
legend("topright",legend=c(sapply(round(sset[lsset],digit=2),function(d)as.expression(bquote(delta==.(d)))),expression(delta<.6)),col=c(cols[lsset],NA),fill=c(rep(0,4),alpha(colorbest,.5)),border=c(rep(NA,4),alpha(colorbest,.5)),angle=c(rep(NA,4),60),density=c(rep(NA,4),25),pch=20)
dev.off()


    hdr.boxplot.2d(posterior$best$inf_r,log10(posterior$best$sat_r),prob=seq(75,90,10),shadecols=alpha(colorbest,.9),xlim=c(0,1),ylim=c(-1,3),xlab=expression(nu[r]),ylab="revert steepness (log10)")
    hdr.boxplot.2d(posterior$best$inf,log10(posterior$best$sat),prob=seq(75,90,10),shadecols=alpha(colorbest,.9),xlim=c(0,1),ylim=c(-1,3),xlab=expression(nu[r]),ylab="steepness (log10)")

    hdr.boxplot.2d(posterior$max_infect150$inf_r,log10(posterior$max_infect150$sat_r),prob=seq(20,100,10),shadecols=alpha(colorbest,.9),xlim=c(0,1),ylim=c(-1,3),xlab=expression(nu[r]),ylab="revert steepness (log10)")
    hdr.boxplot.2d(posterior$max_infect150$inf,log10(posterior$max_infect150$sat),prob=seq(20,100,10),shadecols=alpha(colorbest,.9),xlim=c(0,1),ylim=c(-1,3),xlab=expression(nu[r]),ylab="steepness (log10)")

    ##GRaph old results (to chekc maybe to find sense of p_ind)
######

load("allresults.bin")
allresults=allresults[-which(allresults$max_infect <4),]
allresults$distances=(1-bdistance(allresults$time_max) + bdistance(allresults$max_infect))/2

png("5e8779c56517890001536101/figures/exploringIndividualLearning.png",width=1000,height=1000,pointsize=25)
par(mfrow=c(3,3))
    clrssat=colorRampPalette(c("blue","red"))(10)
for(p_i in unique(allresults$pind)){
    test=allresults[allresults$pind == p_i ,]
    a=lapply(unique(test$sat),function(s){ sapply(unique(test$inf),function(i)quantile(test$distances[test$inf ==i & test$sat== s])[2:4])})
    b=lapply(unique(test$sat),function(s){ sapply(unique(test$inf),function(i)quantile(test$time_max[test$inf ==i & test$sat== s])[2:4])})
    c=lapply(unique(test$sat),function(s){ sapply(unique(test$inf),function(i)quantile(test$max_inf[test$inf ==i & test$sat== s])[2:4])})
    plotQuartiles(a,xlab="inp",ylab=expression(delta),main=paste("proba individual learning=",p_i))
    plotQuartiles(b,xlab="inp",ylab=expression(tau),main=paste("proba individual learning=",p_i))
    plotQuartiles(c,xlab="inp",ylab=expression(I[max]),main=paste("proba individual learning=",p_i))
    legend("topleft",legend=paste0("steep=10^",round(log10(unique(test$sat)[c(1,5,10)]))),col=clrssat[c(1,5,10)],lty=1,cex=.6)
}
dev.off()


