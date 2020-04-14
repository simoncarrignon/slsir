m<-function(t,Rm,Rt,mm,mt,s)(Rm*(2*s-1+(t-mt)/Rt))+mm 

source("abmEpi.R")

myred=rgb(r=213,g=94,b=0,maxColorValue=255)
mygreen=rgb(r=0,g=158,b=115,maxColorValue=255)



plotQuartiles <- function(allres,ylim=NULL,...){
    if(is.null(ylim))ylim=range(allres)
    plot(1,1,type="n",ylim=ylim,xlim=range(test$inf),...)
    lapply(1:length(allres),function(i){
           ai=allres[[i]]
           lines(unique(test$inf),ai[2,],col=clrssat[i],pch=20,type="b",cex=.8)
           arrows(unique(test$inf),ai[1,],unique(test$inf),ai[3,],col=clrssat[i],angle=90,code=3,length=.01,lwd=1)
})
}

color.gradient <- function(x, colors=c(myred,"yellow",mygreen), colsteps=100) {
    return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
}


bscore <- function(x)(x-min(x))/(max(x)-min(x))

load("allresults.bin")
allresults=allresults[-which(allresults$max_infect <4),]
allresults$scores=(1-bscore(allresults$time_max) + bscore(allresults$max_infect))/2

png("5e8779c56517890001536101/figures/exploringIndividualLearning.png",width=1000,height=1000,pointsize=25)
par(mfrow=c(3,3))
for(p_i in unique(allresults$pind)){
    test=allresults[allresults$pind == p_i ,]
    a=lapply(unique(test$sat),function(s){ sapply(unique(test$inf),function(i)quantile(test$scores[test$inf ==i & test$sat== s])[2:4])})
    b=lapply(unique(test$sat),function(s){ sapply(unique(test$inf),function(i)quantile(test$time_max[test$inf ==i & test$sat== s])[2:4])})
    c=lapply(unique(test$sat),function(s){ sapply(unique(test$inf),function(i)quantile(test$max_inf[test$inf ==i & test$sat== s])[2:4])})
    plotQuartiles(a,xlab="inp",ylab="score",main=paste("proba individual learning=",p_i))
    plotQuartiles(b,xlab="inp",ylab="time to max",main=paste("proba individual learning=",p_i))
    plotQuartiles(c,xlab="inp",ylab="max infected",main=paste("proba individual learning=",p_i))
    legend("topleft",legend=paste0("steep=10^",round(log10(unique(test$sat)[c(1,5,10)]))),col=clrssat[c(1,5,10)],lty=1,cex=.6)
}
dev.off()



png("5e8779c56517890001536101/figures/threeDimensionDistance.png",width=1000,height=400,pointsize=20)
par(mfrow=c(1,3),cex=.5)
test=allresults
plot(test$time_max,test$max_infect,col=alpha(color.gradient(log(test$sat),c(mygreen,myred)),.2),main="steepness",pch=20,xlab="Time Max",ylab="Max Infected")
legend("topright",legend=paste("steepness=10^",round(log10(unique(test$sat)[c(1,5,10)]))),col=clrssat[c(1,5,10)],lty=1)
plot(test$time_max,test$max_infect,col=alpha(color.gradient(test$inf,c(mygreen,myred)),.2),main="inflection point",pch=20,xlab="Time Max",ylab="Max Infected")
legend("topright",legend=paste("point=",round(unique(test$inf)[c(1,5,10)],digit=1)),col=clrssat[c(1,5,10)],lty=1)
plot(test$time_max,test$max_infect,col=alpha(color.gradient(test$pind,c(mygreen,myred)),.2),main="indiv learning",pch=20,xlab="Time Max",ylab="Max Infected")
legend("topright",legend=paste("P indiv learning=",unique(test$pind)),col=clrssat[c(1,5,10)],lty=1)
dev.off()


repetBest=sapply(1:100,function(i){print(i);abmSIR(poptest,1500,p=c(1,.2),di=2,i0=1,inf=.2,sat=1,xsize=xsize,ysize=ysize,visu=F,ap=F,ts=T,file=F,p_i=.5,log=F)$timeseries[,2]})
repetBestDos=parSapply(cl,1:100,function(i){print(i);abmSIR(poptest,1500,p=c(1,.2),di=2,i0=1,inf=.9,sat=1000,xsize=xsize,ysize=ysize,visu=F,ap=F,ts=T,file=F,p_i=.5,log=F)$timeseries[,2]})
neutralbad=sapply(1:10,function(i){print(i);abmSIR(poptest,1500,p=c(1,1),di=2,i0=1,inf=9,sat=1000,xsize=xsize,ysize=ysize,visu=F,ap=F,ts=T,file=F,p_i=1,log=F)$timeseries[,2]})
library(parallel)

save(file="worst.bin",worst)
save(file="best.bin",best)
load("best.bin")
cl <- makeForkCluster(3,outfile="")
neutralGood=parSapply(cl,1:100,function(i){print(i);abmSIR(poptest,1500,p=c(.1,.1),di=2,i0=1,inf=9,sat=1000,xsize=xsize,ysize=ysize,visu=F,ap=F,ts=T,file=F,p_i=1,log=F)$timeseries[,2]})
neutralBad=parSapply(cl,1:100,function(i){n=sample(nrow(best),1);print(paste(i,n));abmSIR(poptest,1500,p=c(1,1),di=2,i0=1,inf=best$inf[n],pind=best$pind[n],inf_r=best$inf_r[n],sat_r=best$sat_r[n],sat=best$sat[n],xsize=xsize,ysize=ysize,visu=F,ap=F,ts=T,file=F,p_i=1,log=F)$timeseries[,2]})
repetBest=sapply(1:100,function(i){print(i);abmSIR(poptest,1500,p=c(1,.2),di=2,i0=1,inf=.2,sat=1,xsize=xsize,ysize=ysize,visu=F,ap=F,ts=T,file=F,p_i=.5,log=F)$timeseries[,2]})
stopCluster(cl)

neutralbadList=lapply(1:nrow(neutralbad),function(i)neutralbad[i,])   
neutralGoodList=lapply(1:nrow(neutralGood),function(i)neutralGood[i,])   
repetBestDosList=lapply(1:nrow(repetBestDos),function(i)repetBestDos[i,])   
repetBestList=lapply(1:nrow(repetBest),function(i)repetBest[i,])   

subset=seq(10,1100,5)
acs=seq(1,length(subset),length.out=5)

png("5e8779c56517890001536101/figures/fullTraj.png",width=800,height=400,pointsize=17)
par(mfrow=c(1,2))
plot(1,1,type="n",ylim=c(0,500),xlim=c(1,1100),main="No learning")
legend("topright",legend=c("All bad","All good"),col=c(myred,mygreen),lwd=2)
apply(neutralGood,2,function(i)lines(i,col=alpha(mygreen,.2),lwd=2))
apply(neutralbad,2,function(i)lines(i,col=alpha(myred,.1),lwd=2))
plot(1,1,type="n",ylim=c(0,500),xlim=c(1,1100),main="Learning")
apply(repetBest,2,function(i)lines(i,col=alpha(mygreen,.2),lwd=2))
apply(repetBestDos,2,function(i)lines(i,col=alpha(myred,.1),lwd=2))
legend("topright",legend=c("Revert","No revert"),col=c(myred,mygreen),lwd=2)
dev.off()

png("5e8779c56517890001536101/figures/meanTrajIllu.png",width=800,height=800,pointsize=17)
plot(1,1,type="n",ylim=c(0,500),xlim=c(1,1500),ylab="number of infected people",xlab="time",xaxt="n",yaxt="n")
legend("topright",legend=c("No social distancing","Social distancing"),fill=c(myred,mygreen))
meanGood=apply(neutralGood,1,mean)
meanBad=apply(neutralbad,1,mean)
polygon(c(0:1500),meanGood,col=alpha(mygreen,.5),border=NA)
polygon(c(0:1500),meanBad,col=alpha(myred,.5),border=NA)
lines(meanGood,col=alpha(mygreen,.9),lwd=2)
lines(meanBad,col=alpha(myred,.9),lwd=2)
abline(h=max(meanBad),lty=2,col=myred)
abline(v=which.max(meanBad),lty=2,col=myred)
abline(h=max(meanGood),lty=2,col=mygreen)
abline(v=which.max(meanGood),lty=2,col=mygreen)
#text(paste("max infected=",round(max(meanBad))),x=1500,y=max(meanBad)+10,pos=2)
#mtext(paste("time max=",round(which.max(meanBad))),at=which.max(meanBad),adj=1)
#text(paste("max infected=",round(max(meanGood))),x=1500,y=max(meanGood)+10,pos=2)
#mtext(paste("time max=",round(which.max(meanGood))),at=which.max(meanGood),adj=0)
text(paste("max infected no SD"),x=1500,y=max(meanBad)+10,pos=2)
mtext(paste("time max no SD"),at=which.max(meanBad),adj=1)
text(paste("max infected SD"),x=1500,y=max(meanGood)+10,pos=2)
mtext(paste("time max SD"),at=which.max(meanGood),adj=0)
abline(v=which.max(meanBad),lty=2,col=myred)
abline(h=max(meanGood),lty=2,col=mygreen)
abline(v=which.max(meanGood),lty=2,col=mygreen)
dev.off()



png("5e8779c56517890001536101/figures/fullTrajHDR.png",width=800,height=400,pointsize=17)
par(mfrow=c(1,2))
hdr.boxplot(neutralGoodList[subset],pch=".",outline=F,,prob=c(50,95),space=0,ylab="#infected",border=NA,h=10,col=mygreen,ylim=c(0,500),main="No learning")
legend("topright",legend=c("All bad","All good"),fill=c(myred,mygreen))
par(new=T)
hdr.boxplot(neutralbadList[subset],pch=".",outline=F,,prob=c(50,95),space=0,ylab="#infected",border=NA,h=10,col=myred,ylim=c(0,500))
axis(1,at=acs,subset[acs])


hdr.boxplot(repetBestList[subset],pch=".",outline=F,,prob=c(50,95),space=0,ylab="#infected",border=NA,h=10,col=mygreen,ylim=c(0,500),main="Learning")
legend("topright",legend=c("Revert","No revert"),fill=c(myred,mygreen))
par(new=T)
hdr.boxplot(repetBestDosList[subset],pch=".",outline=F,,prob=c(50,95),space=0,ylab="#infected",border=NA,h=10,col=myred,ylim=c(0,500))
axis(1,at=acs,subset[acs])
dev.off()


##rsync -avz --info=progress2 --exclude "*simu_*.bin" volos:/home/share/simon/abmEpi/fullRandom .

for(name in c("randomSL","bestSL")){

aa=lapply(unlist(lapply(list.dirs(name),list.files,pattern="allresults.bin",full.names=T)),function(f){load(f);return(allresults)})
allresults=do.call("rbind",aa)
allresults=allresults[-which(allresults$max_infect <4),]
allresults$scores=(1-bscore(allresults$time_max) + bscore(allresults$max_infect))/2

nrow(allresults)
png(paste0("5e8779c56517890001536101/figures/alldim",name,".png"),width=1000,height=1200,pointsize=25)
par(mfrow=c(3,2),cex=.5)
test=allresults
plot(test$time_max,test$max_infect,col=alpha(color.gradient(log(test$sat),c(mygreen,"yellow", myred)),.2),main="steepness",pch=20,xlab="Time Max",ylab="Max Infected")
sset=sort(log10(test$sat))
cols=alpha(color.gradient(sset,c(mygreen,"yellow", myred)),1)
lsset=seq(1,length(sset),length.out=4)
legend("topright",legend=paste("steepness=10^",round(sset[lsset])),col=cols[lsset],pch=20)
plot(test$time_max,test$max_infect,col=alpha(color.gradient(log(test$sat_r),c(mygreen,"yellow", myred)),.2),main="revert steepness",pch=20,xlab="Time Max",ylab="Max Infected")
legend("topright",legend=paste("steepness=10^",round(sset[lsset])),col=cols[lsset],pch=20)
plot(test$time_max,test$max_infect,col=alpha(color.gradient(test$inf,c(mygreen,"yellow", myred)),.2),main="inflection point",pch=20,xlab="Time Max",ylab="Max Infected")
sset=sort(test$inf)
cols=alpha(color.gradient(sset,c(mygreen,"yellow", myred)),1)
lsset=seq(1,length(sset),length.out=4)
legend("topright",legend=paste("inp",round(sset[lsset],digit=2)),col=cols[lsset],pch=20)
plot(test$time_max,test$max_infect,col=alpha(color.gradient(test$inf_r,c(mygreen,"yellow", myred)),.2),main="revert inif. point",pch=20,xlab="Time Max",ylab="Max Infected")
legend("topright",legend=paste("inp",round(sset[lsset],digit=2)),col=cols[lsset],pch=20)
plot(test$time_max,test$max_infect,col=alpha(color.gradient(test$pind,c(mygreen,"yellow", myred)),.2),main="individual learning",pch=20,xlab="Time Max",ylab="Max Infected")
sset=sort(test$pind)
cols=alpha(color.gradient(sset,c(mygreen,"yellow", myred)),1)
lsset=seq(1,length(sset),length.out=4)
legend("topright",legend=paste("ind. learn.",round(sset[lsset],digit=2)),col=cols[lsset],pch=20)
plot(test$time_max,test$max_infect,col=alpha(color.gradient(test$scores,c(mygreen,"yellow", myred)),.2),main="score",pch=20,xlab="Time Max",ylab="Max Infected")
sset=sort(test$scores)
cols=alpha(color.gradient(sset,c(mygreen,"yellow", myred)),1)
lsset=seq(1,length(sset),length.out=4)
legend("topright",legend=paste("score",round(sset[lsset],digit=2)),col=cols[lsset],pch=20)
dev.off()



for(d in c(.8,.6,.4,.2)){
    png(paste0("5e8779c56517890001536101/figures/posterior2d_",d*10,"_",name,".png"),width=1100,height=333,pointsize=20)
    best=allresults[allresults$scores<d,]
    par(mfrow=c(1,4))
    hdr.boxplot.2d(best$inf_r,log10(best$sat_r),prob=seq(20,100,10),shadecols=mygreen,xlim=c(0,1),ylim=c(-1,3),xlab="revert inf. point",ylab="rever steepness")
    hdr.boxplot.2d(best$inf,log10(best$sat),prob=seq(20,100,10),shadecols="yellow",xlim=c(0,1),ylim=c(-1,3),xlab="inflexion point",ylab="steepness")
    hdr.boxplot.2d(best$inf,best$pind,prob=seq(20,100,10),shadecols=myred,xlim=c(0,1),ylim=c(0,1),xlab="inflexion point",ylab="individual learning")
    hdr.boxplot.2d(log10(best$sat),best$pind,prob=seq(20,100,10),shadecols=mygreen,xlim=c(-1,3),ylim=c(0,1),xlab="steepness",ylab="individual learning")
    dev.off()
}
}


everyall=lapply(c("randomSL","bestSL") ,function(name){
aa=lapply(unlist(lapply(list.dirs(name),list.files,pattern="allresults.bin",full.names=T)),function(f){load(f);return(allresults)})
allresults=do.call("rbind",aa)
allresults=allresults[-which(allresults$max_infect <4),]
allresults$scores=(1-bscore(allresults$time_max) + bscore(allresults$max_infect))/2
return(allresults)
}
)
names(everyall)=c("randomSL","bestSL") 
    bestAll=lapply(everyall,function(i)i[i$scores<.3,])

name ="randomSL"

aa=lapply(unlist(lapply(list.dirs(name),list.files,pattern="allresults.bin",full.names=T)),function(f){load(f);return(allresults)})
allresults=do.call("rbind",aa)
allresults=allresults[-which(allresults$max_infect <4),]
allresults$scores=(1-bscore(allresults$time_max) + bscore(allresults$max_infect))/2

nrow(allresults)
png(paste0("5e8779c56517890001536101/figures/alldim.png"),width=1000,height=1200,pointsize=25)
par(mfrow=c(3,2),cex=.5)
test=allresults
plot(test$time_max,test$max_infect,col=alpha(color.gradient(log(test$sat),c(mygreen,"yellow", myred)),.2),main="steepness",pch=20,xlab="Time Max",ylab="Max Infected")
sset=sort(log10(test$sat))
cols=alpha(color.gradient(sset,c(mygreen,"yellow", myred)),1)
lsset=seq(1,length(sset),length.out=4)
legend("topright",legend=paste("steepness=10^",round(sset[lsset])),col=cols[lsset],pch=20)
plot(test$time_max,test$max_infect,col=alpha(color.gradient(log(test$sat_r),c(mygreen,"yellow", myred)),.2),main="revert steepness",pch=20,xlab="Time Max",ylab="Max Infected")
legend("topright",legend=paste("steepness=10^",round(sset[lsset])),col=cols[lsset],pch=20)
plot(test$time_max,test$max_infect,col=alpha(color.gradient(test$inf,c(mygreen,"yellow", myred)),.2),main="inflection point",pch=20,xlab="Time Max",ylab="Max Infected")
sset=sort(test$inf)
cols=alpha(color.gradient(sset,c(mygreen,"yellow", myred)),1)
lsset=seq(1,length(sset),length.out=4)
legend("topright",legend=paste("inp",round(sset[lsset],digit=2)),col=cols[lsset],pch=20)
plot(test$time_max,test$max_infect,col=alpha(color.gradient(test$inf_r,c(mygreen,"yellow", myred)),.2),main="revert inif. point",pch=20,xlab="Time Max",ylab="Max Infected")
legend("topright",legend=paste("inp",round(sset[lsset],digit=2)),col=cols[lsset],pch=20)
plot(test$time_max,test$max_infect,col=alpha(color.gradient(test$pind,c(mygreen,"yellow", myred)),.2),main="individual learning",pch=20,xlab="Time Max",ylab="Max Infected")
sset=sort(test$pind)
cols=alpha(color.gradient(sset,c(mygreen,"yellow", myred)),1)
lsset=seq(1,length(sset),length.out=4)
legend("topright",legend=paste("ind. learn.",round(sset[lsset],digit=2)),col=cols[lsset],pch=20)
plot(test$time_max,test$max_infect,col=alpha(color.gradient(test$scores,c(mygreen,"yellow", myred)),.2),main="score",pch=20,xlab="Time Max",ylab="Max Infected")
sset=sort(test$scores)
cols=alpha(color.gradient(sset,c(mygreen,"yellow", myred)),1)
lsset=seq(1,length(sset),length.out=4)
lapply(seq(.2,.8),function(s){t=min(test$time_max):max(test$time_max);lines(t,m(t=t,Rm=max(test$max_infect)-min(test$max_infect),mm=min(test$max_infect),Rt=max(test$time_max)-min(test$time_max),mt=min(test$time_max),s=s),col="grey",lwd=2)})
legend("topright",legend=paste("score",round(sset[lsset],digit=2)),col=cols[lsset],pch=20)
dev.off()



for(d in c(1,.9,.8,.6,.4,.2)){
    png(paste0("5e8779c56517890001536101/figures/posterior2d_",d*10,".png"),width=1100,height=333,pointsize=20)
    best=allresults[allresults$scores<d,]
    par(mfrow=c(1,4))
    hdr.boxplot.2d(best$inf_r,log10(best$sat_r),prob=seq(20,100,10),shadecols="grey",xlim=c(0,1),ylim=c(-1,3),xlab="revert inf. point",ylab="rever steepness")
    hdr.boxplot.2d(best$inf,log10(best$sat),prob=seq(20,100,10),shadecols="grey",xlim=c(0,1),ylim=c(-1,3),xlab="inflexion point",ylab="steepness")
    hdr.boxplot.2d(best$inf,best$pind,prob=seq(20,100,10),shadecols="grey",xlim=c(0,1),ylim=c(0,1),xlab="inflexion point",ylab="individual learning")
    hdr.boxplot.2d(log10(best$sat),best$pind,prob=seq(20,100,10),shadecols="grey",xlim=c(-1,3),ylim=c(0,1),xlab="steepness",ylab="individual learning")
    dev.off()
}


