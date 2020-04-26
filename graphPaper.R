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


bdistance <- function(x)(x-min(x))/(max(x)-min(x))
bscore <- function(x)(x-min(x))/(max(x)-min(x))



##### Grpah papers
##get data
name="midCurveAllBad"
name="midCurveAllBadBestSLS"
name="midCurve15Good"
#name="testradius"

aa=lapply(unlist(lapply(list.dirs(name),list.files,pattern="allresults.bin",full.names=T)),function(f){load(f);return(allresults)})
allresults=do.call("rbind",aa)
allresults=allresults[-which(allresults$max_infect <4),]
allresults$distances=(1-bdistance(allresults$time_max) + bdistance(allresults$max_infect))/2
allresults$distances2=(1-bdistance(allresults$time_max150) + bdistance(allresults$max_infect150))/2
allresults$distances3=(1-bdistance(allresults$time_max250) + bdistance(allresults$max_infect250))/2

#figure 1
load("cydia/neutralBad.bin")
load("cydia/neutralGood.bin")
load("cydia/neutralGoodG.bin")
load("cydia/repetBest.bin")
load("cydia/repetWorst.bin")

png("5e8779c56517890001536101/figures/meanTrajIllu.png",pointsize=17, res = 100, width = 7, height = 7, units = "in")
par(mar=c(2,2,2,2))
plot(1,1,type="n",ylim=c(0,500),xlim=c(1,1500),ylab="",xlab="",xaxt="n",yaxt="n")
legend("topright",legend=c("No social distancing","Social distancing"),fill=c(myred,mygreen),cex=.95)
mtext("time",1,1)
mtext("number of infected people",2,1)
meanBad=apply(neutralGood,1,mean)
meanGood=apply(neutralGoodG,1,mean)
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
text(paste("max infected no SD"),x=1500,y=max(meanBad)+10,pos=2,cex=.8)
mtext(paste("time max no SD"),at=which.max(meanBad),adj=1,cex=.8)
text(paste("max infected SD"),x=1500,y=max(meanGood)+10,pos=2,cex=.8)
mtext(paste("time max SD"),at=which.max(meanGood),adj=0,cex=.8)
abline(v=which.max(meanBad),lty=2,col=myred)
abline(h=max(meanGood),lty=2,col=mygreen)
abline(v=which.max(meanGood),lty=2,col=mygreen)
dev.off()

### Figure3
png("5e8779c56517890001536101/figures/distribSimulation.png",pointsize=17, res = 100, width = 7, height = 7, units = "in")
plot(allresults$time_max,allresults$max_infect,col=alpha(color.gradient(allresults$distances,c(mygreen,"yellow", myred)),.2),main="distance",pch=20,xlab="Time Max",ylab="Max Infected")
sset=sort(allresults$distances)
cols=alpha(color.gradient(sset,c(mygreen,"yellow", myred)),1)
lsset=seq(1,length(sset),length.out=4)
legend("topright",legend=sapply(round(sset[lsset],digit=2),function(d)as.expression(bquote(delta==.(d)))),col=cols[lsset],pch=20)
dev.off()

### Figure4a
png("5e8779c56517890001536101/figures/distribSimulationGreyZonedS30.png",pointsize=17, res = 100, width = 7, height = 7, units = "in")
plot(allresults$time_max,allresults$max_infect,col=alpha(color.gradient(allresults$distances,c(mygreen,"yellow", myred)),.2),main="distance",pch=20,xlab="Time Max",ylab="Max Infected")
sset=sort(allresults$distances)
cols=alpha(color.gradient(sset,c(mygreen,"yellow", myred)),1)
lsset=seq(1,length(sset),length.out=4)
v=m(t=t,Rm=max(allresults$max_infect)-min(allresults$max_infect),mm=min(allresults$max_infect),Rt=max(allresults$time_max)-min(allresults$time_max),mt=min(allresults$time_max),s=.3)
t=400:1500
polygon(c(100,t,1510),c(0,v,0),col=alpha("dark blue",.5),lwd=2,angle=60,density=5)
legend("topright",legend=c(sapply(round(sset[lsset],digit=2),function(d)as.expression(bquote(delta==.(d)))),expression(delta<.25)),col=c(cols[lsset],NA),fill=c(rep(0,4),alpha("dark blue",.5)),border=c(rep(NA,4),alpha("dark blue",.5)),angle=c(rep(NA,4),60),density=c(rep(NA,4),25),pch=20)
#points(best1000$time_max,best1000$max_infect,col="green",main="distance",pch=20,xlab="Time Max",ylab="Max Infected")
#points(best500$time_max,best500$max_infect,col="red",main="distance",pch=20,xlab="Time Max",ylab="Max Infected")
#points(best250$time_max,best250$max_infect,col="dark blue",main="distance",pch=20,xlab="Time Max",ylab="Max Infected")
dev.off()


### Figure4b
png("5e8779c56517890001536101/figures/distribSimulationGreyZonedS60.png",pointsize=17, res = 100, width = 7, height = 7, units = "in")
plot(allresults$time_max,allresults$max_infect,col=alpha(color.gradient(allresults$distances,c(mygreen,"yellow", myred)),.2),main="distance",pch=20,xlab="Time Max",ylab="Max Infected")
sset=sort(allresults$distances)
cols=alpha(color.gradient(sset,c(mygreen,"yellow", myred)),1)
lsset=seq(1,length(sset),length.out=4)
v=m(t=t,Rm=max(allresults$max_infect)-min(allresults$max_infect),mm=min(allresults$max_infect),Rt=max(allresults$time_max)-min(allresults$time_max),mt=min(allresults$time_max),s=.6)
polygon(c(100,t,1010),c(0,v,0),col=alpha("dark blue",.5),lwd=2,angle=60,density=5)
legend("topright",legend=c(sapply(round(sset[lsset],digit=2),function(d)as.expression(bquote(delta==.(d)))),expression(delta<.6)),col=c(cols[lsset],NA),fill=c(rep(0,4),alpha("dark blue",.5)),border=c(rep(NA,4),alpha("dark blue",.5)),angle=c(rep(NA,4),60),density=c(rep(NA,4),25),pch=20)
dev.off()

dev.new()
####
##take 500 best simulation
best=allresults[order(allresults$distances),][1:500,]
best250=allresults[order(allresults$distances),][1:250,]
best1000=allresults[order(allresults$distances),][1:10000,]
worst=allresults[order(allresults$distances,decreasing=T),][1:500,]
save(file="worst.bin",worst)
save(file="best.bin",best)

### FIGURE 5
png(paste0("5e8779c56517890001536101/figures/switchfunctions.png"),pointsize=17, res = 100, width = 13, height = 7, units = "in")
    par(mfrow=c(1,2))
    x=seq(0,1,.01)
    plot(x,sig(x),type="n",ylim=c(0,1),xlim=c(0,1),xlab="proportion of infected people of the same age", ylab="probability to switch behavior",main="P(NA->A)")
    for(i in 1:nrow(best)) lines(x,sig(x,b=best$inf[i],a=best$sat[i]),ylim=c(0,1),xlim=c(0,1),col=alpha("dark blue",.1),lwd=2) 
    #lines(x,sig(x,b=mean(best$inf),a=mean(best$sat)),ylim=c(0,1),xlim=c(0,1),col=alpha("black",.6),lwd=2) 
    plot(x,sig(x),type="n",ylim=c(0,1),xlim=c(0,1),xlab="proportion of infected people of the same age", ylab="probability to switch behavior",main="P(A->NA)")
    for(i in 1:nrow(best)) lines(x,sig(x,b=best$inf_r[i],a=best$sat_r[i]),ylim=c(0,1),xlim=c(0,1),col=alpha("dark blue",.1),lwd=2) 
    #lines(x,sig(x,b=mean(best$inf_r),a=mean(best$sat_r)),ylim=c(0,1),xlim=c(0,1),col=alpha("black",.6),lwd=2) 
dev.off()




#### FIGURE6

png(paste0("5e8779c56517890001536101/figures/meanPosteriors.png"),pointsize=17, res = 100, width = 7, height = 7, units = "in")
par(mfrow=c(1,1))
plot(1,1,type="n",ylim=c(0,500),xlim=c(1,1500),main="",ylab="#infected",xlab="time")
legend("topright",legend=c("No learning (all NA)","No learning (all A)","Worst simulations","Best Simulations"),col=c(myred,mygreen),lwd=2,lty=c(1,1,2,2))
lines(apply(neutralGoodG,1,mean),col=alpha(mygreen,.9),lwd=2)
lines(apply(neutralGood,1,mean),col=alpha(myred,.9),lwd=2)
lines(apply(repetBest,1,mean),col=alpha(mygreen,.9),lwd=2,lty=2)
lines(apply(repetWorst,1,mean),col=alpha(myred,.9),lwd=2,lty=2)
dev.off()


##FIGURE 7

png("5e8779c56517890001536101/figures/posterior_pind.png",pointsize=17, res = 100, width = 7, height = 7, units = "in")
 plot2dens(prior=allresults$pind,A=best$pind,cols=c(P="white",A=alpha("dark blue",.5),NA),from=0,to=1,main="posterior proba individual learning",xlab="proba individual learning")
dev.off()

png("5e8779c56517890001536101/figures/posterior_sl_rad.png",pointsize=17, res = 100, width = 7, height = 7, units = "in")
 plot2dens(prior=allresults$sl_rad,A=best$sl_rad,cols=c(P="white",A=alpha("dark blue",.5),NA),from=1,to=100,main="posterior proba radius social learning",xlab="radius social learning")
dev.off()

png("5e8779c56517890001536101/figures/posterior_sat_r.png",pointsize=17, res = 100, width = 7, height = 7, units = "in")
 plot2dens(prior=log10(allresults$sat_r),A=log10(best$sat_r),cols=c(P="white",A=alpha("dark blue",.5),NA),from=-1,to=3,main="posterior revert steepness",xlab="revert steepness (log10)",xaxt="n")
dev.off()

png("5e8779c56517890001536101/figures/posterior_sat.png",pointsize=17, res = 100, width = 7, height = 7, units = "in")
 plot2dens(prior=log10(allresults$sat),A=log10(best$sat),cols=c(P="white",A=alpha("dark blue",.5),NA),from=-1,to=3,xlab="steepness (log10)",main="posterior steepness",xaxt="n")
dev.off()

png("5e8779c56517890001536101/figures/posterior_inf_r.png",pointsize=17, res = 100, width = 7, height = 7, units = "in")
 plot2dens(prior=(allresults$inf_r),A=(best$inf_r),cols=c(P="white",A=alpha("dark blue",.5),NA),from=0,to=1,main="posterior revert inflection point",xlab="revert inflection point",xaxt="n")
dev.off()

png("5e8779c56517890001536101/figures/posterior_inf.png",pointsize=17, res = 100, width = 7, height = 7, units = "in")
 plot2dens(prior=(allresults$inf),A=(best$inf),cols=c(P="white",A=alpha("dark blue",.5),NA),from=0,to=1,main="posterior inflection point",xlab="inflection point",xaxt="n")
dev.off()



#### FIGURE 8
for(d in c(1,.9,.8,.6,.4,.3,.2)){
    png(paste0("5e8779c56517890001536101/figures/posterior2d_",d*10,".png"),pointsize=30, res = 100, width = 25, height = 6, units = "in")
    best=allresults[allresults$distances<d,]
    par(mfrow=c(1,6))
    par(mar=c(4,4,1,1))
    hdr.boxplot.2d(best$inf_r,log10(best$sat_r),prob=seq(20,100,10),shadecols=alpha("dark blue",.9),xlim=c(0,1),ylim=c(-1,3),xlab="revert inf. point",ylab="revert steepness (log10)")
    hdr.boxplot.2d(best$inf,log10(best$sat),prob=seq(20,100,10),shadecols=alpha("dark blue",.9),xlim=c(0,1),ylim=c(-1,3),xlab="inflexion point",ylab="steepness (log10)")
    hdr.boxplot.2d(best$inf,best$pind,prob=seq(20,100,10),shadecols=alpha("dark blue",.9),xlim=c(0,1),ylim=c(0,1),xlab="inflexion point",ylab="individual learning")
    hdr.boxplot.2d(log10(best$sat),best$pind,prob=seq(20,100,10),shadecols=alpha("dark blue",.9),xlim=c(-1,3),ylim=c(0,1),xlab="steepness (log10)",ylab="individual learning")
    hdr.boxplot.2d(best$sl_rad,best$pind,prob=seq(20,100,10),shadecols=alpha("dark blue",.9),xlim=c(1,100),ylim=c(0,1),xlab="radius social learning",ylab="individual learning")
    hdr.boxplot.2d(best$sl_rad,best$inf,prob=seq(20,100,10),shadecols=alpha("dark blue",.9),xlim=c(1,100),ylim=c(0,1),xlab="radius social learning",ylab="infection points")
    dev.off()
}




######

load("allresults.bin")
allresults=allresults[-which(allresults$max_infect <4),]
allresults$distances=(1-bdistance(allresults$time_max) + bdistance(allresults$max_infect))/2

png("5e8779c56517890001536101/figures/exploringIndividualLearning.png",width=1000,height=1000,pointsize=25)
par(mfrow=c(3,3))
for(p_i in unique(allresults$pind)){
    test=allresults[allresults$pind == p_i ,]
    a=lapply(unique(test$sat),function(s){ sapply(unique(test$inf),function(i)quantile(test$distances[test$inf ==i & test$sat== s])[2:4])})
    b=lapply(unique(test$sat),function(s){ sapply(unique(test$inf),function(i)quantile(test$time_max[test$inf ==i & test$sat== s])[2:4])})
    c=lapply(unique(test$sat),function(s){ sapply(unique(test$inf),function(i)quantile(test$max_inf[test$inf ==i & test$sat== s])[2:4])})
    plotQuartiles(a,xlab="inp",ylab="distance",main=paste("proba individual learning=",p_i))
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

load("best.bin")
cl <- makeForkCluster(3,outfile="")
neutralGood=parSapply(cl,1:100,function(i){print(i);abmSIR(poptest,1500,p=c(.1,.1),di=2,i0=1,inf=9,sat=1000,xsize=xsize,ysize=ysize,visu=F,ap=F,ts=T,file=F,p_i=1,log=F)$timeseries[,2]})
neutralBad=parSapply(cl,1:100,function(i){n=sample(nrow(best),1);print(paste(i,n));abmSIR(poptest,1500,p=c(1,1),di=2,i0=1,inf=best$inf[n],pind=best$pind[n],inf_r=best$inf_r[n],sat_r=best$sat_r[n],sat=best$sat[n],xsize=xsize,ysize=ysize,visu=F,ap=F,ts=T,file=F,p_i=1,log=F)$timeseries[,2]})
repetBest=sapply(1:100,function(i){print(i);abmSIR(poptest,1500,p=c(1,.2),di=2,i0=1,inf=.2,sat=1,xsize=xsize,ysize=ysize,visu=F,ap=F,ts=T,file=F,p_i=.5,log=F)$timeseries[,2]})
stopCluster(cl)

load("cydia/neutralBad.bin")
load("cydia/neutralGood.bin")
load("cydia/neutralGoodG.bin")
load("cydia/repetBest.bin")
load("cydia/repetWorst.bin")
neutralbadList=lapply(1:nrow(neutralBad),function(i)neutralBad[i,])   
neutralGoodList=lapply(1:nrow(neutralGood),function(i)neutralGood[i,])   
repetBestList=lapply(1:nrow(repetBest),function(i)repetBest[i,])   
repetWorstList=lapply(1:nrow(repetWorst),function(i)repetWorst[i,])   

subset=seq(10,1100,5)
acs=seq(1,length(subset),length.out=5)

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


load("cydia/repetBest.bin") 
load("cydia/repetWorst.bin") 
load("cydia/neutralBad.bin") 
load("cydia/neutralGood.bin") 



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

allexpe=sapply(c("randomSL","bestSL","testradius"),function(name){

aa=lapply(unlist(lapply(list.dirs(name),list.files,pattern="allresults.bin",full.names=T)),function(f){load(f);return(allresults)})
allresults=do.call("rbind",aa)
allresults=allresults[-which(allresults$max_infect <4),]
allresults$distances=(1-bdistance(allresults$time_max) + bdistance(allresults$max_infect))/2
nrow(allresults)
return(allresults)
})

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
plot(test$time_max,test$max_infect,col=alpha(color.gradient(test$distances,c(mygreen,"yellow", myred)),.2),main="distance",pch=20,xlab="Time Max",ylab="Max Infected")
sset=sort(test$distances)
cols=alpha(color.gradient(sset,c(mygreen,"yellow", myred)),1)
lsset=seq(1,length(sset),length.out=4)
legend("topright",legend=paste("distance",round(sset[lsset],digit=2)),col=cols[lsset],pch=20)
dev.off()



for(d in c(.8,.6,.4,.2)){
    png(paste0("5e8779c56517890001536101/figures/posterior2d_",d*10,"_",name,".png"),width=1100,height=333,pointsize=20)
    best=allresults[allresults$distances<d,]
    par(mfrow=c(1,6))
    hdr.boxplot.2d(best$inf_r,log10(best$sat_r),prob=seq(20,100,10),shadecols=mygreen,xlim=c(0,1),ylim=c(-1,3),xlab="revert inf. point",ylab="rever steepness")
    hdr.boxplot.2d(best$inf,log10(best$sat),prob=seq(20,100,10),shadecols="yellow",xlim=c(0,1),ylim=c(-1,3),xlab="inflexion point",ylab="steepness")
    hdr.boxplot.2d(best$inf,best$pind,prob=seq(20,100,10),shadecols=myred,xlim=c(0,1),ylim=c(0,1),xlab="inflexion point",ylab="individual learning")
    hdr.boxplot.2d(log10(best$sat_r),best$pind,prob=seq(20,100,10),shadecols=mygreen,xlim=c(-1,3),ylim=c(0,1),xlab="steepness",ylab="individual learning")
    hdr.boxplot.2d(best$sl_rad,best$pind,prob=seq(20,100,10),shadecols=myred,xlim=c(1,50),ylim=c(0,1),xlab="radius social learning",ylab="individual learning")
    hdr.boxplot.2d(best$sl_rad,log10(best$sat),prob=seq(20,100,10),shadecols=myred,xlim=c(0,100),ylim=c(-1,3),xlab="radius social learning",ylab="steepness")
    dev.off()
}
}


everyall=lapply(c("randomSL","bestSL") ,function(name){
aa=lapply(unlist(lapply(list.dirs(name),list.files,pattern="allresults.bin",full.names=T)),function(f){load(f);return(allresults)})
allresults=do.call("rbind",aa)
allresults=allresults[-which(allresults$max_infect <4),]
allresults$distances=(1-bdistance(allresults$time_max) + bdistance(allresults$max_infect))/2
return(allresults)
}
)
names(everyall)=c("randomSL","bestSL") 
    bestAll=lapply(everyall,function(i)i[i$distances<.3,])

name ="randomSL"

aa=lapply(unlist(lapply(list.dirs(name),list.files,pattern="allresults.bin",full.names=T)),function(f){load(f);return(allresults)})
allresults=do.call("rbind",aa)
allresults=allresults[-which(allresults$max_infect <4),]
allresults$distances=(1-bscore(allresults$time_max) + bscore(allresults$max_infect))/2
allresults$distances2=(1-bscore(allresults$time_max250) + bscore(allresults$max_infect250))/2

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
plot(test$time_max,test$max_infect,col=alpha(color.gradient(test$distances,c(mygreen,"yellow", myred)),.2),main="distance",pch=20,xlab="Time Max",ylab="Max Infected")
sset=sort(test$distances)
cols=alpha(color.gradient(sset,c(mygreen,"yellow", myred)),1)
lsset=seq(1,length(sset),length.out=4)
lapply(c(.3,.6),function(s){t=200:1000;lines(t,m(t=t,Rm=max(test$max_infect)-min(test$max_infect),mm=min(test$max_infect),Rt=max(test$time_max)-min(test$time_max),mt=min(test$time_max),s=s),col="grey",lwd=2)})
legend("topright",legend=paste("distance",round(sset[lsset],digit=2)),col=cols[lsset],pch=20)
dev.off()

layout(matrix(c(1,1,1,1,2,3,4,5),nrow=4,ncol=2),width=c(2,1))  



    png(paste0("5e8779c56517890001536101/figures/posterior2d_inf_r_sat_r.png"),pointsize=17, res = 100, width = 7, height = 7, units = "in")
    hdr.boxplot.2d(best$inf_r,log10(best$sat_r),prob=seq(20,100,10),shadecols=alpha("dark blue",.9),xlim=c(0,1),ylim=c(-1,3),xlab="revert inf. point",ylab="revert steepness")
    dev.off()
    png(paste0("5e8779c56517890001536101/figures/posterior2d_inf_sat_r.png"),pointsize=17, res = 100, width = 7, height = 7, units = "in")
    hdr.boxplot.2d(best$inf,log10(best$sat_r),prob=seq(20,100,10),shadecols=alpha("dark blue",.9),xlim=c(0,1),ylim=c(-1,3),xlab="inf. point",ylab="revert steepness")
dev.off()


for(d in c(1,.9,.8,.6,.4,.3)){
    png(paste0("5e8779c56517890001536101/figures/posterior2d_",d*10,".png"),width=1600,height=333,pointsize=20)
    best=allresults[allresults$distances<d,]
    par(mfrow=c(1,5))
    hdr.boxplot.2d(best$inf_r,log10(best$sat_r),prob=seq(20,100,10),shadecols="grey",xlim=c(0,1),ylim=c(-1,3),xlab="revert inf. point",ylab="rever steepness")
    hdr.boxplot.2d(best$inf,log10(best$sat),prob=seq(20,100,10),shadecols="grey",xlim=c(0,1),ylim=c(-1,3),xlab="inflexion point",ylab="steepness")
    hdr.boxplot.2d(best$inf,best$pind,prob=seq(20,100,10),shadecols="grey",xlim=c(0,1),ylim=c(0,1),xlab="inflexion point",ylab="individual learning")
    hdr.boxplot.2d(log10(best$sat),best$pind,prob=seq(20,100,10),shadecols="grey",xlim=c(-1,3),ylim=c(0,1),xlab="steepness",ylab="individual learning")
    hdr.boxplot.2d(best$sl_rad,best$pind,prob=seq(20,100,10),shadecols="grey",xlim=c(1,100),ylim=c(0,1),xlab="inflexion point",ylab="individual learning")
    dev.off()
}

