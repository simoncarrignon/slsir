m<-function(t,Rm,Rt,mm,mt,s)(Rm*(2*s-1+(t-mt)/Rt))+mm 

source("abmEpi.R")
library(RColorBrewer)

color_class=rev(brewer.pal(5,"PuBu"))
colorbest=color_class[1]
myred=rgb(r=213,g=94,b=0,maxColorValue=255)
myred="red"
mygreen=rgb(r=0,g=158,b=115,maxColorValue=255)
mygreen=colorbest



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
name="testradius"
name="burnin100BestSLSFixed"

aa=lapply(unlist(lapply(list.dirs(name),list.files,pattern="allresults.bin",full.names=T)),function(f){load(f);return(allresults)})
allresults=do.call("rbind",aa)
#allresults=allresults[allresults$inf <.02 & allresults$sat >10,]
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

png("5e8779c56517890001536101/figures/meanTrajIllu.png",pointsize=17, res = 300, width = 7, height = 7, units = "in")
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
text(bquote(paste("max infected (",I[max],") no SD")),x=1500,y=max(meanBad)+10,pos=2,cex=.8)
mtext(bquote(paste("time max (",tau,") no SD")),at=which.max(meanBad),adj=1,cex=.8)
text(bquote(paste("max infected (",I[max],") SD")),x=1500,y=max(meanGood)+10,pos=2,cex=.8)
mtext(bquote(paste("time max (",tau,") SD")),at=which.max(meanGood),adj=0,cex=.8)
abline(v=which.max(meanBad),lty=2,col=myred)
abline(h=max(meanGood),lty=2,col=mygreen)
abline(v=which.max(meanGood),lty=2,col=mygreen)
dev.off()


##############FIGURE 2

    png("5e8779c56517890001536101/figures/illuSig.png",pointsize=25, res = 300, width = 15, height = 7, units = "in")
    par(mfrow=c(1,2))
    par(mar=c(4,4,1,1))
    x=seq(0,1,.01)   
    inp=seq(0,1,.05)   
    clrs=colorRampPalette(c("purple4","yellow"))(length(inp))
    #plot(x,sig(x),type="n",ylim=c(0,1),xlim=c(0,1),ylab=expression(sig(x,st=10,inp)),main=expression(inp %in%  "(" * list(0,1) * ")") )
    plot(x,sig(x),type="n",ylim=c(0,1),xlim=c(0,1),xlab="% infected",ylab=expression(P[switch]),main="" )
    for(i in 1:length(inp)) lines(x,sig(x,b=inp[i],a=4),ylim=c(0,1),xlim=c(0,1),col=clrs[i],lwd=2) 
    leg=seq(1,length(inp),length.out=5)
    legend("bottomright",legend=rev(sapply(inp[leg],function(d)as.expression(bquote(.(d))))),title=expression(nu),col=rev(clrs[leg]),lwd=2,cex=.8)

    tstp=25
    stp=rev(seq(-3,3,length.out=tstp))
    clrs=colorRampPalette(c("purple4","yellow"))(tstp)
    #plot(x,sig(x),type="n",ylim=c(0,1),xlim=c(0,1),ylab="P(B->G)~sig(x,st,inp=.5)",main=bquote(stp %in% "(" * list(10^.(stp[1]),10^.(stp[tstp])) * ")") )
    plot(x,sig(x),type="n",ylim=c(0,1),xlim=c(0,1),xlab="% infected",ylab=expression(P[switch]),main="" )
    for(i in rev(1:length(stp))) lines(x,sig(x,a=10^stp[i]),ylim=c(0,1),xlim=c(0,1),col=clrs[i],lwd=2) 
    leg=seq(1,tstp,length.out=5)
    #legend("bottomright",legend=sapply(round(stp[leg]),function(d)as.expression(bquote(kappa==.(10^d)))),col=clrs[leg],lwd=2)
    legend("bottomright",legend=sapply(round(stp[leg]),function(d)as.expression(bquote(.(10^d)))),col=clrs[leg],lwd=2,title=expression(kappa),cex=.8)
    dev.off()

### Figure3
png("5e8779c56517890001536101/figures/distribSimulation.png",pointsize=25, res = 300, width = 11, height = 11, units = "in")
rank=rank(allresults$distances)
orderscol=rank
color_class=rev(brewer.pal(5,"PuBu"))
colorbest=color_class[1]
orderscol=rep(alpha(color_class[5],.4),length(orderscol))
nvl=c(1000,2500,12500,62500)
for(i in (length(nvl):1)) orderscol[rank<=nvl[i]]=alpha(color_class[i],.4)
#plot(allresults$time_max,allresults$max_infect,bg=orderscol,main=expression(I[max]*" wrt "*tau),pch=21,xlab=expression(tau),ylab=expression(I[max]),lwd=.1)
plot(allresults$time_max,allresults$max_infect,bg=orderscol,main="",pch=21,xlab=expression(tau),ylab=expression(I[max]),lwd=.1)
legend("topright",legend=c(paste("<",nvl),"all"),pt.bg=color_class,pch=21,col=alpha(1,.5),pt.lwd=.1,title="Rank")
dev.off()

rank=rank(allresults$max_infect150)
orderscol=rank
color_class=rev(brewer.pal(5,"PuBu"))
colorbest=color_class[1]
orderscol=rep(alpha(color_class[5],.6),length(orderscol))
#nvl=c(100,2000,10000,62500)
for(i in (length(nvl):1)) orderscol[rank<=nvl[i]]=alpha(color_class[i],.6)
plot(allresults$time_max150,allresults$max_infect150,bg=orderscol,main="Maximum number of infected individual and time to reach this maximum for all simulations",pch=21,xlab=expression(tau),ylab=expression(I[max]),lwd=.1)
legend("topright",legend=c(paste("<",nvl),"all"),pt.bg=color_class,pch=21,pt.lwd=.1,title="Rank")

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

dev.new()
####
##take 500 best simulation
best=allresults[order(allresults$distances),][1:1000,]
best250=allresults[order(allresults$distances),][1:250,]
best1000=allresults[order(allresults$distances),][1:10000,]
worst=allresults[order(allresults$distances,decreasing=T),][1:1000,]
save(file="worst.bin",worst)
save(file="best.bin",best)
saveRDS(file="worst.bin",posterior$worst)
saveRDS(file="best.bin",posterior$best)
n=1000
posterior=list(
               best=allresults[order(allresults$distances),][1:n,],
               worst=allresults[order(allresults$distances,decreasing=T),][1:n,],
               time_max=allresults[order(allresults$time_max,decreasing=T),][1:n,],
               time_max150=allresults[order(allresults$time_max150,decreasing=T),][1:n,],
               time_max250=allresults[order(allresults$time_max250,decreasing=T),][1:n,],
               wtime_max150=allresults[order(allresults$time_max150,decreasing=F),][1:n,],
               wtime_max250=allresults[order(allresults$time_max250,decreasing=F),][1:n,],
               max_infect=allresults[order(allresults$max_infect,decreasing=F),][1:n,],
               max_infect150=allresults[order(allresults$max_infect150,decreasing=F),][1:n,],
               max_infect250=allresults[order(allresults$max_infect250,decreasing=F),][1:n,],
               wmax_infect150=allresults[order(allresults$max_infect150,decreasing=T),][1:n,],
               wmax_infect250=allresults[order(allresults$max_infect250,decreasing=T),][1:n,]
               )

### FIGURE 5
png(paste0("5e8779c56517890001536101/figures/switchfunctions.png"),pointsize=17, res = 100, width = 13, height = 7, units = "in")
    par(mfrow=c(1,2))
    x=seq(0,1,.01)
    plot(x,sig(x),type="n",ylim=c(0,1),xlim=c(0,1),xlab="proportion of infected people of the same age", ylab="probability to switch behavior",main="P(NA->A)")
    for(i in 1:nrow(posterior$best)) lines(x,sig(x,b=posterior$best$inf[i],a=posterior$best$sat[i]),ylim=c(0,1),xlim=c(0,1),col=alpha(colorbest,.1),lwd=2) 
    #lines(x,sig(x,b=mean(posterior$best$inf),a=mean(posterior$best$sat)),ylim=c(0,1),xlim=c(0,1),col=alpha("black",.6),lwd=2) 
    plot(x,sig(x),type="n",ylim=c(0,1),xlim=c(0,1),xlab="proportion of infected people of the same age", ylab="probability to switch behavior",main="P(A->NA)")
    for(i in 1:nrow(posterior$best)) lines(x,sig(x,b=posterior$best$inf_r[i],a=posterior$best$sat_r[i]),ylim=c(0,1),xlim=c(0,1),col=alpha(colorbest,.1),lwd=2) 
    #lines(x,sig(x,b=mean(best$inf_r),a=mean(best$sat_r)),ylim=c(0,1),xlim=c(0,1),col=alpha("black",.6),lwd=2) 
dev.off()
 
##FIGURE 7
mar=c(4,3,2,1)

colors=c(P="NA",A=alpha(colorbest,1),B=alpha(colorbest,1))
png("5e8779c56517890001536101/figures/posterior_pind.png",pointsize=30, res = 100, width = 10, height = 10, units = "in")
par(mar=mar)
 plot2dens(prior=1-allresults$pind,A=1-posterior$best$pind,B=1-posterior$max_infect150$pind,cols=colors,from=0,to=1,main="posterior proba social learning",xlab="proba social learning")
dev.off()

png("5e8779c56517890001536101/figures/posterior_sl_rad.png",pointsize=30, res = 100, width = 10, height = 10, units = "in")
par(mar=mar)
 plot2dens(prior=allresults$sl_rad,A=posterior$best$sl_rad,B=posterior$max_infect150$sl_rad,cols=colors,from=1,to=100,main="",xlab="radius social learning")
dev.off()

png("5e8779c56517890001536101/figures/posterior_sat_r.png",pointsize=30, res = 100, width = 10, height = 10, units = "in")
par(mar=mar)
 plot2dens(prior=log10(allresults$sat_r),A=log10(posterior$best$sat_r),B=log10(posterior$time_max150$sat_r),cols=colors,from=-1,to=3,main="",xlab=expression(kappa[r]),log=T)
dev.off()

png("5e8779c56517890001536101/figures/posterior_sat.png",pointsize=30, res = 100, width = 10, height = 10, units = "in")
par(mar=mar)
 plot2dens(prior=log10(allresults$sat),A=log10(posterior$best$sat),B=log10(posterior$max_infect150$sat),cols=colors,from=-1,to=3,xlab=expression(kappa),main="",xaxt="n",log=T)
dev.off()

png("5e8779c56517890001536101/figures/posterior_inf_r.png",pointsize=30, res = 100, width = 10, height = 10, units = "in")
par(mar=mar)
 plot2dens(prior=(allresults$inf_r),A=(posterior$best$inf_r),B=(posterior$max_infect150$inf_r),cols=colors,from=0,to=1,main="posterior revert inflection point",xlab=expression(nu[r]),xaxt="n")
legend("topleft",legend=c("prior","final time","150 timesteps"),fill=c(0,colorbest,colorbest),density=c(NA,NA,20))
dev.off()

png("5e8779c56517890001536101/figures/posterior_inf.png",pointsize=30, res = 100, width = 10, height = 10, units = "in")
par(mar=mar)
 plot2dens(prior=(allresults$inf),A=(posterior$best$inf),B=(posterior$max_infect150$inf),cols=colors,from=0,to=1,main="posterior inflection point",xlab=expression(nu),xaxt="n")
dev.off()


##FIGURE 7BIS
mar=c(4,3,2,1)

colors=c(P="NA",A=alpha(colorbest,1),B=alpha(colorbest,1))
png("5e8779c56517890001536101/figures/worst_pind.png",pointsize=30, res = 100, width = 10, height = 10, units = "in")
par(mar=mar)
 plot2dens(prior=1-posterior$worst$pind,A=1-posterior$best$pind,B=1-posterior$max_infect150$pind,cols=colors,from=0,to=1,main="",xlab="proba social learning")
dev.off()

png("5e8779c56517890001536101/figures/worst_sl_rad.png",pointsize=30, res = 100, width = 10, height = 10, units = "in")
par(mar=mar)
 plot2dens(prior=posterior$worst$sl_rad,A=posterior$best$sl_rad,B=posterior$max_infect150$sl_rad,cols=colors,from=1,to=100,main="posterior proba radius social learning",xlab="radius social learning")
dev.off()

png("5e8779c56517890001536101/figures/worst_sat_r.png",pointsize=30, res = 100, width = 10, height = 10, units = "in")
par(mar=mar)
 plot2dens(prior=log10(posterior$worst$sat_r),A=log10(posterior$best$sat_r),B=log10(posterior$time_max150$sat_r),cols=colors,from=-1,to=3,main="posterior revert steepness",xlab=expression(kappa[r]),log=T)
dev.off()

png("5e8779c56517890001536101/figures/worst_sat.png",pointsize=30, res = 100, width = 10, height = 10, units = "in")
par(mar=mar)
 plot2dens(prior=log10(posterior$worst$sat),A=log10(posterior$best$sat),B=log10(posterior$max_infect150$sat),cols=colors,from=-1,to=3,xlab=expression(kappa),main="posterior steepness",xaxt="n",log=T)
dev.off()

png("5e8779c56517890001536101/figures/worst_inf_r.png",pointsize=30, res = 100, width = 10, height = 10, units = "in")
par(mar=mar)
 plot2dens(prior=(posterior$worst$inf_r),A=(posterior$best$inf_r),B=(posterior$max_infect150$inf_r),cols=colors,from=0,to=1,main="",xlab=expression(nu[r]),xaxt="n")
legend("topleft",legend=c("prior","final time","150 timesteps"),fill=c(0,colorbest,colorbest),density=c(NA,NA,20))
dev.off()

png("5e8779c56517890001536101/figures/worst_inf.png",pointsize=30, res = 100, width = 10, height = 10, units = "in")
par(mar=mar)
 plot2dens(prior=(posterior$worst$inf),A=(posterior$best$inf),B=(posterior$max_infect150$inf),cols=colors,from=0,to=1,main="",xlab=expression(nu),xaxt="n")
dev.off()


#### FIGURE 8
for(d in nvl){
    png(paste0("5e8779c56517890001536101/figures/posterior2d_",d,".png"),pointsize=30, res = 100, width = 25, height = 6, units = "in")
    best=allresults[rank(allresults$distances)<d,]
    par(mfrow=c(1,6))
    par(mar=c(4,4,1,1))
    hdr.boxplot.2d(best$inf_r,log10(best$sat_r),prob=seq(20,100,10),shadecols=alpha(colorbest,.9),xlim=c(0,1),ylim=c(-1,3),xlab=expression(nu[r]),ylab=expression(kappa[r]))
    hdr.boxplot.2d(best$inf,log10(best$sat),prob=seq(20,100,10),shadecols=alpha(colorbest,.9),xlim=c(0,1),ylim=c(-1,3),xlab=expression(nu[r]),ylab=expression(kappa))
    hdr.boxplot.2d(best$inf,best$pind,prob=seq(20,100,10),shadecols=alpha(colorbest,.9),xlim=c(0,1),ylim=c(0,1),xlab=expression(nu[r]),ylab="individual learning")
    hdr.boxplot.2d(log10(best$sat),best$pind,prob=seq(20,100,10),shadecols=alpha(colorbest,.9),xlim=c(-1,3),ylim=c(0,1),xlab=expression(kappa),ylab="individual learning")
    hdr.boxplot.2d(best$sl_rad,best$pind,prob=seq(20,100,10),shadecols=alpha(colorbest,.9),xlim=c(1,100),ylim=c(0,1),xlab="radius social learning",ylab="individual learning")
    hdr.boxplot.2d(best$sl_rad,best$inf,prob=seq(20,100,10),shadecols=alpha(colorbest,.9),xlim=c(1,100),ylim=c(0,1),xlab="radius social learning",ylab=expression(nu))
    dev.off()
}

    hdr.boxplot.2d(posterior$best$inf_r,log10(posterior$best$sat_r),prob=seq(75,90,10),shadecols=alpha(colorbest,.9),xlim=c(0,1),ylim=c(-1,3),xlab=expression(nu[r]),ylab="revert steepness (log10)")
    hdr.boxplot.2d(posterior$best$inf,log10(posterior$best$sat),prob=seq(75,90,10),shadecols=alpha(colorbest,.9),xlim=c(0,1),ylim=c(-1,3),xlab=expression(nu[r]),ylab="steepness (log10)")

    hdr.boxplot.2d(posterior$max_infect150$inf_r,log10(posterior$max_infect150$sat_r),prob=seq(20,100,10),shadecols=alpha(colorbest,.9),xlim=c(0,1),ylim=c(-1,3),xlab=expression(nu[r]),ylab="revert steepness (log10)")
    hdr.boxplot.2d(posterior$max_infect150$inf,log10(posterior$max_infect150$sat),prob=seq(20,100,10),shadecols=alpha(colorbest,.9),xlim=c(0,1),ylim=c(-1,3),xlab=expression(nu[r]),ylab="steepness (log10)")




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
plot(test$time_max,test$max_infect,col=alpha(color.gradient(log(test$sat),c(mygreen,myred)),.2),main=expression(kappa),pch=20,xlab=expression(tau),ylab=expression(I[max]))
legend("topright",legend=paste("steepness=10^",round(log10(unique(test$sat)[c(1,5,10)]))),col=clrssat[c(1,5,10)],lty=1)
plot(test$time_max,test$max_infect,col=alpha(color.gradient(test$inf,c(mygreen,myred)),.2),main=expression(nu),pch=20,xlab=expression(tau),ylab=expression(I[max]))
legend("topright",legend=paste("point=",round(unique(test$inf)[c(1,5,10)],digit=1)),col=clrssat[c(1,5,10)],lty=1)
plot(test$time_max,test$max_infect,col=alpha(color.gradient(test$pind,c(mygreen,myred)),.2),main="indiv learning",pch=20,xlab=expression(tau),ylab=expression(I[max]))
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
plot(test$time_max,test$max_infect,col=alpha(color.gradient(log(test$sat),c(mygreen,"yellow", myred)),.2),main=expression(kappa),pch=20,xlab=expression(tau),ylab=expression(I[max]))
sset=sort(log10(test$sat))
cols=alpha(color.gradient(sset,c(mygreen,"yellow", myred)),1)
lsset=seq(1,length(sset),length.out=4)
legend("topright",legend=paste("steepness=10^",round(sset[lsset])),col=cols[lsset],pch=20)
plot(test$time_max,test$max_infect,col=alpha(color.gradient(log(test$sat_r),c(mygreen,"yellow", myred)),.2),main=expression(kappa[r]),pch=20,xlab=expression(tau),ylab=expression(I[max]))
legend("topright",legend=paste("steepness=10^",round(sset[lsset])),col=cols[lsset],pch=20)
plot(test$time_max,test$max_infect,col=alpha(color.gradient(test$inf,c(mygreen,"yellow", myred)),.2),main=expression(nu),pch=20,xlab=expression(tau),ylab=expression(I[max]))
sset=sort(test$inf)
cols=alpha(color.gradient(sset,c(mygreen,"yellow", myred)),1)
lsset=seq(1,length(sset),length.out=4)
legend("topright",legend=paste("inp",round(sset[lsset],digit=2)),col=cols[lsset],pch=20)
plot(test$time_max,test$max_infect,col=alpha(color.gradient(test$inf_r,c(mygreen,"yellow", myred)),.2),main="revert inif. point",pch=20,xlab=expression(tau),ylab=expression(I[max]))
legend("topright",legend=paste("inp",round(sset[lsset],digit=2)),col=cols[lsset],pch=20)
plot(test$time_max,test$max_infect,col=alpha(color.gradient(test$pind,c(mygreen,"yellow", myred)),.2),main="individual learning",pch=20,xlab=expression(tau),ylab=expression(I[max]))
sset=sort(test$pind)
cols=alpha(color.gradient(sset,c(mygreen,"yellow", myred)),1)
lsset=seq(1,length(sset),length.out=4)
legend("topright",legend=paste("ind. learn.",round(sset[lsset],digit=2)),col=cols[lsset],pch=20)
plot(test$time_max,test$max_infect,col=alpha(color.gradient(test$distances,c(mygreen,"yellow", myred)),.2),main="distance",pch=20,xlab=expression(tau),ylab=expression(I[max]))
sset=sort(test$distances)
cols=alpha(color.gradient(sset,c(mygreen,"yellow", myred)),1)
lsset=seq(1,length(sset),length.out=4)
legend("topright",legend=paste("distance",round(sset[lsset],digit=2)),col=cols[lsset],pch=20)
dev.off()



for(d in c(.8,.6,.4,.2)){
    png(paste0("5e8779c56517890001536101/figures/posterior2d_",d*10,"_",name,".png"),width=1100,height=333,pointsize=20)
    best=allresults[allresults$distances<d,]
    par(mfrow=c(1,6))
    hdr.boxplot.2d(best$inf_r,log10(best$sat_r),prob=seq(20,100,10),shadecols=mygreen,xlim=c(0,1),ylim=c(-1,3),xlab=expression(nu[r]),ylab="rever steepness")
    hdr.boxplot.2d(best$inf,log10(best$sat),prob=seq(20,100,10),shadecols="yellow",xlim=c(0,1),ylim=c(-1,3),xlab=expression(nu[r]),ylab=expression(kappa))
    hdr.boxplot.2d(best$inf,best$pind,prob=seq(20,100,10),shadecols=myred,xlim=c(0,1),ylim=c(0,1),xlab=expression(nu[r]),ylab="individual learning")
    hdr.boxplot.2d(log10(best$sat_r),best$pind,prob=seq(20,100,10),shadecols=mygreen,xlim=c(-1,3),ylim=c(0,1),xlab=expression(kappa),ylab="individual learning")
    hdr.boxplot.2d(best$sl_rad,best$pind,prob=seq(20,100,10),shadecols=myred,xlim=c(1,50),ylim=c(0,1),xlab="radius social learning",ylab="individual learning")
    hdr.boxplot.2d(best$sl_rad,log10(best$sat),prob=seq(20,100,10),shadecols=myred,xlim=c(0,100),ylim=c(-1,3),xlab="radius social learning",ylab=expression(kappa))
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
plot(test$time_max,test$max_infect,col=alpha(color.gradient(log(test$sat),c(mygreen,"yellow", myred)),.2),main=expression(kappa),pch=20,xlab=expression(tau),ylab=expression(I[max]))
sset=sort(log10(test$sat))
cols=alpha(color.gradient(sset,c(mygreen,"yellow", myred)),1)
lsset=seq(1,length(sset),length.out=4)
legend("topright",legend=paste("steepness=10^",round(sset[lsset])),col=cols[lsset],pch=20)
plot(test$time_max,test$max_infect,col=alpha(color.gradient(log(test$sat_r),c(mygreen,"yellow", myred)),.2),main=expression(kappa[r]),pch=20,xlab=expression(tau),ylab=expression(I[max]))
legend("topright",legend=paste("steepness=10^",round(sset[lsset])),col=cols[lsset],pch=20)
plot(test$time_max,test$max_infect,col=alpha(color.gradient(test$inf,c(mygreen,"yellow", myred)),.2),main=expression(nu),pch=20,xlab=expression(tau),ylab=expression(I[max]))
sset=sort(test$inf)
cols=alpha(color.gradient(sset,c(mygreen,"yellow", myred)),1)
lsset=seq(1,length(sset),length.out=4)
legend("topright",legend=paste("inp",round(sset[lsset],digit=2)),col=cols[lsset],pch=20)
plot(test$time_max,test$max_infect,col=alpha(color.gradient(test$inf_r,c(mygreen,"yellow", myred)),.2),main="revert inif. point",pch=20,xlab=expression(tau),ylab=expression(I[max]))
legend("topright",legend=paste("inp",round(sset[lsset],digit=2)),col=cols[lsset],pch=20)
plot(test$time_max,test$max_infect,col=alpha(color.gradient(test$pind,c(mygreen,"yellow", myred)),.2),main="individual learning",pch=20,xlab=expression(tau),ylab=expression(I[max]))
sset=sort(test$pind)
cols=alpha(color.gradient(sset,c(mygreen,"yellow", myred)),1)
lsset=seq(1,length(sset),length.out=4)
legend("topright",legend=paste("ind. learn.",round(sset[lsset],digit=2)),col=cols[lsset],pch=20)
plot(test$time_max,test$max_infect,col=alpha(color.gradient(test$distances,c(mygreen,"yellow", myred)),.2),main="distance",pch=20,xlab=expression(tau),ylab=expression(I[max]))
sset=sort(test$distances)
cols=alpha(color.gradient(sset,c(mygreen,"yellow", myred)),1)
lsset=seq(1,length(sset),length.out=4)
lapply(c(.3,.6),function(s){t=200:1000;lines(t,m(t=t,Rm=max(test$max_infect)-min(test$max_infect),mm=min(test$max_infect),Rt=max(test$time_max)-min(test$time_max),mt=min(test$time_max),s=s),col="grey",lwd=2)})
legend("topright",legend=paste("distance",round(sset[lsset],digit=2)),col=cols[lsset],pch=20)
dev.off()

layout(matrix(c(1,1,1,1,2,3,4,5),nrow=4,ncol=2),width=c(2,1))  



    png(paste0("5e8779c56517890001536101/figures/posterior2d_inf_r_sat_r.png"),pointsize=17, res = 100, width = 7, height = 7, units = "in")
    hdr.boxplot.2d(best$inf_r,log10(best$sat_r),prob=seq(20,100,10),shadecols=alpha(colorbest,.9),xlim=c(0,1),ylim=c(-1,3),xlab=expression(nu[r]),ylab=expression(kappa[r]))
    dev.off()
    png(paste0("5e8779c56517890001536101/figures/posterior2d_inf_sat_r.png"),pointsize=17, res = 100, width = 7, height = 7, units = "in")
    hdr.boxplot.2d(best$inf,log10(best$sat_r),prob=seq(20,100,10),shadecols=alpha(colorbest,.9),xlim=c(0,1),ylim=c(-1,3),xlab=expression(nu),ylab=expression(kappa[r]))
dev.off()


for(d in c(1,.9,.8,.6,.4,.3)){
    png(paste0("5e8779c56517890001536101/figures/posterior2d_",d*10,".png"),width=1600,height=333,pointsize=20)
    best=allresults[allresults$distances<d,]
    par(mfrow=c(1,5))
    hdr.boxplot.2d(best$inf_r,log10(best$sat_r),prob=seq(20,100,10),shadecols="grey",xlim=c(0,1),ylim=c(-1,3),xlab=expression(nu[r]),ylab=expression(kappa[r]))
    hdr.boxplot.2d(best$inf,log10(best$sat),prob=seq(20,100,10),shadecols="grey",xlim=c(0,1),ylim=c(-1,3),xlab=expression(nu[r]),ylab=expression(kappa))
    hdr.boxplot.2d(best$inf,best$pind,prob=seq(20,100,10),shadecols="grey",xlim=c(0,1),ylim=c(0,1),xlab=expression(nu[r]),ylab="individual learning")
    hdr.boxplot.2d(log10(best$sat),best$pind,prob=seq(20,100,10),shadecols="grey",xlim=c(-1,3),ylim=c(0,1),xlab=expression(kappa),ylab="individual learning")
    hdr.boxplot.2d(best$sl_rad,best$pind,prob=seq(20,100,10),shadecols="grey",xlim=c(1,100),ylim=c(0,1),xlab=expression(nu[r]),ylab="individual learning")
    dev.off()
}

png("5e8779c56517890001536101/figures/meanPosteriors.png",pointsize=17, res = 300, width = 7, height = 7, units = "in")
par(mar=c(4,4,2,2))
plot(1,1,type="n",ylim=c(0,500),xlim=c(1,1500),ylab="",xlab="",xaxt="n",yaxt="n")
legend("topright",legend=c("No Learning (all NA)","No Learning (all A)","Worst Simulations","Best Simulations"),col=c(myred,mygreen),lty=c(1,1,2,2),lwd=2,cex=.95)
axis(1)
axis(2)
mtext("time",1,3)
mtext("number of infected people",2,3)
meanBad=apply(neutralGood,1,mean)
meanGood=apply(neutralGoodG,1,mean)
meanWorst=apply(repetWorst,1,mean)
meanBest=apply(repetBest,1,mean)

lines(meanGood,col=alpha(mygreen,.9),lwd=2)
lines(meanBad,col=alpha(myred,.9),lwd=2)
lines(meanBest,col=alpha(mygreen,.9),lwd=2,lty=2)
lines(meanWorst,col=alpha(myred,.9),lwd=2,lty=2)
dev.off()

#' plot posteriors distribution against priors
#'@param A : a vector with posterior
#'@param B : a vector with posterior 
#'@param prior : a vector with posterior 
#' @export plot2dens 
plot2dens <- function(A=NULL,B=NULL,C=NULL,from=NULL,to=NULL,prior=NULL,cols=c(alpha("red",.8),alpha("blue",.8),alpha("yellow",.8)),hdr=F,yaxt=NULL,log=F,...){

    denseP=NULL
    denseA=NULL
    denseB=NULL
    denseC=NULL
    if(!is.null(prior))prior=prior[!is.na(prior)]
    if(is.null(yaxt))yaxt="n"
    if(is.null(from))from=min(A,B,prior)
    if(is.null(to))to=max(A,B,prior)
    if(!is.null(A))denseA=density(A,from=from,to=to)
    if(!is.null(B))denseB=density(B,from=from,to=to)
    if(!is.null(C))denseC=density(C,from=from,to=to)
    if(length(prior)==2)denseP=density(runif(100000,prior[1],prior[2]),from=from,to=to)
    else if(!is.null(prior))denseP=density(prior,from=from,to=to)
    print("donde")

    if(is.null(names(cols)))names(cols)=c("P","A","B")
    rangex=range(denseB$x,denseA$x,denseP$x,denseC$x)
    rangey=range(0,denseB$y,denseA$y,denseP$y,denseC$y)
    stepy=max(rangey)*0.2
    if(hdr)miny=-1*stepy else miny=0
    plot(denseA,ylim=c(miny,max(rangey)),xlim=rangex,type="n",xaxt="n",yaxt=yaxt,ylab="",...)
    if(log){
        xaxp=par()$xaxp
        axis(1,at=seq(xaxp[1],xaxp[2],length.out=xaxp[3]+1),label=10^seq(xaxp[1],xaxp[2],length.out=xaxp[3]+1))
    }
    else axis(1)
    mtext("Density",2,1)
    if(!is.null(C))
        polygon(c(from,denseC$x,to),c(0,denseC$y,0),col=cols["C"],lwd=2)
    if(!is.null(A)){
        #polygon(c(from,denseA$x,to),c(0,denseA$y,0),col="white",border=NA,lwd=2)
        if(hdr){
            hdstaA=hdr(A,prob=c(75,95),lambda=0.9)
            hdrA=hdstaA$hdr
            print(hdrA)
            polygon(c(hdrA[1,1],hdrA[1,1],hdrA[1,2],hdrA[1,2]),c(-.2*stepy,-.9*stepy,-.9*stepy,-.2*stepy),col=cols["A"],lwd=1)

            polygon(c(hdrA[2,1],hdrA[2,1],hdrA[2,2],hdrA[2,2]),c(-.2*stepy,-.9*stepy,-.9*stepy,-.2*stepy),col=cols["A"],lwd=1)
            polygon(c(hdrA[2,1],hdrA[2,1],hdrA[2,2],hdrA[2,2]),c(-.2*stepy,-.9*stepy,-.9*stepy,-.2*stepy),col=cols["A"],lwd=1)
            segments(hdstaA$mode,-.3*stepy,hdstaA$mode,-.8*stepy)
            middle=(-.2*stepy+-.9*stepy)/2
            segments(min(A),middle,hdrA[1,1],middle)
            segments(hdrA[1,2],middle,max(A),middle)
            polygon(c(from,denseA$x,to),c(0,denseA$y,0),col=cols["A"],lwd=2)
        }
        else
            polygon(c(from,denseA$x,to),c(0,denseA$y,0),col=cols["A"],lwd=2)
    }
    if(!is.null(B)){
        polygon(c(from,denseB$x,to),c(0,denseB$y,0),col="white",lwd=2,density=NA,border=0)
        polygon(c(from,denseB$x,to),c(0,denseB$y,0),col=cols["B"],lwd=2,density=5,border=1)
    }
    if(!is.null(prior))
        polygon(c(from,denseP$x,to),c(0,denseP$y,0),col=cols["P"],lwd=2)

    #if(!is.null(A))
    #        polygon(c(from,denseA$x,to),c(0,denseA$y,0),col=NA,lwd=2)

}
