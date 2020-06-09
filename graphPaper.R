m<-function(t,Rm,Rt,mm,mt,s)(Rm*(2*s-1+(t-mt)/Rt))+mm 

source("abmEpi.R")
library(RColorBrewer)
library(scales)
library(hdrcde)

color_class=rev(brewer.pal(5,"PuBu"))
colorbest=color_class[1]
myred=rgb(r=213,g=94,b=0,maxColorValue=255)
myred="red"
mygreen=rgb(r=0,g=158,b=115,maxColorValue=255)
mygreen=colorbest


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
        polygon(c(from,denseB$x,to),c(0,denseB$y,0),col=alpha("white",.6),lwd=2,density=NA,border=0)
        polygon(c(from,denseB$x,to),c(0,denseB$y,0),col=cols["B"],lwd=2,density=5,border=1)
    }
    if(!is.null(prior))
        polygon(c(from,denseP$x,to),c(0,denseP$y,0),col=cols["P"],lwd=2)

    #if(!is.null(A))
    #        polygon(c(from,denseA$x,to),c(0,denseA$y,0),col=NA,lwd=2)

}

plotQuartiles <- function(allres,ylim=NULL,...){
    clrssat=colorRampPalette(c("blue","red"))(10)
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

pdf("5e8779c56517890001536101/figures/meanTrajIllu.pdf",pointsize=17, width = 7, height = 7)
#png("5e8779c56517890001536101/figures/meanTrajIllu.png",pointsize=17, res = 300, width = 7, height = 7, units = "in")
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

    pdf("5e8779c56517890001536101/figures/illuSig.pdf",pointsize=25,  width = 15, height = 7)
    par(mfrow=c(1,2))
    par(mar=c(4,4,1,1))
    x=seq(0,1,.01)   
    inp=seq(0,1,.05)   
    clrs=colorRampPalette(c("purple4","yellow"))(length(inp))
    #plot(x,sig(x),type="n",ylim=c(0,1),xlim=c(0,1),ylab=expression(sig(x,st=10,inp)),main=expression(inp %in%  "(" * list(0,1) * ")") )
    plot(x,sig(x),type="n",ylim=c(0,1),xlim=c(0,1),xlab="% infected",ylab=expression(P[switch]),main="" )
    for(i in 1:length(inp)) lines(x,sig(x,b=inp[i],a=4),ylim=c(0,1),xlim=c(0,1),col=clrs[i],lwd=2) 
    leg=seq(1,length(inp),length.out=5)
    legend("bottomright",legend=rev(sapply(inp[leg],function(d)as.expression(bquote(.(d))))),title=expression(nu),col=rev(clrs[leg]),lwd=2,cex=.8,bg="white")

    tstp=25
    stp=rev(seq(-3,3,length.out=tstp))
    clrs=colorRampPalette(c("purple4","yellow"))(tstp)
    #plot(x,sig(x),type="n",ylim=c(0,1),xlim=c(0,1),ylab="P(B->G)~sig(x,st,inp=.5)",main=bquote(stp %in% "(" * list(10^.(stp[1]),10^.(stp[tstp])) * ")") )
    plot(x,sig(x),type="n",ylim=c(0,1),xlim=c(0,1),xlab="% infected",ylab=expression(P[switch]),main="" )
    for(i in rev(1:length(stp))) lines(x,sig(x,a=10^stp[i]),ylim=c(0,1),xlim=c(0,1),col=clrs[i],lwd=2) 
    leg=seq(1,tstp,length.out=5)
    #legend("bottomright",legend=sapply(round(stp[leg]),function(d)as.expression(bquote(kappa==.(10^d)))),col=clrs[leg],lwd=2)
    legend("bottomright",legend=sapply(round(stp[leg]),function(d)as.expression(bquote(.(10^d)))),col=clrs[leg],lwd=2,title=expression(kappa),cex=.8,bg="white")
    dev.off()

### Figure3
subresults=allresults[sample(nrow(allresults),20000),]
pdf("5e8779c56517890001536101/figures/distribSimulation.pdf",pointsize=25,  width = 11, height = 11)
#png("5e8779c56517890001536101/figures/distribSimulation.png",pointsize=25, res = 300, width = 11, height = 11, units = "in")
rank=rank(subresults$distances)
orderscol=rank
color_class=rev(brewer.pal(5,"PuBu"))
colorbest=color_class[1]
orderscol=rep(alpha(color_class[5],.4),length(orderscol))
nvl=c(1000,2500,5000,10000)
for(i in (length(nvl):1)) orderscol[rank<=nvl[i]]=alpha(color_class[i],.4)
#plot(subresults$time_max,subresults$max_infect,bg=orderscol,main=expression(I[max]*" wrt "*tau),pch=21,xlab=expression(tau),ylab=expression(I[max]),lwd=.1)
plot(subresults$time_max,subresults$max_infect,bg=orderscol,main="",pch=21,xlab=expression(tau),ylab=expression(I[max]),lwd=.1)
legend("topright",legend=c(paste("<",nvl),"all"),pt.bg=color_class,pch=21,col=alpha(1,.5),pt.lwd=.1,title="Rank")
dev.off()

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
pdf(paste0("5e8779c56517890001536101/figures/switchfunctions.pdf"),pointsize=17,  width = 13, height = 7)
#png(paste0("5e8779c56517890001536101/figures/switchfunctions.png"),pointsize=17, res = 100, width = 13, height = 7, units = "in")
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

pdf("5e8779c56517890001536101/figures/posterior_pind.pdf",pointsize=30,  width = 10, height = 10)
par(mar=mar)
 plot2dens(prior=1-allresults$pind,A=1-posterior$best$pind,B=1-posterior$max_infect250$pind,cols=colors,from=0,to=1,main="posterior proba social learning",xlab="proba social learning")
dev.off()

pdf("5e8779c56517890001536101/figures/posterior_sl_rad.pdf",pointsize=30,  width = 10, height = 10)
par(mar=mar)
 plot2dens(prior=allresults$sl_rad,A=posterior$best$sl_rad,B=posterior$max_infect150$sl_rad,cols=colors,from=1,to=100,main="",xlab="radius social learning")
dev.off()

pdf("5e8779c56517890001536101/figures/posterior_sat_r.pdf",pointsize=30,  width = 10, height = 10)
par(mar=mar)
 plot2dens(prior=log10(allresults$sat_r),A=log10(posterior$best$sat_r),B=log10(posterior$time_max150$sat_r),cols=colors,from=-1,to=3,main="",xlab=expression(kappa[r]),log=T)
dev.off()

pdf("5e8779c56517890001536101/figures/posterior_sat.pdf",pointsize=30,  width = 10, height = 10)
par(mar=mar)
 plot2dens(prior=log10(allresults$sat),A=log10(posterior$best$sat),B=log10(posterior$max_infect150$sat),cols=colors,from=-1,to=3,xlab=expression(kappa),main="",xaxt="n",log=T)
dev.off()

pdf("5e8779c56517890001536101/figures/posterior_inf_r.pdf",pointsize=30,  width = 10, height = 10)
par(mar=mar)
 plot2dens(prior=(allresults$inf_r),A=(posterior$best$inf_r),B=(posterior$max_infect150$inf_r),cols=colors,from=0,to=1,main="posterior revert inflection point",xlab=expression(nu[r]),xaxt="n")
legend("topleft",legend=c("prior","final time","150 timesteps"),fill=c(0,colorbest,colorbest),density=c(NA,NA,20))
dev.off()

pdf("5e8779c56517890001536101/figures/posterior_inf.pdf",pointsize=30,  width = 10, height = 10)
par(mar=mar)
 plot2dens(prior=(allresults$inf),A=(posterior$best$inf),B=(posterior$max_infect150$inf),cols=colors,from=0,to=1,main="posterior inflection point",xlab=expression(nu),xaxt="n")
dev.off()


##FIGURE 7BIS
mar=c(4,3,2,1)

colors=c(P="NA",A=alpha(colorbest,1),B=alpha(colorbest,1))
pdf("5e8779c56517890001536101/figures/worst_pind.pdf",pointsize=30,  width = 10, height = 10)
par(mar=mar)
 plot2dens(prior=1-posterior$worst$pind,A=1-posterior$best$pind,B=1-posterior$max_infect150$pind,cols=colors,from=0,to=1,main="",xlab="proba social learning")
dev.off()

pdf("5e8779c56517890001536101/figures/worst_sl_rad.pdf",pointsize=30,  width = 10, height = 10)
par(mar=mar)
 plot2dens(prior=posterior$worst$sl_rad,A=posterior$best$sl_rad,B=posterior$max_infect150$sl_rad,cols=colors,from=1,to=100,main="posterior proba radius social learning",xlab="radius social learning")
dev.off()

pdf("5e8779c56517890001536101/figures/worst_sat_r.pdf",pointsize=30,  width = 10, height = 10)
par(mar=mar)
 plot2dens(prior=log10(posterior$worst$sat_r),A=log10(posterior$best$sat_r),B=log10(posterior$time_max150$sat_r),cols=colors,from=-1,to=3,main="posterior revert steepness",xlab=expression(kappa[r]),log=T)
dev.off()

pdf("5e8779c56517890001536101/figures/worst_sat.pdf",pointsize=30,  width = 10, height = 10)
par(mar=mar)
 plot2dens(prior=log10(posterior$worst$sat),A=log10(posterior$best$sat),B=log10(posterior$max_infect150$sat),cols=colors,from=-1,to=3,xlab=expression(kappa),main="posterior steepness",xaxt="n",log=T)
dev.off()

pdf("5e8779c56517890001536101/figures/worst_inf_r.pdf",pointsize=30,  width = 10, height = 10)
par(mar=mar)
 plot2dens(prior=(posterior$worst$inf_r),A=(posterior$best$inf_r),B=(posterior$max_infect150$inf_r),cols=colors,from=0,to=1,main="",xlab=expression(nu[r]),xaxt="n")
legend("topleft",legend=c("prior","final time","150 timesteps"),fill=c(0,colorbest,colorbest),density=c(NA,NA,20))
dev.off()

pdf("5e8779c56517890001536101/figures/worst_inf.pdf",pointsize=30,  width = 10, height = 10)
par(mar=mar)
 plot2dens(prior=(posterior$worst$inf),A=(posterior$best$inf),B=(posterior$max_infect150$inf),cols=colors,from=0,to=1,main="",xlab=expression(nu),xaxt="n")
dev.off()


#### FIGURE 8
for(d in nvl){
    pdf(paste0("5e8779c56517890001536101/figures/posterior2d_",d,".pdf"),pointsize=30,  width = 25, height = 6)
    best=allresults[rank(allresults$distances)<d,]
    par(mfrow=c(1,6))
    par(mar=c(4,4,1,1))
    hdr.boxplot.2d(best$inf_r,log10(best$sat_r),prob=seq(20,100,10),shadecols=alpha(colorbest,.9),xlim=c(0,1),ylim=c(-1,3),xlab=expression(nu[r]),ylab=expression(kappa[r]))
    hdr.boxplot.2d(best$inf,log10(best$sat),prob=seq(20,100,10),shadecols=alpha(colorbest,.9),xlim=c(0,1),ylim=c(-1,3),xlab=expression(nu[r]),ylab=expression(kappa))
    hdr.boxplot.2d(best$inf,best$pind,prob=seq(20,100,10),shadecols=alpha(colorbest,.9),xlim=c(0,1),ylim=c(0,1),xlab=expression(nu[r]),ylab="p")
    hdr.boxplot.2d(log10(best$sat),best$pind,prob=seq(20,100,10),shadecols=alpha(colorbest,.9),xlim=c(-1,3),ylim=c(0,1),xlab=expression(kappa),ylab="p")
    hdr.boxplot.2d(best$sl_rad,best$pind,prob=seq(20,100,10),shadecols=alpha(colorbest,.9),xlim=c(1,100),ylim=c(0,1),xlab="r",ylab="p")
    hdr.boxplot.2d(best$sl_rad,best$inf,prob=seq(20,100,10),shadecols=alpha(colorbest,.9),xlim=c(1,100),ylim=c(0,1),xlab="r",ylab=expression(nu))
    dev.off()
}





save(file="worst.bin",worst)
save(file="best.bin",best)



best=readRDS("best.bin")
worst=readRDS("worst.bin")
xsize=ysize=100
library(parallel)
cl <- makeForkCluster(6,outfile="")
neutralGood=parSapply(cl,1:500,function(i){
                      print(i);
                         simu=abmSIR(500,1500,p=c(1,1),di=2,i0=1,recovery=c(8,14)*25,speed=c(1,.2),xsize=xsize,ysize=ysize,
                                     inf=9,
                                     sat=1000,
                                     p_i=1,
                                     strategy="best",
                                     ts=T,ap=F,visu=F
                                     )$timeseries[,2]

})
save(file="neutralGood.bin",neutralGood)

neutralBad=parSapply(cl,1:500,function(i){
                     print(i);
                     simu=abmSIR(500,1500,p=c(.1,.1),di=2,i0=1,recovery=c(8,14)*25,speed=c(1,.2),xsize=xsize,ysize=ysize,
                                 inf=9,
                                 sat=1000,
                                 p_i=1,
                                 strategy="best",
                                 ts=T,ap=F,visu=F
                                 )$timeseries[,2]

})
save(file="neutralBad.bin",neutralBad)


repetBest=parSapply(cl,1:500,function(i){
                    j=sample(nrow(best),1);print(paste(i,j));
                    simu=abmSIR(500,1500,p=c(1,.1),di=2,i0=1,recovery=c(8,14)*25,speed=c(1,.2),xsize=xsize,ysize=ysize,
                                inf=best$inf[j],
                                sat=best$sat[j],
                                inf_r=best$inf_r[j],
                                sat_r=best$sat_r[j],
                                p_i=best$pind[j],
                                sl_rad=best$sl_rad[j],
                                strategy="best",
                                ts=T,ap=F,visu=F
                                )$timeseries[,2]
})
save(file="repetBest.bin",repetBest)

repetWorst=parSapply(cl,1:500,function(i){
                    j=sample(nrow(worst),1);print(paste(i,j));
                    simu=abmSIR(500,1500,p=c(1,.1),di=2,i0=1,recovery=c(8,14)*25,speed=c(1,.2),xsize=xsize,ysize=ysize,
                                inf=worst$inf[j],
                                sat=worst$sat[j],
                                inf_r=worst$inf_r[j],
                                sat_r=worst$sat_r[j],
                                p_i=worst$pind[j],
                                sl_rad=worst$sl_rad[j],
                                strategy="best",
                                ts=T,ap=F,visu=F
                                )$timeseries[,2]
})
save(file="repetWorst.bin",repetWorst)


stopCluster(cl)
pop=generatePopulation(N=500,xsize=xsize,ysize=ysize,recovery=c(8,14)*25,speed=c(1,.2),behavior=rep(G,500))
cl <- makeForkCluster(6,outfile="")
neutralGoodG=parSapply(cl,1:500,function(i){
                      print(i);
                         simu=abmSIR(1500,p=c(.1,.1),di=2,i0=1,xsize=xsize,ysize=ysize,
                                     pop=pop,
                                     inf=9,
                                     sat=1000,
                                     p_i=1,
                                     strategy="random",
                                     ts=T,ap=F,visu=F
                                     )$timeseries[,2]

})
save(file="neutralGoodG.bin",neutralGoodG)
stopCluster(cl)

>>>>>>> 4749fe46c3b6bb2ea33162ed569370e5c9f464d1

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

pdf("5e8779c56517890001536101/figures/meanPosteriors150.pdf",pointsize=17,  width = 7, height = 7)
par(mar=c(4,4,2,2))
plot(1,1,type="n",ylim=c(0,500),xlim=c(1,1500),ylab="",xlab="",xaxt="n",yaxt="n")
abline(v=150,lwd=2,col="gray")
par(xpd=NA)
text(150,-20,"t=150",col="1",pos=1)
par(xpd=F)
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

