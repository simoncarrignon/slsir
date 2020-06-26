ior
##### Grpah papers
##get data
library(RColorBrewer)
modelnames=c("midCurveAllBad","midCurveAllBadBestSLS","midCurve15Good","midCurve15GoodBestSLS","midCurve100GoodBestSLS","burnin100BestSLS")
names(modelnames)=modelnames
modelcol=brewer.pal(length(modelnames),name="Dark2")
names(modelcol)=modelnames
modelsubnames=sapply(modelnames,function(mod)paste0(strsplit(mod,"")[[1]][-(1:8)],collapse=""))
modelsubnames=LETTERS[1:6]


#name="testradius"

allexpe=lapply(modelnames,function(name){
aa=lapply(unlist(lapply(list.dirs(name),list.files,pattern="allresults.bin",full.names=T)),function(f){load(f);return(allresults)})
allresults=do.call("rbind",aa)
if(sum(which(allresults$max_infect <4)>2))allresults=allresults[-which(allresults$max_infect <4),]
allresults$distances=(1-bdistance(allresults$time_max) + bdistance(allresults$max_infect))/2
allresults$distances2=(1-bdistance(allresults$time_max150) + bdistance(allresults$max_infect150))/2
allresults$distances3=(1-bdistance(allresults$time_max250) + bdistance(allresults$max_infect250))/2
return(allresults)})
nsims=sapply(allexpe,nrow)

#Bayes factors
allscores=c("distances","distances2","distances3","time_max","time_max150","time_max250","max_infect","max_infect150","max_infect250")
names(allscores)=allscores
tablesscores=sapply(allscores,function(as){
                    allscores=lapply(allexpe,"[[",as)
                    allnames=rep(names(allscores),lengths(allscores))
                    table(allnames[order(unlist(allscores),decreasing=F)][1:(10000)])
})
normalized=apply(tablesscores,2,"/",nsims[rownames(tablesscores)])
dis=normalized[,1]
##calculate bayes factor, by row, uie mat_rat[i,j] = dis[i]/dis[j]
mat_rat=sapply(1:length(dis),function(i)sapply(1:length(dis),function(j)dis[i]/dis[j]))  
colnames(mat_rat)=names(dis)
rownames(mat_rat)=names(dis)

### Figure3 (for all
par(mfrow=c(2,2))
lapply(allexpe,function(allresults){
       #png("5e8779c56517890001536101/figures/distribSimulation.png",pointsize=17, res = 100, width = 7, height = 7, units = "in")
       plot(allresults$time_max,allresults$max_infect,col=alpha(color.gradient(allresults$distances,c(mygreen,"yellow", myred)),.2),main="distance",pch=20,xlab="Time Max",ylab="Max Infected")
       sset=sort(allresults$distances)
       cols=alpha(color.gradient(sset,c(mygreen,"yellow", myred)),1)
       lsset=seq(1,length(sset),length.out=4)
       legend("topright",legend=sapply(round(sset[lsset],digit=2),function(d)as.expression(bquote(delta==.(d)))),col=cols[lsset],pch=20)
})
#dev.off()



## distribution distances/metrcis/4models
pdf(paste0("5e8779c56517890001536101/figures/distributionMetricsAllModels_fulltime.pdf"),pointsize=25, width = 13, height = 5)
par(mfrow=c(1,3),oma=c(0,2,0,0),mar=c(4,1,1,0))
for(v in c(1,4,7)){
    var=allscores[v]
        xlimits=range(sapply(allexpe,"[[",var))
        alld=lapply(allexpe,function(i,var)density(i[[var]]),var=var)
        ylimits=range(sapply(alld,"[[","y"))
        t=""
        if(v==1)t=expression(delta)
        if(v==4)t=expression(tau)
        if(v==7)t=expression(I[max])
        plot(1,1,ylim=ylimits,xlim=xlimits,type="n",main=t,yaxt="n",xlab=t,ylab="density")
        lapply(names(alld),function(i)lines(alld[[i]],col=modelcol[i],lwd=2))

}
legend("topleft",legend=modelsubnames,col=modelcol,lwd=3,lty=1)
mtext("density",2,0,outer=T)
dev.off()
png(paste0("distributionMetricsAllModels_150.png"),pointsize=15, res = 100, width = 13, height = 5, units = "in")
par(mfrow=c(1,3))
for(var in allscores[c(7)+1]){
        xlimits=range(sapply(allexpe,"[[",var))
        alld=lapply(allexpe,function(i,var)density(i[[var]]),var=var)
        ylimits=range(sapply(alld,"[[","y"))
        plot(1,1,ylim=ylimits,xlim=xlimits,type="n",main=var,yaxt="n",xlab=var,ylab="density")
        lapply(names(alld),function(i)lines(alld[[i]],col=modelcol[i],lwd=3))

}
legend("topright",legend=modelsubnames,col=modelcol,lwd=3,lty=1)
dev.off()
png(paste0("distributionMetricsAllModels_250.png"),pointsize=15, res = 100, width = 13, height = 5, units = "in")
par(mfrow=c(1,3))
for(var in allscores[c(1,4,7)+2]){
        xlimits=range(sapply(allexpe,"[[",var))
        alld=lapply(allexpe,function(i,var)density(i[[var]]),var=var)
        ylimits=range(sapply(alld,"[[","y"))
        plot(1,1,ylim=ylimits,xlim=xlimits,type="n",main=var,yaxt="n",xlab=var,ylab="density")
        lapply(names(alld),function(i)lines(alld[[i]],col=modelcol[i],lwd=3))

}
legend("topright",legend=modelsubnames,col=modelcol,lwd=3,lty=1)
dev.off()
sapply(allexpe,function(i)quantile(i[[allscores[4]]],prob=.5)) 


##take 1000 best simulation
n=1000
listallposterior=lapply(allexpe,function(e)
                        {
                            list(
                                 best=e[order(e$distances),][1:n,],
                                 best=e[order(e$distances),][1:n,],
                                 worst=e[order(e$distances,decreasing=T),][1:n,],
                                 time_max=e[order(e$time_max,decreasing=T),][1:n,],
                                 time_max150=e[order(e$time_max150,decreasing=T),][1:n,],
                                 time_max250=e[order(e$time_max250,decreasing=T),][1:n,],
                                 max_infect=e[order(e$max_infect,decreasing=F),][1:n,],
                                 max_infect150=e[order(e$max_infect150,decreasing=F),][1:n,],
                                 max_infect250=e[order(e$max_infect250,decreasing=F),][1:n,],
                                 wtime_max=e[order(e$time_max,decreasing=F),][1:n,],
                                 wtime_max150=e[order(e$time_max150,decreasing=F),][1:n,],
                                 wtime_max250=e[order(e$time_max250,decreasing=F),][1:n,],
                                 wmax_infect=e[order(e$max_infect,decreasing=T),][1:n,],
                                 wmax_infect150=e[order(e$max_infect150,decreasing=T),][1:n,],
                                 wmax_infect250=e[order(e$max_infect250,decreasing=T),][1:n,]
                                 )
                        }
)

### FIGURE 5
for(mname in names(listallposterior)){
    png(paste0("all_switchfunctions_",mname,".png"),pointsize=17, res = 100, width = 7, height = 13, units = "in")
    x=seq(0,1,.01)
    par(mfrow=c(length(listallposterior[[mname]]),2),oma=c(3,3,2,0))
    par(mar=rep(.5,4))
    for(mod in names(listallposterior[[mname]])){
        best=listallposterior[[mname]][[mod]]
        plot(x,sig(x),type="n",ylim=c(0,1),xlim=c(0,1),xlab="% same age infected people", ylab="proba. to switch",main="",xaxt="n",yaxt="n")
        modsubname=paste0(strsplit(mod,"")[[1]][-(1:8)],collapse="")
        axis(2)
        mtext(modsubname,2.5,2,cex=.9)
        a=axis(1,outer=T)
        for(i in 1:nrow(best)) lines(x,sig(x,b=best$inf[i],a=best$sat[i]),ylim=c(0,1),xlim=c(0,1),col=alpha(modelcol[mod],.1),lwd=2) 
        plot(x,sig(x),type="n",ylim=c(0,1),xlim=c(0,1),xlab="% same age infected people", ylab="proba. to switch",main="",xaxt="n",yaxt="n")
        for(i in 1:nrow(best)) lines(x,sig(x,b=best$inf_r[i],a=best$sat_r[i]),ylim=c(0,1),xlim=c(0,1),col=alpha(modelcol[mod],.1),lwd=2) 
    }
    axis(1,outer=F)
    mtext(mname,3,1,outer=T)
    dev.off()
}


##FIGURE 7

modsubname=sapply(modelnames,function(mod)paste0(strsplit(mod,"")[[1]][-(1:8)],collapse=""))
png("",pointsize=17, res = 100, width = 7, height = 7, units = "in")
parameters=colnames(allexpe$midCurve15Good)[1:6]
par(mfrow=c(3,length(parameters)),mar=c(1,2,2,1))
for(best in listallposterior[6:8]){
for(var in parameters){
    p=unlist(sapply(allexpe,"[[",var))
    xlimits=range(sapply(best,"[[",var),p)
    alld=lapply(best,function(i,var)density(i[[var]],from=xlimits[1],to=xlimits[2]),var=var)
    pd=density(p)
    if(var %in% c("sat_r","sat")){
        print(var)
        xlimits=c(-1,3)
        alld=lapply(best,function(i,var)density(log10(i[[var]]),from=xlimits[1],to=xlimits[2]),var=var)
        pd=density(log10(p),from=xlimits[1],to=xlimits[2])
    }
    ylimits=range(c(sapply(alld,"[[","y"),pd$y))
    plot(1,1,ylim=ylimits,xlim=xlimits,type="n",main=var,yaxt="n",xlab=var,ylab="density")
    lines(pd,lwd=3)
    lapply(names(alld),function(i)lines(alld[[i]],col=modelcol[i],lwd=3))
if(var=="inf_r")legend("topleft",legend=modsubname,col=modelcol,lwd=3,lty=1)

}
}
dev.off()

#### better separation use those::
parameters=colnames(allexpe$midCurve15Good)[1:6]
for(n in names(listallposterior)){
    sbest=listallposterior[[n]]
    png(paste0("allposterior_metric_",n,".png"),width=1200,height=480,pointsize=15)
    par(mfrow=c(2,length(parameters)),mar=c(1,2,2,1),oma=c(3,3,0,0))
    sub=list(1:2,3:6)
    for( s in sub){
        best=sbest[s]
        for(var in parameters){
            p=unlist(sapply(allexpe,"[[",var))
            xlimits=range(sapply(best,"[[",var),p)
            alld=lapply(best,function(i,var)density(i[[var]],from=xlimits[1],to=xlimits[2]),var=var)
            pd=density(p)
            if(var %in% c("sat_r","sat")){
                print(var)
                xlimits=c(-1,3)
                alld=lapply(best,function(i,var)density(log10(i[[var]]),from=xlimits[1],to=xlimits[2]),var=var)
                pd=density(log10(p),from=xlimits[1],to=xlimits[2])
            }
            ylimits=range(c(sapply(alld,"[[","y"),pd$y))
            plot(1,1,ylim=ylimits,xlim=xlimits,type="n",main=var,yaxt="n",xlab=var,ylab="density")
            lines(pd,lwd=3)
            lapply(names(alld),function(i)lines(alld[[i]],col=modelcol[i],lwd=3))
            if(var=="sl_rad")legend("bottom",legend=modelsubnames[s],col=modelcol[s],lwd=3,lty=1)
        }

        }
        mtext("t0 : 15-100% A",2,1,0.3,outer=T) 
        mtext("t0 : 0% NA",2,1,0.7,outer=T) 
        dev.off()
}









#### FIGURE 8
for(n in names(listallposterior)){
    sbest=listallposterior[[n]]
    sbest=sbest[-5]


    png(paste0("posterior2d_",n,"_inf_r_sat_r.png"),pointsize=15, res = 100, width = 6, height = 6, units = "in")
    par(mar=c(1,1,1,1),oma=c(3,3,3,3),xpd=T)
    par(mfrow=c(3,2))
    for(mod in names(sbest)){
        best=sbest[[mod]]
        if(mod=="burnin100BestSLS")plot.new()
        hdr.boxplot.2d(best$inf_r,log10(best$sat_r),prob=seq(20,100,10),shadecols=alpha(modelcol[mod],1),xlim=c(0,1),ylim=c(-1,3),xlab="revert inf. point",ylab="revert steepness (log10)")
    }
    mtext("inf_r",1,2,outer=T)
    mtext("sat_r",2,2,outer=T)

    mtext("random SL",3,1,.3,outer=T)
    mtext("best SL",3,1,.7,outer=T)

    mtext("wait100A",4,1,.15,outer=T)
    mtext("15%A",4,1,.5,outer=T)
    mtext("0%AD",4,1,,.85,outer=T)
    dev.off()

    png(paste0("posterior2d_",n,"_inf_sat.png"),pointsize=15, res = 100, width = 6, height = 6, units = "in")
    par(mar=c(1,1,1,1),oma=c(3,3,3,3),xpd=T)
    par(mfrow=c(3,2))
    for(mod in names(sbest)){
        best=listallposterior$allbest[[mod]]
        if(mod=="burnin100BestSLS")plot.new()
        hdr.boxplot.2d(best$inf,log10(best$sat),prob=seq(20,100,10),shadecols=alpha(modelcol[mod],.9),xlim=c(0,1),ylim=c(-1,3),xlab="inflexion point",ylab="steepness (log10)")
    }
    mtext("inf",1,2,outer=T)
    mtext("sat",2,2,outer=T)
                                     
    mtext("random SL",3,1,.3,outer=T)
    mtext("best SL",3,1,.7,outer=T)

    mtext("wait100A",4,1,.15,outer=T)
    mtext("15%A",4,1,.5,outer=T)
    mtext("0%AD",4,1,,.85,outer=T)
                                 
    dev.off()

    png(paste0("posterior2d_",n,"_inf_pind.png"),pointsize=15, res = 100, width = 6, height = 6, units = "in")
    par(mar=c(1,1,1,1),oma=c(3,3,3,3),xpd=T)
    par(mfrow=c(3,2))
    for(mod in names(sbest)){
        best=listallposterior$allbest[[mod]]
        if(mod=="burnin100BestSLS")plot.new()
        hdr.boxplot.2d(best$inf,best$pind,prob=seq(20,100,10),shadecols=alpha(modelcol[mod],.9),xlim=c(0,1),ylim=c(0,1),xlab="inflexion point",ylab="individual learning")
    }
    mtext("inf",1,2,outer=T)
    mtext("pind",2,2,outer=T)

    mtext("random SL",3,1,.3,outer=T)
    mtext("best SL",3,1,.7,outer=T)

    mtext("wait100A",4,1,.15,outer=T)
    mtext("15%A",4,1,.5,outer=T)
    mtext("0%AD",4,1,,.85,outer=T)
    dev.off()

    png(paste0("posterior2d_",n,"_sat_pind.png"),pointsize=15, res = 100, width = 6, height = 6, units = "in")

    par(mar=c(1,1,1,1),oma=c(3,3,3,3),xpd=T)
    par(mfrow=c(3,2))
    for(mod in names(sbest)){
        best=listallposterior$allbest[[mod]]
        if(mod=="burnin100BestSLS")plot.new()
        hdr.boxplot.2d(log10(best$sat),best$pind,prob=seq(20,100,10),shadecols=alpha(modelcol[mod],.9),xlim=c(-1,3),ylim=c(0,1),xlab="steepness (log10)",ylab="individual learning")
    }
    mtext("sat",1,2,outer=T)
    mtext("pind",2,2,outer=T)

    mtext("random SL",3,1,.3,outer=T)
    mtext("best SL",3,1,.7,outer=T)
                                
    mtext("wait100A",4,1,.15,outer=T)
    mtext("15%A",4,1,.5,outer=T)
    mtext("0%AD",4,1,,.85,outer=T)
                                 
    dev.off()

    png(paste0("posterior2d_",n,"_sl_rad_pind.png"),pointsize=15, res = 100, width = 6, height = 6, units = "in")
    par(mar=c(1,1,1,1),oma=c(3,3,3,3),xpd=T)
    par(mfrow=c(3,2))
    for(mod in names(sbest)){
        best=listallposterior$allbest[[mod]]
        if(mod=="burnin100BestSLS")plot.new()
        hdr.boxplot.2d(best$sl_rad,best$pind,prob=seq(20,100,10),shadecols=alpha(modelcol[mod],.9),xlim=c(1,100),ylim=c(0,1),xlab="radius social learning",ylab="individual learning")
    }
    mtext("sl_rad",1,2,outer=T)
    mtext("pind",2,2,outer=T)
    mtext("random SL",3,1,.3,outer=T)
    mtext("best SL",3,1,.7,outer=T)
    
    mtext("wait100A",4,1,.15,outer=T)
    mtext("15%A",4,1,.5,outer=T)
    mtext("0%AD",4,1,,.85,outer=T)
    dev.off()

    png(paste0("posterior2d_",n,"_sl_rad_inf.png"),pointsize=15, res = 100, width = 6, height = 6, units = "in")
    par(mar=c(1,1,1,1),oma=c(3,3,3,3),xpd=T)
    par(mfrow=c(3,2))
    for(mod in names(sbest)){
        best=listallposterior$allbest[[mod]]
        if(mod=="burnin100BestSLS")plot.new()
        hdr.boxplot.2d(best$sl_rad,best$inf,prob=seq(20,100,10),shadecols=alpha(modelcol[mod],.9),xlim=c(1,100),ylim=c(0,1),xlab="radius social learning",ylab="infection points")
    }
    mtext("sl_rad",1,2,outer=T)
    mtext("inf",2,2,outer=T)

    mtext("random SL",3,1,.3,outer=T)
    mtext("best SL",3,1,.7,outer=T)
                                
    mtext("wait100A",4,1,.15,outer=T)
    mtext("15%A",4,1,.5,outer=T)
    mtext("0%AD",4,1,,.85,outer=T)
                                 
    dev.off()

    png(paste0("posterior2d_",n,"_sl_rad_sat.png"),pointsize=15, res = 100, width = 6, height = 6, units = "in")
par(xpd=NA)
    par(mar=c(1,1,1,1),oma=c(3,3,3,3),xpd=T)
    par(mfrow=c(3,2))
    for(mod in names(sbest)){
        best=listallposterior$allbest[[mod]]
        if(mod=="burnin100BestSLS")plot.new()
        hdr.boxplot.2d(best$sl_rad,log10(best$sat),prob=seq(20,100,10),shadecols=alpha(modelcol[mod],.9),xlim=c(1,100),ylim=c(-1,3),xlab="sl_rad",ylab="steepness (log10)")
    }
    mtext("sl_rad",1,2,outer=T)
    mtext("sat",2,2,outer=T)

    mtext("random SL",3,1,.3,outer=T)
    mtext("best SL",3,1,.7,outer=T)
                                
    mtext("wait100A",4,1,.15,outer=T)
    mtext("15%A",4,1,.5,outer=T)
    mtext("0%AD",4,1,,.85,outer=T)
                                 
    dev.off()
}



par(mfrow=c(2,4),oma=c(0,4,0,0))
 
plot(density(log10(listallposterior$allbest$midCurveAllBad$sat),from=-1,to=3),ylim=c(0,.45),main="sat",xlab="sat")
lines(density(log10(listallposterior$allbest_max_infect150$midCurveAllBad$sat)),col="red")
lines(density(log10(allexpe$midCurveAllBad$sat)),col="grey")
legend("topright",legend=c("final time","150ts", "prior"),lwd=1,col=c(1,2,"grey"))
mtext("normal model",2,5)

plot(density(listallposterior$allbest$midCurveAllBad$inf,from=0,to=1),ylim=c(0,2),main="inf",xlab="inf")
lines(density(listallposterior$allbest_max_infect150$midCurveAllBad$inf,from=0,to=1),col="red")
lines(density(allexpe$midCurveAllBad$inf),col="grey")

plot(density(log10(listallposterior$allbest$midCurveAllBad$sat_r),from=-1,to=3),ylim=c(0,.7),main="sat_r",xlab="sat_r")
lines(density(log10(listallposterior$allbest_max_infect150$midCurveAllBad$sat_r)),col="red")
lines(density(log10(allexpe$midCurveAllBad$sat_r)),col="grey")

plot(density(listallposterior$allbest$midCurveAllBad$inf_r,from=0,to=1),ylim=c(0,6),main="inf_r",xlab="inf_r")
lines(density(listallposterior$allbest_max_infect150$midCurveAllBad$inf_r,from=0,to=1),col="red")
lines(density(allexpe$midCurveAllBad$inf_r),col="grey")

plot(density(log10(listallposterior$allbest$midCurveAllBad$sat),from=-1,to=3),ylim=c(0,.45),main="sat",xlab="sat")
lines(density(log10(listallposterior$allbest_max_infect150$midCurveAllBad$sat)),col="red")
lines(density(log10(allexpe$midCurveAllBad$sat)),col="grey")
legend("topright",legend=c("final time","150ts", "prior"),lwd=1,col=c(1,2,"grey"))
mtext("waitn 100ts ",2,5)

plot(density(listallposterior$allbest$midCurveAllBad$inf,from=0,to=1),ylim=c(0,2),main="inf",xlab="inf")
lines(density(listallposterior$allbest_max_infect150$midCurveAllBad$inf,from=0,to=1),col="red")
lines(density(allexpe$midCurveAllBad$inf),col="grey")

plot(density(log10(listallposterior$allbest$midCurveAllBad$sat_r),from=-1,to=3),ylim=c(0,.7),main="sat_r",xlab="sat_r")
lines(density(log10(listallposterior$allbest_max_infect150$midCurveAllBad$sat_r)),col="red")
lines(density(log10(allexpe$midCurveAllBad$sat_r)),col="grey")

plot(density(listallposterior$allbest$midCurveAllBad$inf_r,from=0,to=1),ylim=c(0,6),main="inf_r",xlab="inf_r")
lines(density(listallposterior$allbest_max_infect150$midCurveAllBad$inf_r,from=0,to=1),col="red")
lines(density(allexpe$midCurveAllBad$inf_r),col="grey")



par(mfrow=c(1,2))
plot(1,1,ylim=c(0,500),xlim=c(0,500),type="n",main="#infected")
 for(f in list.files("testingSimHigSL/" ,pattern="allG_*",full.names=T)){ 
     load(f)
     lines(simu$timeseries[,2],col="green")
}
 for(f in list.files("testingSimHigSL/" ,pattern="allB_*",full.names=T)){ 
     load(f)
     lines(simu$timeseries[,2],col="red")
}
plot(1,1,ylim=c(0,500),xlim=c(0,500),type="n",main="#Adherent")
 for(f in list.files("testingSimHigSL/" ,pattern="allG_*",full.names=T)){ 
     load(f)
     lines(simu$timeseries[,5],col="green")
}
 for(f in list.files("testingSimHigSL/" ,pattern="allB_*",full.names=T)){ 
     load(f)
     lines(simu$timeseries[,5],col="red")
}

duocolrs=alpha(c(myred,colorbest),.8)
currange=range=list(dim1=c(0,1),dim2=c(-1,3))
dimlab=list(dim1=expression(nu),dim2=expression(kappa))

### Cmpound marginal with 2d  learning

for(i in names(allexpe)){
    posterior=listallposterior[[i]]
    worst=allexpe[[i]]

    pdf(paste0("5e8779c56517890001536101/figures/2d_marginal_learning_best_vs_worst_",i,".pdf"),pointsize=30,  width = 10, height = 10)
    expnames=c(distribA="Worst",distribB="Best")
    marginAndJoin(
                  distribA=list(dim1=posterior$worst$inf,dim2=posterior$worst$sat),
                  distribB=list(dim1=posterior$best$inf,dim2=posterior$best$sat),
                  dimlab=list(dim1=expression(nu),dim2=expression(kappa)),
                  expnames=expnames, range=currange,cols=duocolrs,log="y",
                  main=""
                  )
    dev.off()


    pdf(paste0("5e8779c56517890001536101/figures/2d_marginal_learning_best_vs_worst_150_",i,".pdf"),pointsize=30,  width = 10, height = 10)
    expnames=c(distribA="Worst 150",distribB="Best 150")
    marginAndJoin(
                  distribA=list(dim1=posterior$wmax_infect150$inf,dim2=posterior$wmax_infect150$sat),
                  distribB=list(dim1=posterior$max_infect150$inf,dim2=posterior$max_infect150$sat),
                  dimlab=list(dim1=expression(nu),dim2=expression(kappa)),
                  expnames=expnames, range=currange,cols=duocolrs,log="y",
                  main=""
                  )
    dev.off()


    ### COmpound marginal with 2d revert learning

    pdf(paste0("5e8779c56517890001536101/figures/2d_marginal_revert_learning_best_vs_worst_",i,".pdf"),pointsize=30,  width = 10, height = 10)
    expnames=c(distribA="Worst",distribB="Best")
    marginAndJoin(
                  distribA=list(dim1=posterior$worst$inf_r,dim2=posterior$worst$sat_r),
                  distribB=list(dim1=posterior$best$inf_r,dim2=posterior$best$sat_r),
                  dimlab=list(dim1=expression(nu[r]),dim2=expression(kappa[r])),
                  expnames=expnames, range=currange,cols=duocolrs,log="y",
                  main=""
                  )
    dev.off()


    pdf(paste0("5e8779c56517890001536101/figures/2d_marginal_revert_learning_best_vs_worst_150_",i,".pdf"),pointsize=30,  width = 10, height = 10)
    expnames=c(distribA="Worst 150",distribB="Best 150")
    marginAndJoin(
                  distribA=list(dim1=posterior$wmax_infect150$inf_r,dim2=posterior$wmax_infect150$sat_r),
                  distribB=list(dim1=posterior$max_infect150$inf_r,dim2=posterior$max_infect150$sat_r),
                  dimlab=list(dim1=expression(nu[r]),dim2=expression(kappa[r])),
                  expnames=expnames, range=currange,cols=duocolrs,log="y",
                  main=""
                  )
    dev.off()


    ################ best vs prior

    ### Cmpound marginal with 2d  learning

    pdf(paste0("5e8779c56517890001536101/figures/2d_marginal_learning_best_vs_prior_",i,".pdf"),pointsize=30,  width = 10, height = 10)
    expnames=c(distribA="Prior",distribB="Best")
    marginAndJoin(
                  distribA=list(dim1=worst$inf,dim2=worst$sat),
                  distribB=list(dim1=posterior$best$inf,dim2=posterior$best$sat),
                  dimlab=list(dim1=expression(nu),dim2=expression(kappa)),
                  expnames=expnames, range=currange,cols=duocolrs,log="y",
                  main=""
                  )
    dev.off()


    pdf(paste0("5e8779c56517890001536101/figures/2d_marginal_learning_best_vs_prior_150_",i,".pdf"),pointsize=30,  width = 10, height = 10)
    expnames=c(distribA="Prior",distribB="Best 150")
    marginAndJoin(
                  distribA=list(dim1=worst$inf,dim2=worst$sat),
                  distribB=list(dim1=posterior$max_infect150$inf,dim2=posterior$max_infect150$sat),
                  dimlab=list(dim1=expression(nu),dim2=expression(kappa)),
                  expnames=expnames, range=currange,cols=duocolrs,log="y",
                  main=""
                  )
    dev.off()


    ### COmpound marginal with 2d revert learning

    pdf(paste0("5e8779c56517890001536101/figures/2d_marginal_revert_learning_best_vs_prior_",i,".pdf"),pointsize=30,  width = 10, height = 10)
    expnames=c(distribA="Prior",distribB="Best")
    marginAndJoin(
                  distribA=list(dim1=worst$inf_r,dim2=worst$sat_r),
                  distribB=list(dim1=posterior$best$inf_r,dim2=posterior$best$sat_r),
                  dimlab=list(dim1=expression(nu[r]),dim2=expression(kappa[r])),
                  expnames=expnames, range=currange,cols=duocolrs,log="y",
                  main=""
                  )
    dev.off()


    pdf(paste0("5e8779c56517890001536101/figures/2d_marginal_revert_learning_best_vs_prior_150_",i,".pdf"),pointsize=30,  width = 10, height = 10, fillOddEven = F)
    expnames=c(distribA="Prior",distribB="Best 150")
    marginAndJoin(
                  distribA=list(dim1=worst$inf_r,dim2=worst$sat_r),
                  distribB=list(dim1=posterior$max_infect150$inf_r,dim2=posterior$max_infect150$sat_r),
                  dimlab=list(dim1=expression(nu[r]),dim2=expression(kappa[r])),
                  expnames=expnames, range=currange,cols=duocolrs,log="y",
                  main=""
                  )
    dev.off()
}

##FIGURE 7
for(i in names(allexpe)){
    posterior=listallposterior[[i]]
mar=c(4,3,2,1)

colors=c(P="NA",A=alpha(colorbest,1),B=alpha(colorbest,1))

pdf(paste0("5e8779c56517890001536101/figures/posterior_pind_",i,".pdf"),pointsize=30,  width = 10, height = 10)
par(mar=mar)
 plot2dens(prior=1-allresults$pind,A=1-posterior$best$pind,B=1-posterior$max_infect250$pind,cols=colors,from=0,to=1,main="posterior proba social learning",xlab="proba social learning")
dev.off()

pdf(paste0("5e8779c56517890001536101/figures/posterior_sl_rad_",i,".pdf"),pointsize=30,  width = 10, height = 10)
par(mar=mar)
 plot2dens(prior=allresults$sl_rad,A=posterior$best$sl_rad,B=posterior$max_infect150$sl_rad,cols=colors,from=1,to=100,main="",xlab="radius social learning")
dev.off()

pdf(paste0("5e8779c56517890001536101/figures/posterior_sat_r_",i,".pdf"),pointsize=30,  width = 10, height = 10)
par(mar=mar)
 plot2dens(prior=log10(allresults$sat_r),A=log10(posterior$best$sat_r),B=log10(posterior$time_max150$sat_r),cols=colors,from=-1,to=3,main="",xlab=expression(kappa[r]),log=T)
dev.off()

pdf(paste0("5e8779c56517890001536101/figures/posterior_sat_",i,".pdf"),pointsize=30,  width = 10, height = 10)
par(mar=mar)
 plot2dens(prior=log10(allresults$sat),A=log10(posterior$best$sat),B=log10(posterior$max_infect150$sat),cols=colors,from=-1,to=3,xlab=expression(kappa),main="",xaxt="n",log=T)
dev.off()

pdf(paste0("5e8779c56517890001536101/figures/posterior_inf_r_",i,".pdf"),pointsize=30,  width = 10, height = 10)
par(mar=mar)
 plot2dens(prior=(allresults$inf_r),A=(posterior$best$inf_r),B=(posterior$max_infect150$inf_r),cols=colors,from=0,to=1,main="posterior revert inflection point",xlab=expression(nu[r]),xaxt="n")
legend("topleft",legend=c("prior","final time","150 timesteps"),fill=c(0,colorbest,colorbest),density=c(NA,NA,20))
dev.off()

pdf(paste0("5e8779c56517890001536101/figures/posterior_inf_",i,".pdf"),pointsize=30,  width = 10, height = 10)
par(mar=mar)
 plot2dens(prior=(allresults$inf),A=(posterior$best$inf),B=(posterior$max_infect150$inf),cols=colors,from=0,to=1,main="posterior inflection point",xlab=expression(nu),xaxt="n")
dev.off()
}
