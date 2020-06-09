
source("abmEpi.R")

poptest=generatePopulation(500,recovery=2500) 
neutral=replicate(50,abmSIR(poptest,2500,speed=c(1,.2),p=c(1,1),di=2,i0=1,visu=F)$timeseries[,2])
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
    neutral=apply(replicate(50,abmSIR(poptest,1500,speed=.6,p=c(.2,.2),di=2,i0=1,visu=F)$timeseries[,2]),1,mean)
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
    tstep=1000
    par(mfrow=c(2,3))
    par(mar=c(4,4,1,1))
    for(i in 1:length(inpoint)){
        plot(1:tstep,neutral,type="l",xlim=c(1,tstep),ylim=c(0,500),main=paste("inpoint=",inpoint[i]),ylab="#infected",xlab="",lwd=2)
        for(u in 1:length(probas)){
            lines(1:tstep,means[[u]][[i]],lty=1,col=clrsprobs[u])
        }
        legend("topright",legend=c("neutral",paste0("PiB=",1/probas[c(1,2,4,10)],"PiG")),lty=c(1,rep(1,4)),col=c(1,clrsprobs[c(1,2,4,10)]),lwd=2)
    }
    par(mar=c(4,4,1,1))
    for(i in (1:length(probas))[c(1,5,9)]){
        plot(neutral,type="l",xlim=c(1,tstep),ylim=c(0,500),main=bquote(P*i[B]==.(1/probas[i])*P*i[G]),xlab="#time",ylab="#infected",lwd=2)
        for(u in 1:length(inpoint)){
            lines(1:tstep,means[[i]][[u]],lty=1,col=clrsinp[u])
        }
        legend("topright",legend=c("neutral",paste0("inpoininpoint=",inpoint)),lty=c(1,rep(1,3)),col=c(1,clrsinp),lwd=2)
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
	poptest[, "behavior"]=B
    a=abmSIR(poptest,1500,p=c(1,.2),di=2,i0=1,inf=.9,sat=5,inf_r=.9,sat_r=5,xsize=xsize,ysize=ysize,visu=T,ap=F,ts=T)
    #neutral=lapply(1:100,function(j){print(j);abmSIR(poptest,1000,p=c(1,1),di=2,i0=1,visu=F,inf=1,sat=20,xsize=xsize,ysize=ysize)}$timeseries[,2])
    #baseline=mean(lapply(neutral,max))
    timeA=mean(sapply(neutral,sapply,getTimeMaxTotal))
    timeB=mean(sapply(neutral,sapply,getTimeMaxInfected))
    par(mfrow=c(1,2))
    image(t(sapply(allres,sapply,function(e)timeA-mean(sapply(e,sapply,getTimeMaxTotal)))))
    image(t(sapply(allres,sapply,function(e)timeB-mean(sapply(e,sapply,getTimeMaxTotal)))))

    dev.new()
    load("allres.bin")
    infected=lapply(allres,lapply,lapply,function(i)i$timeseries[,2])
    table_i=lapply(infected,lapply,function(i)do.call("cbind",i))
    mean_i=lapply(table_i,lapply,function(i)apply(i,1,mean))

    load("oldbin/neutral.bin")
    neutral=apply(do.call("cbind",lapply(neutral,function(u)u$timeseries[,2])),1,mean)
    neutralVar=apply(do.call("cbind",lapply(neutral,function(u)u$timeseries[,2])),1,var)
    load("neutralBest")
    neutralBest=apply(do.call("cbind",lapply(neutralBest,function(u)u$timeseries[,2])),1,mean)
    neutralBestVar=apply(do.call("cbind",lapply(neutralBest,function(u)u$timeseries[,2])),1,var)

    pdf("twofirstdim.pdf",width=14,height=5)
    tstep=1000
    par(mfrow=c(2,3))
    par(mar=c(4,4,1,1))
    for(i in 1:length(inpoint)){
        plot(1:tstep,neutralBest,type="l",xlim=c(1,tstep),ylim=c(0,500),main=paste("inpoint=",inpoint[i]),ylab="#infected",xlab="",lwd=2)
        for(u in 1:length(probas)){
            lines(1:tstep,mean_i[[u]][[i]],lty=1,col=clrsprobs[u])
        }
        legend("topright",legend=c("neutral",paste0("PiB=",1/probas[c(1,2,4,10)],"PiG")),lty=c(1,rep(1,4)),col=c(1,clrsprobs[c(1,2,4,10)]),lwd=2)
    }
    par(mar=c(4,4,1,1))
    for(i in (1:length(probas))[c(1,5,9)]){
        plot(neutral,type="l",xlim=c(1,tstep),ylim=c(0,500),main=bquote(P*i[B]==.(1/probas[i])*P*i[G]),xlab="#time",ylab="#infected",lwd=2)
        for(u in 1:length(inpoint)){
            lines(1:tstep,mean_i[[i]][[u]],lty=1,col=clrsinp[u])
        }
        legend("topright",legend=c("neutral",paste0("inpoininpoint=",inpoint)),lty=c(1,rep(1,3)),col=c(1,clrsinp),lwd=2)
    }
    dev.off()

    inpoint=sort(1/2^seq(1,10))
    probas=c(0.1,.4,.8)
    clrsprobs=colorRampPalette(c("dark green","yellow"))(length(probas))
    clrsinp=colorRampPalette(c("blue","red"))(length(inpoint))


    pdf("revert.pdf",width=14,height=5)
    load("revert.bin")
    infected=lapply(revert,lapply,lapply,function(i)i$timeseries[,2])
    table_i=lapply(infected,lapply,function(i)do.call("cbind",i))
    mean_i=lapply(table_i,lapply,function(i)apply(i,1,mean))
    tstep=1000
    par(mfrow=c(2,3))
    par(mar=c(4,4,1,1))
    for(i in (1:length(inpoint))[c(1,5,9)]){
        plot(1:tstep,neutral,type="l",xlim=c(1,tstep),ylim=c(0,500),main=paste("inpoint=",inpoint[i]),ylab="#infected",xlab="",lwd=2)
        for(u in 1:length(probas)){
            lines(1:tstep,mean_i[[u]][[i]],lty=1,col=clrsprobs[u])
        }
        legend("topright",legend=c("neutral",paste0("PiB=",1/probas,"PiG")),lty=c(1,rep(1,4)),col=c(1,clrsprobs),lwd=2)
    }
    par(mar=c(4,4,1,1))
    for(i in (1:length(probas))){
        plot(neutral,type="l",xlim=c(1,tstep),ylim=c(0,500),main=bquote(P*i[B]==.(1/probas[i])*P*i[G]),xlab="#time",ylab="#infected",lwd=2)
        for(u in 1:length(inpoint)){
            lines(1:tstep,mean_i[[i]][[u]],lty=1,col=clrsinp[u])
        }
        legend("topright",legend=c("neutral",paste0("inpoininpoint=",inpoint[c(1,4,8,10)])),lty=c(1,rep(1,4)),col=c(1,clrsinp[c(1,4,8,10)]),lwd=2)
    }
    dev.off()

    pdf("norevert.pdf",width=14,height=5)
    load("norevert.bin")
    infected=lapply(norevert,lapply,lapply,function(i)i$timeseries[,2])
    table_i=lapply(infected,lapply,function(i)do.call("cbind",i))
    mean_i=lapply(table_i,lapply,function(i)apply(i,1,mean))
    tstep=1000
    par(mfrow=c(2,3))
    par(mar=c(4,4,1,1))
    for(i in (1:length(inpoint))[c(1,5,9)]){
        plot(1:tstep,neutral,type="l",xlim=c(1,tstep),ylim=c(0,500),main=paste("inpoint=",inpoint[i]),ylab="#infected",xlab="",lwd=2)
        for(u in 1:length(probas)){
            lines(1:tstep,mean_i[[u]][[i]],lty=1,col=clrsprobs[u])
        }
        legend("topright",legend=c("neutral",paste0("PiB=",1/probas,"PiG")),lty=c(1,rep(1,4)),col=c(1,clrsprobs),lwd=2)
    }
    par(mar=c(4,4,1,1))
    for(i in (1:length(probas))){
        plot(neutral,type="l",xlim=c(1,tstep),ylim=c(0,500),main=bquote(P*i[B]==.(1/probas[i])*P*i[G]),xlab="#time",ylab="#infected",lwd=2)
        for(u in 1:length(inpoint)){
            lines(1:tstep,mean_i[[i]][[u]],lty=1,col=clrsinp[u])
        }
        legend("topright",legend=c("neutral",paste0("inpoininpoint=",inpoint[c(1,4,8,10)])),lty=c(1,rep(1,4)),col=c(1,clrsinp[c(1,4,8,10)]),lwd=2)
    }
    dev.off()

    load("revertPsoc.bin")
    infected=lapply(revertPsoc,lapply,lapply,function(i)i$timeseries[,2])
    table_i=lapply(infected,lapply,function(i)do.call("cbind",i))
    mean_i=lapply(table_i,lapply,function(i)apply(i,1,mean))
    tstep=1000
    par(mfrow=c(2,3))
    par(mar=c(4,4,1,1))
    for(i in (1:length(inpoint))[c(1,5,9)]){
        plot(1:tstep,neutral,type="l",xlim=c(1,tstep),ylim=c(0,500),main=paste("inpoint=",inpoint[i]),ylab="#infected",xlab="",lwd=2)
        for(u in 1:length(probas)){
            lines(1:tstep,mean_i[[u]][[i]],lty=1,col=clrsprobs[u])
        }
    }
    par(mar=c(4,4,1,1))
    for(i in (1:length(probas))){
        plot(neutral,type="l",xlim=c(1,tstep),ylim=c(0,500),main=bquote(P*i[B]==.(1/probas[i])*P*i[G]),xlab="#time",ylab="#infected",lwd=2)
        for(u in 1:length(inpoint)){
            lines(1:tstep,mean_i[[i]][[u]],lty=1,col=clrsinp[u])
        }
        legend("topright",legend=c("neutral",paste0("inpoininpoint=",inpoint[c(1,4,8,10)])),lty=c(1,rep(1,4)),col=c(1,clrsinp[c(1,4,8,10)]),lwd=2)
    }


    load("norevert.bin")
    infected=lapply(norevert,lapply,lapply,function(i)i$timeseries[,2])
    table_i=lapply(infected,lapply,function(i)do.call("cbind",i))
    mean_i=lapply(table_i,lapply,function(i)apply(i,1,mean))
    tstep=2500
    par(mfrow=c(2,3))
    par(mar=c(4,4,1,1))
    for(i in (1:length(inpoint))[c(1,5,9)]){
        plot(neutral,type="l",xlim=c(1,tstep),ylim=c(0,500),main=paste("inpoint=",inpoint[i]),ylab="#infected",xlab="",lwd=2)
        for(u in 1:length(probas)){
            lines(1:tstep,mean_i[[u]][[i]],lty=1,col=clrsprobs[u])
        }
        legend("topright",legend=c("neutral",paste0("PiB=",1/probas,"PiG")),lty=c(1,rep(1,4)),col=c(1,clrsprobs),lwd=2)
    }
    par(mar=c(4,4,1,1))
    for(i in (1:length(probas))){
        plot(neutral,type="l",xlim=c(1,tstep),ylim=c(0,500),main=bquote(P*i[B]==.(1/probas[i])*P*i[G]),xlab="#time",ylab="#infected",lwd=2)
        for(u in 1:length(inpoint)){
            lines(1:tstep,mean_i[[i]][[u]],lty=1,col=clrsinp[u])
        }
        legend("topright",legend=c("neutral",paste0("inpoininpoint=",inpoint[c(1,4,8,10)])),lty=c(1,rep(1,4)),col=c(1,clrsinp[c(1,4,8,10)]),lwd=2)
    }

    allExpe=c("norevert","revert","revertLessSteep","norevertLessSteep" )
    names(allExpe)=allExpe
    allMeans=lapply( allExpe,function(i)
                    {
                        load(paste0("oldbin/",i,".bin"))
                        infected=lapply(get(i),lapply,lapply,lapply,function(r)r$timeseries[,2])
                        table_i=lapply(infected,lapply,lapply,function(r)do.call("cbind",r))
                        mean_i=lapply(table_i,lapply,lapply,function(r)apply(r,1,mean))
                    }
    )

    allVar=lapply( allExpe,function(i)
                    {
                        load(paste0("oldbin/",i,".bin"))
                        infected=lapply(get(i),lapply,lapply,lapply,function(r)r$timeseries[,2])
                        table_i=lapply(infected,lapply,lapply,function(r)do.call("cbind",r))
                        var_i=lapply(table_i,lapply,lapply,function(r)apply(r,1,var))
                    }
    )
    

    pdf("steep.pdf")
    inpoint=c(0.01,.5)
    probas=c(0.1,.4)
    psoc=c(0)
    tstep=2500
        clrs=colorRampPalette(c("blue","red"))(length(psoc))
        par(mfrow=c(2,2))
        for(u in c(1,2)){
        for( i  in c(1,2)){
        plot(neutral,type="l",xlim=c(1,1500),ylim=c(0,500),main=bquote(P*i[B]==.(1/probas[i])*P*i[G] ~ inp == .(inpoint[u]) ~ stp == 20),xlab="#time",ylab="#infected",lwd=2)
        for(p_soc in 1:length(psoc)){
        lines(1:tstep,allMeans[["revert"]][[p_soc]][[i]][[u]],lty=1,col=clrs[p_soc],lwd=2)
        lines(1:tstep,allMeans[["norevert"]][[p_soc]][[i]][[u]],lty=2,col=clrs[p_soc],lwd=2)
        }
        legend("topright",legend=c("neutral","revert","norevert",paste0("p_social=",psoc)),lty=c(1,1,2,rep(1,3)),col=c(1,1,1,clrs),lwd=2)
        }
        }
        dev.off()

    pdf("lessSteep.pdf")
    inpoint=c(0.01,.5)
    probas=c(0.1,.4)
    psoc=c(0)
    tstep=2500
        clrs=colorRampPalette(c("blue","red"))(length(psoc))
        par(mfrow=c(2,2))
        for(u in c(1,2)){
        for( i  in c(1,2)){
        plot(neutral,type="l",xlim=c(1,1500),ylim=c(0,500),main=bquote(P*i[B]==.(1/probas[i])*P*i[G] ~ inp == .(inpoint[u]) ~ stp == 5),xlab="#time",ylab="#infected",lwd=2)
        for(p_soc in 1:length(psoc)){
        lines(1:tstep,allMeans[["revertLessSteep"]][[p_soc]][[i]][[u]],lty=1,col=clrs[p_soc],lwd=2)
        lines(1:tstep,allMeans[["norevertLessSteep"]][[p_soc]][[i]][[u]],lty=2,col=clrs[p_soc],lwd=2)
        }
        legend("topright",legend=c("neutral","revert","norevert",paste0("p_social=",psoc)),lty=c(1,1,2,rep(1,3)),col=c(1,1,1,clrs),lwd=2)
        }
        }
        dev.off()




    pdf("lessSteepVar.pdf")
    inpoint=c(0.01,.5)
    probas=c(0.1,.4)
    psoc=c(0,.1,.5)
    tstep=2500
        clrs=colorRampPalette(c("blue","red"))(length(psoc))
        par(mfrow=c(2,2))
        for(u in c(1,2)){
        for( i  in c(1,2)){
        plot(neutralBestVar,type="l",xlim=c(1,1500),ylim=range(allVar),main=bquote(P*i[B]==.(1/probas[i])*P*i[G] ~ inp == .(inpoint[u]) ~ stp == 5),xlab="#time",ylab="#infected(var)",lwd=2)
        for(p_soc in 1:length(psoc)){
        lines(1:tstep,allVar[["revertLessSteep"]][[p_soc]][[i]][[u]],lty=1,col=clrs[p_soc],lwd=2)
        lines(1:tstep,allVar[["norevertLessSteep"]][[p_soc]][[i]][[u]],lty=2,col=clrs[p_soc],lwd=2)
        }
        legend("topright",legend=c("neutral","revert","norevert",paste0("p_social=",psoc)),lty=c(1,1,2,rep(1,3)),col=c(1,1,1,clrs),lwd=2)
        }
        }
        dev.off()


    pdf("steepVar.pdf")
    inpoint=c(0.01,.5)
    probas=c(0.1,.4)
    psoc=c(0,.1,.5)
    tstep=2500
        clrs=colorRampPalette(c("blue","red"))(length(psoc))
        par(mfrow=c(2,2))
        for(u in c(1,2)){
        for( i  in c(1,2)){
        plot(neutral,type="l",xlim=c(1,1500),ylim=range(allVar),main=bquote(P*i[B]==.(1/probas[i])*P*i[G] ~ inp == .(inpoint[u]) ~ stp == 20),xlab="#time",ylab="#infected",lwd=2)
        for(p_soc in 1:length(psoc)){
        lines(1:tstep,allVar[["revert"]][[p_soc]][[i]][[u]],lty=1,col=clrs[p_soc],lwd=2)
        lines(1:tstep,allVar[["norevert"]][[p_soc]][[i]][[u]],lty=2,col=clrs[p_soc],lwd=2)
        }
        legend("topright",legend=c("neutral","revert","norevert",paste0("p_social=",psoc)),lty=c(1,1,2,rep(1,3)),col=c(1,1,1,clrs),lwd=2)
        }
        }
        dev.off()




}

getTimeMaxTotal <- function(pop)min(which(pop[,2]+pop[,3]==max(pop[,2]+pop[,3])))
getTimeMaxInfected <- function(pop)min(which(pop[,2]==max(pop[,2])))


printSigmoid  <- function(){


    pdf("5e8779c56517890001536101/figures/sigmoid.pdf",width=9,height=5)
    par(mfrow=c(1,2))
    x=seq(0,1,.01)   
    inp=seq(0,1,.1)   
    clrs=colorRampPalette(c("dark green","yellow"))(length(inp))
    plot(x,sig(x),type="n",ylim=c(0,1),xlim=c(0,1),ylab="P(B->G)~sig(x,st=10,inp)",main=expression(inp %in%  "(" * list(0,1) * ")") )
    for(i in 1:length(inp)) lines(x,sig(x,b=inp[i],a=4),ylim=c(0,1),xlim=c(0,1),col=clrs[i],lwd=2) 
    leg=seq(1,length(inp),length.out=3)
    legend("bottomright",legend=paste("inp=",inp[leg]),col=clrs[leg],lwd=2)

    tstp=20
    stp=seq(-3,3,length.out=tstp)   
    clrs=colorRampPalette(c("dark green","yellow"))(tstp)
    plot(x,sig(x),type="n",ylim=c(0,1),xlim=c(0,1),ylab="P(B->G)~sig(x,st,inp=.5)",main=bquote(stp %in% "(" * list(10^.(stp[1]),10^.(stp[tstp])) * ")") )
    for(i in 1:length(stp)) lines(x,sig(x,a=10^stp[i]),ylim=c(0,1),xlim=c(0,1),col=clrs[i],lwd=2) 
    leg=seq(1,tstp,length.out=4)
    legend("bottomright",legend=paste("stp=10^",round(stp[leg])),col=clrs[leg],lwd=2)
    dev.off()

}



    inf=seq(0,1,length.out=10)
    #inf_r=seq(0,1,length.out=10)
    sat=10^seq(-3,3,length.out=10)
    #sat_r=10^seq(-3,3,length.out=10)
    pind=c(.01,.5,.99)
    #allparameter=expand.grid(inf=inf,inf_r=inf_r,sat=sat,sat_r=sat_r,pind=pind)
    allparameter=expand.grid(inf=inf,sat=sat,pind=pind)
    xsize=ysize=100
    poptest=generatePopulation(500,recovery=c(8,14)*25,speed=c(1,.2),xsize=xsize,ysize=ysize,behavior=rep(G,500) )
    old <- Sys.time() 
    a=abmSIR(poptest,500,p=c(1,.2),di=2,i0=1,inf=.8,sat=10,xsize=xsize,ysize=ysize,visu=T,ap=F,ts=T,file=F,p_i=.5)
    print(Sys.time()-old )

    color.gradient <- function(x, colors=c("red","yellow","green"), colsteps=100) {
              return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
    }


    zscore <- function(x)(x-mean(x))/sd(x)
    bscore <- function(x)(x-min(x))/(max(x)-min(x))
    load("exnores.bin")
    load("allresults.bin")
    allresults=allresults[-which(allresults$max_infect <4),]
    allresults$scores=(1-bscore(allresults$time_max) + bscore(allresults$max_infect))/2

    par(mfrow=c(2,2))
    test=allresults[allresults$pind == .99,]
    plot(test$sat,test$max_infect,col=alpha(color.gradient(test$inf),.2),log="x",pch=20)
    plot(test$sat,test$time_max,col=alpha(color.gradient(test$inf),.2),log="x",pch=20)
    test=allresults


    png("5e8779c56517890001536101/figures/threeDimensionDistance.png",width=1000,height=400,pointsize=20)
    par(mfrow=c(1,3),cex=.5)
	test=allresults
    plot(test$time_max,test$max_infect,col=alpha(color.gradient(log(test$sat),c("blue","red")),.2),main="steepness",pch=20,xlab="Time Max",ylab="Max Infected")
    legend("topright",legend=paste("steepness=10^",round(log10(unique(test$sat)[c(1,5,10)]))),col=clrssat[c(1,5,10)],lty=1)
    plot(test$time_max,test$max_infect,col=alpha(color.gradient(test$inf,c("blue","red")),.2),main="inflection point",pch=20,xlab="Time Max",ylab="Max Infected")
    legend("topright",legend=paste("point=",round(unique(test$inf)[c(1,5,10)],digit=1)),col=clrssat[c(1,5,10)],lty=1)
    plot(test$time_max,test$max_infect,col=alpha(color.gradient(test$pind,c("blue","red")),.2),main="indiv learning",pch=20,xlab="Time Max",ylab="Max Infected")
    legend("topright",legend=paste("P indiv learning=",unique(test$pind)),col=clrssat[c(1,5,10)],lty=1)
    dev.off()

png("5e8779c56517890001536101/figures/alldim.png",width=1000,height=1200,pointsize=25)
    par(mfrow=c(3,2),cex=.5)
test=allresults
    plot(test$time_max,test$max_infect,col=alpha(color.gradient(log(test$sat),c("blue","yellow", "red")),.2),main="steepness",pch=20,xlab="Time Max",ylab="Max Infected")
	sset=sort(log10(test$sat))
	cols=alpha(color.gradient(sset,c("blue","yellow", "red")),1)
	lsset=seq(1,length(sset),length.out=4)
    legend("topright",legend=paste("steepness=10^",round(sset[lsset])),col=cols[lsset],pch=20)
    plot(test$time_max,test$max_infect,col=alpha(color.gradient(log(test$sat_r),c("blue","yellow", "red")),.2),main="revert steepness",pch=20,xlab="Time Max",ylab="Max Infected")
    legend("topright",legend=paste("steepness=10^",round(sset[lsset])),col=cols[lsset],pch=20)
    plot(test$time_max,test$max_infect,col=alpha(color.gradient(test$inf,c("blue","yellow", "red")),.2),main="inflection point",pch=20,xlab="Time Max",ylab="Max Infected")
	sset=sort(test$inf)
	cols=alpha(color.gradient(sset,c("blue","yellow", "red")),1)
	lsset=seq(1,length(sset),length.out=4)
    legend("topright",legend=paste("inp",round(sset[lsset],digit=2)),col=cols[lsset],pch=20)
    plot(test$time_max,test$max_infect,col=alpha(color.gradient(test$inf_r,c("blue","yellow", "red")),.2),main="revert inif. point",pch=20,xlab="Time Max",ylab="Max Infected")
    legend("topright",legend=paste("inp",round(sset[lsset],digit=2)),col=cols[lsset],pch=20)
    plot(test$time_max,test$max_infect,col=alpha(color.gradient(test$pind,c("blue","yellow", "red")),.2),main="individual learning",pch=20,xlab="Time Max",ylab="Max Infected")
	sset=sort(test$pind)
	cols=alpha(color.gradient(sset,c("blue","yellow", "red")),1)
	lsset=seq(1,length(sset),length.out=4)
    legend("topright",legend=paste("ind. learn.",round(sset[lsset],digit=2)),col=cols[lsset],pch=20)
    plot(test$time_max,test$max_infect,col=alpha(color.gradient(test$scores,c("blue","yellow", "red")),.2),main="score",pch=20,xlab="Time Max",ylab="Max Infected")
	sset=sort(test$scores)
	cols=alpha(color.gradient(sset,c("blue","yellow", "red")),1)
	lsset=seq(1,length(sset),length.out=4)
    legend("topright",legend=paste("score",round(sset[lsset],digit=2)),col=cols[lsset],pch=20)
dev.off()


    par(mfrow=c(1,3),cex=.5)
    plot(bscore(test$time_max),bscore(test$max_infect),col=alpha(color.gradient(log(test$sat),c("blue","red")),.2),main="steepness",pch=20)
    plot(bscore(test$time_max),bscore(test$max_infect),col=alpha(color.gradient(test$inf,c("blue","red")),.2),main="point",pch=20)
    plot(bscore(test$time_max),bscore(test$max_infect),col=alpha(color.gradient(test$pind,c("blue","red")),.2),main="indiv learning",pch=20)

    plot(test$sat,test$scores,log="x")



    test=allresults[allresults$score > 600 & allresults$max_infect < 250,]
    test=allresults[allresults$score <1,]
    par(mfrow=c(1,3),cex=.5)
    plot(test$time_max,test$max_infect,col=alpha(color.gradient(log(allresults$sat),c("blue","red")),.2),main="steepness",pch=20)
    plot(test$time_max,test$max_infect,col=alpha(color.gradient(allresults$inf,c("blue","red")),.2),main="point",pch=20)
    plot(test$time_max,test$max_infect,col=alpha(color.gradient(allresults$pind,c("blue","red")),.2),main="indiv learning",pch=20)
    plot(test$sat,test$inf,log="x",col=alpha(color.gradient(test$pind,c("blue","red")),.2),main="indiv learning",pch=20)


repetBest=sapply(1:100,function(i){print(i);abmSIR(poptest,1500,p=c(1,.2),di=2,i0=1,inf=.2,sat=1,xsize=xsize,ysize=ysize,visu=F,ap=F,ts=T,file=F,p_i=.5,log=F)$timeseries[,2]})
repetBestDos=parSapply(cl,1:100,function(i){print(i);abmSIR(poptest,1500,p=c(1,.2),di=2,i0=1,inf=.9,sat=1000,xsize=xsize,ysize=ysize,visu=F,ap=F,ts=T,file=F,p_i=.5,log=F)$timeseries[,2]})
neutralbad=sapply(1:10,function(i){print(i);abmSIR(poptest,1500,p=c(1,1),di=2,i0=1,inf=9,sat=1000,xsize=xsize,ysize=ysize,visu=F,ap=F,ts=T,file=F,p_i=1,log=F)$timeseries[,2]})
cl <- makeForkCluster(3,outfile="")
neutralGood=parSapply(cl,1:100,function(i){print(i);abmSIR(poptest,1500,p=c(.1,.1),di=2,i0=1,inf=9,sat=1000,xsize=xsize,ysize=ysize,visu=F,ap=F,ts=T,file=F,p_i=1,log=F)$timeseries[,2]})
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
legend("topright",legend=c("All bad","All good"),col=c("red","green"),lwd=2)
apply(neutralGood,2,function(i)lines(i,col=alpha("green",.2),lwd=2))
apply(neutralbad,2,function(i)lines(i,col=alpha("red",.1),lwd=2))
plot(1,1,type="n",ylim=c(0,500),xlim=c(1,1100),main="Learning")
apply(repetBest,2,function(i)lines(i,col=alpha("green",.2),lwd=2))
apply(repetBestDos,2,function(i)lines(i,col=alpha("red",.1),lwd=2))
legend("topright",legend=c("Revert","No revert"),col=c("red","green"),lwd=2)
dev.off()



png("5e8779c56517890001536101/figures/fullTrajHDR.png",width=800,height=400,pointsize=17)
par(mfrow=c(1,2))
hdr.boxplot(neutralGoodList[subset],pch=".",outline=F,,prob=c(50,95),space=0,ylab="#infected",border=NA,h=10,col="green",ylim=c(0,500),main="No learning")
legend("topright",legend=c("All bad","All good"),fill=c("red","green"))
par(new=T)
hdr.boxplot(neutralbadList[subset],pch=".",outline=F,,prob=c(50,95),space=0,ylab="#infected",border=NA,h=10,col="red",ylim=c(0,500))
axis(1,at=acs,subset[acs])


hdr.boxplot(repetBestList[subset],pch=".",outline=F,,prob=c(50,95),space=0,ylab="#infected",border=NA,h=10,col="green",ylim=c(0,500),main="Learning")
legend("topright",legend=c("Revert","No revert"),fill=c("red","green"))
par(new=T)
hdr.boxplot(repetBestDosList[subset],pch=".",outline=F,,prob=c(50,95),space=0,ylab="#infected",border=NA,h=10,col="red",ylim=c(0,500))
axis(1,at=acs,subset[acs])
dev.off()


par(mfrow=c(1,2))
hdr.boxplot(neutralGoodList[subset],pch=".",outline=F,,prob=c(50,95),space=0,ylab="#infected",border=NA,h=10,col="green",ylim=c(0,500))
par(new=T)
hdr.boxplot(repetBestList[subset],pch=".",outline=F,,prob=c(50,95),space=0,ylab="#infected",border=NA,h=10,col="green",ylim=c(0,500))
axis(1,at=acs,subset[acs])


hdr.boxplot(neutralbadList[subset],pch=".",outline=F,,prob=c(50,95),space=0,ylab="#infected",border=NA,h=10,col="red",ylim=c(0,500))
par(new=T)
hdr.boxplot(repetBestDosList[subset],pch=".",outline=F,,prob=c(50,95),space=0,ylab="#infected",border=NA,h=10,col="red",ylim=c(0,500))
axis(1,at=acs,subset[acs])

rsync -avz --info=progress2 --exclude "*simu_*.bin" volos:/home/share/simon/abmEpi/fullRandom .

aa=lapply(unlist(lapply(list.dirs("fullRandom"),list.files,pattern="allresults.bin",full.names=T)),function(f){load(f);return(allresults)})
allresults=do.call("rbind",aa)
allresults=allresults[-which(allresults$max_infect <4),]
allresults$scores=(1-bscore(allresults$time_max) + bscore(allresults$max_infect))/2


