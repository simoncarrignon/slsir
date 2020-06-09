require(RColorBrewer)

visualize <- function(allpop,timeseries,xsize,ysize,file=F){
    if(file)png(sprintf("frame_%04d.png",nrow(timeseries)),width=1200,height=900)
    pop=c()
    if(is.null(dim(allpop)))pop=allpop[[length(allpop)]]
    else pop=allpop

    if(is.null(dim(allpop)))layout(matrix(c(1,1,1,2,3,4),nrow=3,ncol=2),width=c(2,1))  
    else layout(matrix(c(1,1,2,3),nrow=2,ncol=2),width=c(2,1))  
    N=nrow(pop)
    plot(pop[,"x"],pop[,"y"],pch=20+pop[,"behavior"],col=pop[,"behavior"]+3,bg=pop[,"health"]-1,ylim=c(0,ysize),lwd=.8,xlim=c(0,xsize),xlab="",ylab="",cex=2,xaxt="n",yaxt="n")
    plot(1:nrow(timeseries),timeseries[,2],col="black",type="l",lwd=2,ylim=c(0,N),xlab="time",ylab="# infected")
    lines(1:nrow(timeseries),timeseries[,3],col="red",type="l",lwd=2)
    lines(1:nrow(timeseries),timeseries[,1],col="green",type="l",lwd=2)


    if(is.null(dim(allpop))){
        allheal=sapply(allpop,function(i){table(factor(i[,"health"],levels=1:3),i[,"ages"])[2,]})
        allbeh=sapply(allpop,function(i){table(factor(i[,"behavior"],levels=1:2),i[,"ages"])[2,]})
        agecol=brewer.pal(6,"Pastel1")
        barplot(allheal,border=NA,space=0,col=agecol)
        barplot(allbeh,border=NA,space=0,col=agecol)
        legend("toplef",legend=c("0-18","18-25","26-34","35-54","55-64","65+"),fill=agecol,title="age category",cex=1)
    }
else{
    plot(1:nrow(timeseries),timeseries[,4],col=4,type="l",lwd=2,ylim=c(0,N),xlab="time",ylab="# behavior")
    lines(1:nrow(timeseries),timeseries[,5],col=5,type="l",lwd=2)
}
    if(file)dev.off()
}

