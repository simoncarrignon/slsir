
groupedresults=readRDS("result_goodbis_n1000_p0.5.rds")
groupedresults=readRDS("result_gis_n1000_p0.2.rds")
pdf("response_threelayout.pdf",width=12,height=7)
par(mfrow=c(1,2))
boxplot(apply(groupedresults,1,function(o)sapply(o,function(i)max(i[[2]]))),ylab="max infected",xlab="layout",outline=F)
boxplot(apply(groupedresults,1,function(o)sapply(o,function(i)which.max(i[[2]]))),ylab="time to max infected",xlab="layout",outline=F)
dev.off()
#Final result after 25 time step:
pdf("network.pdf",width=12,height=4)

par(mfrow=c(1,3))
par(mar=c(0,0,0,0))
for(i in 1:3){
    plot(0,0,xlim=c(-1.5,1.5),ylim=c(-1.5,1.5),type="n",ann=F,axes=F)
    plot(graphs[[i]],layout=st_coordinates(scenarios[[i]]),add=T,rescale=F)
}
dev.off()

plot(1,1,xlim=c(0,100),ylim=c(0,800),type="n")
for(i in 3:1)sapply(lapply(groupedresults[i,],"[[",2),lines,col=categorical_pal(3)[i])
legend("topright",legend=paste0("layout ",1:3),lty=1,col=categorical_pal(3))

groupedresults=readRDS("result_goodbis_n1000_p0.5.rds")
groupedresults=readRDS("result_limited_n500_p0.6.rds")

par(mfrow=c(2,3));
for(n in c(300,3000)){
threepro=lapply(c(2,4,6),function(p){
groupedresults=readRDS(paste0("result_standard_n",n,"_p0.",p,".rds"))
hpeak=sapply(groupedresults,function(o)sapply(o,function(i)max(i[[2]])))
started=apply(hpeak,2,function(i)which(i>10))
tpeak=sapply(groupedresults,function(o)sapply(o,function(i)which.max(i[[2]])))
severity= hpeak/max(unlist(hpeak))*(1-tpeak/max(tpeak))
return(list(tpeak=tpeak,hpeak=hpeak,severity=severity))
})

binded=lapply(names(threepro[[1]]),function(n)do.call("cbind",lapply(threepro,"[[",n)))
names(binded)=names(threepro[[1]])
binded=lapply(binded,function(i)i[,c(1,4,7,2,5,8,3,6,9)])
cls=rep(adjustcolor(categorical_pal(3),.6),c(3,3,3))
boxplot(binded[["hpeak"]],main="peak's height",col=cls,axes=F);boxplot(binded[["tpeak"]],main="time to peak",col=cls,axes=F);boxplot(binded[["severity"]],main="severity",col=cls,,axes=F)
}

pdf("response_threelayout.pdf",width=12,height=7)
par(mfrow=c(1,2))
boxplot(sapply(groupedresults,function(o)sapply(o,function(i)max(i[[2]]))),ylab="max infected",xlab="layout",outline=F)
boxplot(sapply(groupedresults,function(o)sapply(o,function(i)which.max(i[[2]]))),ylab="time to max infected",xlab="layout",outline=F)
dev.off()
#Final result after 25 time step:
pdf("network.pdf",width=12,height=4)

par(mfrow=c(1,3))
par(mar=c(0,0,0,0))
for(i in 1:3){
    plot(0,0,xlim=c(-1.5,1.5),ylim=c(-1.5,1.5),type="n",ann=F,axes=F)
    plot(graphs[[i]],layout=st_coordinates(scenarios[[i]]),add=T,rescale=F)
}
dev.off()

plot(1,1,xlim=c(0,200),ylim=c(0,700),type="n",xlab="time",ylab="#contagious individuals")
for(i in 3:1)sapply(lapply(groupedresults[[i]],"[[",2),lines,col=adjustcolor(categorical_pal(3)[i],.1),lwd=2)
legend("topright",legend=paste0("layout ",1:3),lty=1,col=adjustcolor(categorical_pal(3),.8),lwd=2)
