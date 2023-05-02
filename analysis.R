
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
for(i in 3:1)sapply(lapply(groupedresults[i,],"[[",2),lines,col=i+1)
legend("topright",legend=paste0("layout ",1:3),lty=1,col=2:4)

groupedresults=readRDS("result_goodbis_n1000_p0.5.rds")
groupedresults=readRDS("result_limited_n50_p0.2.rds")
groupedresultsb=readRDS("result_limitedbis_n50_p0.2.rds")
hpeak=sapply(groupedresults,function(o)sapply(o,function(i)max(i[[2]])))
hpeakb=sapply(groupedresultsb,function(o)sapply(o,function(i)max(i[[2]])))
hpeak=rbind(hpeak,hpeakb)
tpeak=sapply(groupedresults,function(o)sapply(o,function(i)which.max(i[[2]])))
tpeakb=sapply(groupedresultsb,function(o)sapply(o,function(i)which.max(i[[2]])))
tpeak=rbind(tpeak,tpeakb)
severity= hpeak/max(hpeak)*(1-tpeak/max(tpeak))
par(mfrow=c(1,3));boxplot(hpeak,main="peak's height");boxplot(tpeak,main="time to peak");boxplot(severity,main="severity")

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

plot(1,1,xlim=c(0,50),ylim=c(0,100),type="n")
for(i in 3:1)sapply(lapply(groupedresults[[i]],"[[",2),lines,col=adjustcolor(i+1,.1),lwd=2)
legend("topright",legend=paste0("layout ",1:3),lty=1,col=2:4)
boxplot(sapply(1:3,function(n)sapply(groupedresults[[n]],function(i)max(i[[2]])/which.max(i[[2]]))))
