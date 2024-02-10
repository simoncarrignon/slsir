
outfold='out3'
allres=list.files(outfold,pattern="subex*")
n.house=round(10^seq(1,2,length.out=15))
d.og=round(10^seq(-2,-.7,length.out=15),digits=3)
allparam=expand.grid(round(10^seq(1,2,length.out=15)),round(10^seq(-2,-.7,length.out=15),digits=3))
allrep=c()
for(i in 1:nrow(allparam))
   allrep=rbind(allrep,cbind(allparam[rep(i,500),],t(sapply(readRDS(file.path(outfold,list.files(outfold,pattern=paste0("subex_",i,"_.*")))),function(o)o[900,]) )))

image(log(tapply(allrep[,4],allrep[,c(1,2)],mean)))
image(log(tapply(allrep[,4],allrep[,c(1,2)],mean)),xaxt="n",yaxt="n",xlab="house per groups",ylab="q.og")
axis(1,at=seq(0,1,length.out=5),labels =seq(unique(allparam[,1])[1],unique(allparam[,1])[length(unique(allparam[,1]))],length.out=5))
axis(2,at=seq(0,1,length.out=5),labels =seq(unique(allparam[,2])[1],unique(allparam[,2])[length(unique(allparam[,2]))],length.out=5))
axis(1)
image()
image(log(tapply(allrep[,4],allrep[,c(1,2)],mean)),xaxt="n",yaxt="n",xlab="house per groups",ylab="q.og")
matres=tapply(allrep[,4],allrep[,c(1,2)],function(i)sum(i==0))
image(matres,x=n.house,y=d.og,log="xy")
plot(1,1,type="n",xlim=c(1,length(d.og)),ylim=c(0,1),xlab="p.og",ylab="% extinction",xaxt="n")
clrs=heat.colors(length(n.house))
for(i in 1:length(n.house)){
    lines(matres[i,]/500,col=clrs[i],lwd=4)
}
tick=seq(1,15,length.out=5)
axis(1,at=tick,label= paste0("10^",round(seq(-2,-.7,length.out=15),digit=1))[tick])
axis(1,at=tick,label= round(d.og,digit=2)[tick])
legend("topright",legend=n.house[seq(1,length(n.house),length.out=5)],lwd=2,col=clrs[seq(1,length(n.house),length.out=5)],title="house per groups")


n.group=unique(round(10^seq(.7,2.25,length.out=30)))
d.oq=10^seq(-3,-.6,length.out=30)
allparam=expand.grid(n.group=n.group,q=d.oq)

allresults=readRDS(file.path("Ngroup.RDS"))

allres=c()
for(i in 1:nrow(allparam)){
    allres=rbind(allres,cbind(allparam[i,],t(allresults[[i]])))
}

matres=tapply(allres[,4],allres[,c(1,2)],function(i)sum(i==0))
pdf("heatmapNgroupVpog.pdf")
image(matres,x=n.group,y=d.oq,log="xy",main="% of exctintion (dark red=100&, yellow=0%) with 10 house per groups",ylab="proba link group",xlab="# of groups")
dev.off()
pdf("ngroupVspog.pdf")
plot(1,1,type="n",xlim=c(1,length(d.oq)),ylim=c(0,1),xlab="p.og",ylab="% extinction",xaxt="n")
clrs=.colors(length(n.house))
for(i in rev(1:length(n.house))){
    lines(matres[i,]/1000,col=clrs[i],lwd=4)
}
d.oq=10^seq(-3,-.6,length.out=25)
n.house=unique(round(10^seq(.7,2.25,length.out=25)))
allparam=expand.grid(n.house=nhouse,q=d.oq)
tick=seq(1,length(d.oq),length.out=5)
axis(1,at=tick,label= round(10^seq(-3,-.6,length.out=25)[tick],digit=3))
legend("bottomleft",legend=n.house[seq(1,length(n.house),length.out=5)],lwd=2,col=clrs[seq(1,length(n.house),length.out=5)],title="house per groups")
dev.off()



#House per clus vs DOQ

allparam=expand.grid(unique(round(10^seq(.7,2.25,length.out=30))),10^seq(-3,-.6,length.out=30))

houseperclus=unique(round(10^seq(.7,2.25,length.out=30)))
d.oq=10^seq(-3,-.6,length.out=30)
allparam=expand.grid(houseperclus=houseperclus,q=d.oq)

allresults=readRDS(file.path("final.RDS"))

allres=c()
for(i in 1:nrow(allparam)){
    allres=rbind(allres,cbind(allparam[i,],t(allresults[[i]])))
}

matres=tapply(allres[,4],allres[,c(1,2)],function(i)sum(i==0))
pdf("heatmapHousePergroupVpog.pdf")
image(matres,x=houseperclus,y=d.oq,log="xy",main="% of exctintions (dark red=100%, yellow=0%) for 10 groups",xlab="house/groups",ylab="proba link groups")
dev.off()

pdf("pergroupVspog.pdf")
plot(1,1,type="n",xlim=c(1,length(d.oq)),ylim=c(0,1),xlab="p.og",ylab="% extinction",xaxt="n")
clrs=topo.colors(length(n.house))
for(i in rev(1:length(n.house))){
    lines(matres[i,]/500,col=clrs[i],lwd=4)
}
    lines(matres[6,]/500,col=1,lwd=4,lty=3)
d.oq.lab=10^seq(-3,-.6,length.out=25)
nhouselab=unique(round(10^seq(.7,2.25,length.out=25)))
tick=seq(1,length(d.oq),length.out=5)
#axis(1,at=tick,label= paste0("10^",seq(-3,-.6,length.out=25))[tick])
axis(1,at=tick,label=round(10^seq(-3,-.6,length.out=25),digit=3)[tick])
legend("topright",legend=nhouselab[seq(1,length(nhouselab),length.out=5)],lwd=2,col=clrs[seq(1,length(nhouselab),length.out=5)],title="house per groups",bg="white")
dev.off()

plot(1,1,type="n",xlim=c(1,length(d.oq)),ylim=c(0,1),xlab="p.og",ylab="% extinction",xaxt="n")
clrs=heat.colors(length(d.oq))
for(i in rev(1:length(d.oq))){
    lines(matres[,i]/500,col=clrs[i],lwd=4)
}
d.oq=10^seq(-3,-.6,length.out=25)
n.house=unique(round(10^seq(.7,2.25,length.out=25)))
allparam=expand.grid(n.house=nhouse,q=d.oq)
tick=seq(1,length(d.oq),length.out=5)
#axis(1,at=tick,label= paste0("10^",seq(-3,-.6,length.out=25))[tick])
axis(1,at=tick,label= round(10^seq(-3,-.6,length.out=25),digit=3)[tick])
legend("topright",legend=n.house[seq(1,length(n.house),length.out=5)],lwd=2,col=clrs[seq(1,length(n.house),length.out=5)],title="house per groups")
dev.off()
