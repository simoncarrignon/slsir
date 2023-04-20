
aaa=readRDS("result_good_n200_p0.25.rds")
pdf("response_threelayout.pdf",width=12,height=7)
par(mfrow=c(1,2))
boxplot(apply(aaa,1,function(o)sapply(o,function(i)max(i[[2]]))),ylab="max infected",xlab="layout",outline=F)
boxplot(apply(aaa,1,function(o)sapply(o,function(i)which.max(i[[2]]))),ylab="time to max infected",xlab="layout",outline=F)
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

plot(1,1,xlim=c(0,100),ylim=c(0,80),type="n")
for(i in 3:1)sapply(lapply(aaa[i,],"[[",2),lines,col=i+1)
