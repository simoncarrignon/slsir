library(sf)
library(igraph)
tripillia=st_read("Ukraine_Shapefile/Megasite_houses.shp")

plot(tripillia[,"House_type"])
trip.graph=toGraph(tripillia)
plot(trip.graph)
plot(density(log(E(trip.graph)$weight)))
degree(trip.graph)
fac= delete.edges(trip.graph,E(trip.graph)[E(trip.graph)$weight>50])
V(fac)$color=mole$memberships[3,]
plot(fac)

plot(st_geometry(tripillia[1:200,]),col="white")
fewhouses=c(1,20,50,35)
plot(st_geometry(tripillia[fewhouses,]),col="red",add=T,lwd=1)
xyc=st_coordinates(st_centroid(tripillia[fewhouses,]))
text(xyc[,1],xyc[,2] ,labels=fewhouses,las=1,cex=3)
didi=st_distance(tripillia[fewhouses,])
rownames(didi)=fewhouses
colnames(didi)=fewhouses

E(trip.graph,P=fewhouses)$weight
E(trip.graph,P=fewhouses[c(1,2,1,3,1,4)])$weight
didi[as.character(fewhouses)[1],]

fac= delete.edges(trip.graph,E(trip.graph)[E(trip.graph)$weight>800])
V(fac)$color=mole$memberships[3,]

V(fac)$state
tripillia$color=V(fac)$state
plot(tripillia[,"color"])
plot(fac)
library(igraph)
library(sf)
#lapply(1:10,function(i){
print(i);
gt=fac
t=0;
crve=c();
n=nrow(tripillia)
while( sum(V(gt)$state<0)<.99*n && t < 1500 && sum(V(gt)$state>0)>0 ){
    gt=diffuse_culture(gt,1);
    ni=sum(V(gt)$state>0); 
    crve=c(crve,ni)
    print(paste(t,ni))
    t=t+1;
    saveRDS(file=paste0("fullgraph_t",t,".RDS"),gt)
    tripillia$color=V(gt)$state
    png(sprintf("tpl_t%04d.png",t),width=800,height=800)
    plot(st_centroid(tripillia[,"color"]),cex=5,pch=20,main=t)
    dev.off()
};

sample(V(trip.graph))
allreas=as.numeric(st_area(tripillia))

V(trip.graph)$type="H"
V(trip.graph)$type[allreas>180]="M"
fac= delete.edges(trip.graph,E(trip.graph)[E(trip.graph)$weight>800])

plot(st_geometry(tripillia),lwd=.2,col="red")
plot(st_geometry(tripillia)[as.numeric(allreas)>180],lwd=3,add=T,col="green")

V(gt)[locals]$color=NULL
u=1
for( meethouse in V(gt)[V(gt)$type=="M"]){
    locals=neighbors(gt, meethouse)
    V(gt)[locals]$comus=paste0(locals$comus,u)
    print(sum(V(gt)$state))
    u=u+1
    print(u)
}


gt=fac
prob_diffuse = .2
t=0

precomprob=list()
for( house in V(gt)[V(gt)$type=="H"]){
    pm=neighbors(gt, house)
    pm=pm[pm$type=="M"]
    probas=sapply(seq_along(pm),function(i)E(gt,P=c(house,pm[i]))$weight)
    precomprob[[as.character(house)]]=probas
}

while( sum(V(gt)$state<0)<.90*n && t < 1500 && sum(V(gt)$state>0)>0 ){
    print(paste("step",t))
    V(gt)$comus=""
    allmeethouse=list()
    for( house in V(gt)[V(gt)$type=="H"]){
        pm=neighbors(gt, house)
        pm=pm[pm$type=="M"]
        probas=precomprob[[as.character(house)]]
        meethouse=sample(pm,1,prob=1/(probas^3))
        allmeethouse[[as.character(meethouse)]]=c(allmeethouse[[as.character(meethouse)]],house)
        V(gt)$comus[as.numeric(house)]=as.character(meethouse)
    }


    for(mt in allmeethouse){
        if(any(V(gt)[mt]$state>0)){
            infected=runif(length(mt))<prob_diffuse
            ni=V(gt)[mt]$state==0
            infect=ni & infected
            V(gt)$state[mt][infect]=10
        }

    }

    V(gt)$state[V(gt)$state==1]=-1 
    V(gt)$state[V(gt)$state>0]=V(gt)$state[V(gt)$state>0]-1 

    tripillia$epicol="red"
    tripillia$epicol=ifelse(V(gt)$state>0,"red","white")
    tripillia$comucol=V(gt)$comus
    png(sprintf("mh_tpl_t%04d.png",t),width=800,height=800)
    plot(st_centroid(tripillia[,c("epicol","comucol")]),cex=2,pch=20,main=t,bgc="black")
    dev.off()
    png(sprintf("mh_tpl_bis_t%04d.png",t),width=800,height=800)
    par(mfrow=c(1,2))
    plot(st_coordinates(st_centroid(tripillia[,c("epicol","comucol")])),axes=F,ann=F,bg=tripillia$epicol,pch=21,cex=2)
    plot(st_coordinates(st_centroid(tripillia[,c("epicol","comucol")])),axes=F,ann=F,bg=as.factor(tripillia$comucol),pch=21,cex=2)
    dev.off()
    print("done")
    t=t+1
}

TRIp.graph(V(trip.graph)[1]
