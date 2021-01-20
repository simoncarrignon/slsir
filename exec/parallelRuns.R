old <- Sys.time() 

args=commandArgs(trailingOnly = TRUE) #pass number of slave by comand line, should be #node - 1 as one node is used by the master 

ns=args[1]#first argument is the number of slave
nsm=args[2]#second argument the number of simulation 
mainfold=args[3] #third argument = name of the folder wher to store the results

if(is.na(mainfold) | mainfold=="") mainfold=Sys.info()['nodename']

fi=0
fold=paste0(mainfold,fi)
while(file.exists(fold)){
    fi=fi+1
    fold=paste0(mainfold,fi)
}

print(paste0("Abc will be stored in folder: ",fold))
dir.create(fold)


source("abmEpi.R")
library(parallel)
xsize=ysize=100

nsim=nsm

inf=runif(nsim,0,1)
inf_r=runif(nsim,0,1)
sat=10^runif(nsim,-1,3)
sat_r=10^runif(nsim,-1,3)
pind=runif(nsim)
sl_rad=sample(100,nsim,replace=T)

#allparameter=expand.grid(inf=inf,inf_r=inf_r,sat=sat,sat_r=sat_r,pind=pind)
allparameter=cbind.data.frame(inf=inf,sat=sat,inf_r=inf_r,sat_r=sat_r,pind=pind,sl_rad=sl_rad)


pg=0
behave=rep(c(G,B),500*c(pg,1-pg))

cl <- makeForkCluster(ns,outfile="")
allsummary=parSapply(cl,1:nsim,function(j)
                     {
                         print(paste("sim #",j,"/",nsim));
                         pop=generatePopulation(N=500,xsize=xsize,ysize=ysize,recovery=c(8,14)*25,speed=c(1,.2),behavior=behave)
                         simu=abmSIR(1500,p=c(1,.1),di=2,i0=1,recovery=c(8,14)*25,speed=c(1,.2),xsize=xsize,ysize=ysize,
                                     pop=pop,
                                     inf=allparameter$inf[j],
                                     sat=allparameter$sat[j],
                                     inf_r=allparameter$inf_r[j],
                                     sat_r=allparameter$sat_r[j],
                                     p_i=allparameter$pind[j],
                                     sl_rad=allparameter$sl_rad[j],
                                     strategy="best",
                                     ts=T,ap=F,visu=F,bt=0
                                     )
                         #save(file=file.path(fold,paste0("simu_",j,".bin")),simu)
                         max_infect=max(simu$timeseries[,2])
                         max_infect150=max(simu$timeseries[1:150,2])
                         max_infect250=max(simu$timeseries[1:250,2])
                         c(
                           max_infect=max_infect,
                           max_infect150=max_infect150,
                           max_infect250=max_infect250,
                           time_max=which.max(simu$timeseries[,2]),
                           time_max2=which.max(simu$timeseries[,2]>=(max_infect/2)),
                           time_max4=which.max(simu$timeseries[,2]>=(max_infect/4)),
                           time_max150=which.max(simu$timeseries[,2]>=(max_infect150)),
                           time_max250=which.max(simu$timeseries[,2]>=(max_infect250)),
                           final_size=sum(simu$timeseries[1500,2:3]),
                           size150=sum(simu$timeseries[150,2:3]),
                           size250=sum(simu$timeseries[250,2:3])
                           )
                     }
)


stopCluster(cl)

allresults=cbind(allparameter,t(allsummary))
save(file=file.path(fold,"allresults.bin"),allresults)
new <- Sys.time() - old # calculate difference
print(new) # print in nice format


