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

#allparameter=expand.grid(inf=inf,inf_r=inf_r,sat=sat,sat_r=sat_r,pind=pind)
allparameter=cbind.data.frame(inf=inf,sat=sat,inf_r=inf_r,sat_r=sat_r,pind=pind)

cl <- makeForkCluster(ns,outfile="")
allsummary=parSapply(cl,1:nsim,function(j)
                     {
                         print(paste("sim #",j,"/",nsim));
                         simu=abmSIR(500,1500,p=c(1,.1),di=2,i0=1,recovery=c(8,14)*25,speed=c(1,.2),xsize=xsize,ysize=ysize,
                                     inf=allparameter$inf[j],
                                     sat=allparameter$sat[j],
                                     inf_r=allparameter$inf_r[j],
                                     sat_r=allparameter$sat_r[j],
                                     p_i=allparameter$pind[j],
                                     ts=T,ap=F,visu=F
                                     )
                         #save(file=file.path(fold,paste0("simu_",j,".bin")),simu)
                         c(time_max=which.max(simu$timeseries[,2]),final_size=sum(simu$timeseries[1500,2:3]),max_infect=max(simu$timeseries[,2]))
                     }
)


stopCluster(cl)

allresults=cbind(allparameter,t(allsummary))
save(file=file.path(fold,"allresults.bin"),allresults)
new <- Sys.time() - old # calculate difference
print(new) # print in nice format


