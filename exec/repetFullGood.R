fold="testingSim"
pg=.5
behave=rep(c(G,B),500*c(pg,1-pg))
xsize=ysize=100

nsim=10
cl <- makeForkCluster(5,outfile="")
allsummary=parLapply(cl,1:10,function(j)
                     {
                         print(paste("sim #",j,"/",nsim));
                         pop=generatePopulation(N=500,xsize=xsize,ysize=ysize,recovery=c(8,14)*25,speed=c(1,.2),behavior=behave)
                         simu=abmSIR(1500,p=c(1,.1),di=2,i0=1,recovery=c(8,14)*25,speed=c(1,.2),xsize=xsize,ysize=ysize,
                                     pop=pop,
                                     inf=.5,
                                     sat=10,
                                     inf_r=.5,
                                     sat_r=100,
                                     p_i=.1,
                                     sl_rad=10,
                                     strategy="random",
                                     ts=T,ap=F,visu=T
                                     )
                         save(file=file.path(fold,paste0("allG_simu_",j,".bin")),simu)
                     }
)
stopCluster(cl)



pg=0
behave=rep(c(G,B),500*c(pg,1-pg))


cl <- makeForkCluster(5,outfile="")
allsummary=parLapply(cl,1:10,function(j)
                     {
                         print(paste("sim #",j,"/",nsim));
                         pop=generatePopulation(N=500,xsize=xsize,ysize=ysize,recovery=c(8,14)*25,speed=c(1,.2),behavior=behave)
                         simu=abmSIR(500,p=c(1,.1),di=2,i0=1,recovery=c(8,14)*25,speed=c(1,.2),xsize=xsize,ysize=ysize,
                                     pop=pop,
                                     inf=.5,
                                     sat=10,
                                     inf_r=.5,
                                     sat_r=10,
                                     p_i=.5,
                                     sl_rad=10,
                                     strategy="random",
                                     ts=T,ap=F,visu=F
                                     )
                         save(file=file.path(fold,paste0("allB_simu_",j,".bin")),simu)
                     }
)
stopCluster(cl)


