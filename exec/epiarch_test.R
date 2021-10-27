library(devtools)
load_all(".")


#Simple exploration of the main function of the package:


poptest=generatePopulation(100,recovery=2500,speed=.8)  #generate a simple population

### position agents around three different subpopulations
v=sample(1:3,nrow(poptest),replace=T) 
poptest[v==1,"x"]=10
poptest[v==1,"y"]=10
poptest[v==2,"x"]=90
poptest[v==2,"y"]=50
poptest[v==3,"x"]=10
poptest[v==3,"y"]=90
poptest[,"x"]=apply(poptest,1,function(l)l["x"]=l["x"] +rnorm(1,0,2))
poptest[,"y"]=apply(poptest,1,function(l)l["y"]=l["y"] +rnorm(1,0,2))
poptest[,"speed"]=rnorm(100,.5,.1)
#######

##create 3 "travelers", that will travel in between populations 
travelers=sample.int(100,3)
poptest[travelers,"speed"]=5


# run simulation
slsirSimu(poptest,2500,p=c(1,1),di=2,i0=1,visu=T)


