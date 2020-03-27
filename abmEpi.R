

#' @param N the number of agents
#' @param tstep the duration of the simualtion
#' @param p the probability for one agent to transmit the disease to another one
#' @param di the distance between to agent under which the disease cna be transmitted
#' @param i0 the number of initial infection
#' @param speed the speed of the agents
abmSIR <- function(pop,tstep,p=1,i0=1,di=2,remi=10,speed=.8,xsize=100,ysize=100,visu=FALSE,inf=.5,sat=10){
    S=1
    I=2
    R=3
    G=2
    B=1

    
    if(is.null(dim(pop))){
        N=pop
        pop=cbind(x=runif(N,0,xsize),y=runif(N,0,ysize)) #generate a population 
        health=rep(S,N) #set all population as Susceptible to be infected
        infect=sample(N,i0) #choose the first random individuals to be infected
        health[infect]=I
        policies=sample(c(G,B),N,replace=T)
        ages=sample(c(G,B),N,replace=T)
        pop=cbind(pop,health=health,policies=policies)
    }
    N=nrow(pop)
    timeseries=c()
    infect=sample(N,i0) #choose the first random individuals to be infected
    pop[,"health"][infect]=I
    for(t in 1:tstep){
        pop[,"x"] = pop[,"x"]+rnorm(N,0,speed)
        pop[,"y"] = pop[,"y"]+rnorm(N,0,speed)
        pop[,"y"][pop[,"y"]>ysize]=ysize
        pop[,"x"][pop[,"x"]>xsize]=xsize
        pop[,"y"][pop[,"y"]<0]=0
        pop[,"x"][pop[,"x"]<0]=0
        
        infected=table(pop[,"health"],pop[,"ages"])
        #print(infected)

        for(i in which(pop[,"health"] == S)){#for each individual S we check if the are close enought to an I individual
            ind=pop[i,]

            ### Policies and Behavioral changes 

            if(ind["policies"] == B){
                group_infection=infected[2,ind["ages"]]/sum(infected[,ind["ages"]]) #compute the percentage of infect people from the same group 
                proba_switch=sig(group_infection,a=10,b=.5)
                if(runif(1)<proba_switch)
                    pop[i,"policies"]=G
            }

            p_ind = p[ind["policies"]]

            ### disease spread
            dist=sqrt(abs(pop[,"x"]-ind["x"])^2+abs(pop[,"y"]-ind["y"])^2) #check the distance of the all other agents
            ni=pop[,"health"][dist<di]==I #find infected neighbours

            if(any(ni)){  #if some agent are close enough 
                if(any(runif(sum(ni))<p_ind))pop[,"health"][i]=I #if one of the neighbours transmit the virus, the agent becomes Infected
            }
        }

        timeseries=rbind(timeseries,c(table(factor(pop[,"health"],levels=1:3)),table(factor(pop[,"policies"],levels=1:2))))#store the ratio S vs I


        ## a few lines of code ot visual the spread
        if(visu){
            par(mfrow=c(1,2))
            plot(pop[,"x"],pop[,"y"],pch=21,bg=pop[,"health"]-1,ylim=c(0,ysize),lwd=.2,xlim=c(0,xsize),xlab="",ylab="")
            plot(1:t,timeseries[,2],col="red",type="l",ylim=c(0,N),xlab="time",ylab="# infected")
            lines(1:t,timeseries[,5],col="blue",type="l",ylim=c(0,N),xlab="time",ylab="# GOOD")
        }
    }
    return(list(timeseries=timeseries,pop=pop))
}


generatePopulation <- function(N,agedistrib,policies){
    B=1
    G=2
    agedistrib=c(.24,.09,.12,.26,.13,.16)
    names(agedistrib)=letters[1:length(agedistrib)]
    ages=rep(names(agedistrib),agedistrib*N)
    if(sum(prop.table(table(ages))-agedistrib))stop()
    pop=cbind(x=runif(N,0,xsize),y=runif(N,0,ysize)) #generate a population 
    health=rep(S,N) #set all population as Susceptible to be infected
    policies=rep(B,N)
    pop=cbind(pop,health=health,policies=policies,ages=as.factor(ages))
    return(pop)
}


sig<-function(x,a=10,b=.5)1/(1+exp(-a*(x-b)))


updateBehavior <- function(
#try one run with visualisation:
 res=abmSIR(200,1000,p=1,di=3,visu=T)

