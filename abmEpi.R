S=1
I=2
R=3
G=2
B=1

sir=c(1:3)
names(sir)=c("S","I","R")

#' @param N the number of agents
#' @param tstep the duration of the simualtion
#' @param p the probability for one agent to transmit the disease to another one
#' @param di the distance between to agent under which the disease cna be transmitted
#' @param i0 the number of initial infection
#' @param speed the speed of the agents
#' @param reco time for agent to Recover
#' @param visu TRUE or FALSE, if output should be plotted
#' @param sat speed of saturation of the sigmoid
#' @param inf inflexion point of the sigmoid
abmSIR <- function(pop,tstep,p=1,i0=1,di=2,reco=10,speed=.8,xsize=100,ysize=100,visu=FALSE,inf=.5,sat=10){
    
    if(is.null(dim(pop))) #if pop is a unique number (ie not preinitialized) 
        pop=generatePopulation(N=pop)

    N=nrow(pop)

    infect=sample(N,i0) #choose the first random individuals to be infected
    pop[,"health"][infect]=I

    timeseries=c() #table to store output
    for(t in 1:tstep){

        ##move the agents 
        pop[,"x"] = pop[,"x"]+rnorm(N,0,speed)
        pop[,"y"] = pop[,"y"]+rnorm(N,0,speed)
        pop[,"y"][pop[,"y"]>ysize]=ysize
        pop[,"x"][pop[,"x"]>xsize]=xsize
        pop[,"y"][pop[,"y"]<0]=0
        pop[,"x"][pop[,"x"]<0]=0
        
        #count effected by agents
        infected=table(pop[,"health"],pop[,"ages"])

        for(i in which(pop[,"health"] == S)){#for each individual S we check if the are close enought to an I individual
            ind=pop[i,]

            ### Policies and Behavioral changes 

            if(ind["behavior"] == B){
                group_infection=infected[2,ind["ages"]]/sum(infected[,ind["ages"]]) #compute the percentage of infect people from the same group 
                proba_switch=sig(group_infection,a=sat,b=inf)
                if(runif(1)<proba_switch)
                    pop[i,"behavior"]=G
            }

            p_ind = p[ind["behavior"]]

            ### disease spread
            dist=sqrt(abs(pop[,"x"]-ind["x"])^2+abs(pop[,"y"]-ind["y"])^2) #check the distance of the all other agents
            ni=pop[,"health"][dist<di]==I #find infected neighbours

            if(any(ni)){  #if some agent are close enough 
                if(any(runif(sum(ni))<p_ind))pop[,"health"][i]=I #if one of the neighbours transmit the virus, the agent becomes Infected
            }
        }

        timeseries=rbind(timeseries,c(table(factor(pop[,"health"],levels=1:3)),table(factor(pop[,"behavior"],levels=1:2))))#store the ratio S vs I

        if(visu)visualize(pop,timeseries)
    }
    return(list(timeseries=timeseries,pop=pop))
}


#' @param N number of agents
#' @param agedistrib a distribution defining the percentage of population represented in different age class
#' @param behavior a distribution of behaviors
#' @param xsize spatial limits (to keep in the model?)
#' @param ysize spatial limits 

generatePopulation <- function(N,agedistrib=NULL,behavior=NULL,xsize=100,ysize=100){
    if(is.null(agedistrib)){
       agedistrib=c(.24,.09,.12,.26,.13,.16) #source: https://www.kff.org/other/state-indicator/distribution-by-age/
       names(agedistrib)=letters[1:length(agedistrib)]
    }
    ages=rep(names(agedistrib),agedistrib*N)
    if(sum(prop.table(table(ages))-agedistrib))stop()

    if(is.null(behavior))behavior=rep(B,N) #by default start with everyone as Bad

    pop=cbind(x=runif(N,0,xsize),y=runif(N,0,ysize)) #generate a population 
    health=rep(S,N) #set all population as Susceptible to be infected
    pop=cbind(pop,health=health,behavior=behavior,ages=as.factor(ages))
    return(pop)
}


#'@param
#'@param a parameter to define the steepness of the sigmoid slope  
#'@param b parameter to define when the slope starts
sig<-function(x,a=10,b=.5)1/(1+exp(-a*(x-b)))


visualize <- function(pop,timeseries){
            par(mfrow=c(1,2))
            plot(pop[,"x"],pop[,"y"],pch=21,bg=pop[,"health"]-1,ylim=c(0,ysize),lwd=.2,xlim=c(0,xsize),xlab="",ylab="")
            plot(1:nrow(timeseries),timeseries[,2],col="red",type="l",ylim=c(0,N),xlab="time",ylab="# infected")
            lines(1:nrow(timeseries),timeseries[,5],col="blue",type="l",ylim=c(0,N),xlab="time",ylab="# GOOD")
}
