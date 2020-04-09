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
#' @param ts count and return the number of users during the run
#' @param ap keep and return the full population for each time step
abmSIR <- function(pop,tstep,p=1,i0=1,di=2,recovery=10,speed=.8,xsize=100,ysize=100,visu=FALSE,inf=.5,sat=10,sat_r=10000,inf_r=1.1,log=F,checkcountact=F,ts=T,ap=F,p_i=1,file=F){

    if(is.null(dim(pop))) #if pop is a unique number (ie not preinitialized) 
        pop=generatePopulation(N=pop,xsize=xsize,ysize=ysize,recovery=recovery,speed=speed)

    N=nrow(pop)

    infect=sample(N,i0) #choose the first random individuals to be infected
    pop[,"health"][infect]=I

    if(checkcountact){
        meancontact=c(0)
        contacts=rep(0,N)
    }
    if(ts)timeseries=c(table(factor(pop[,"health"],levels=1:3)),table(factor(pop[,"behavior"],levels=1:2)))
    if(ap)allpop=list()
    output=list()
    for(t in 1:tstep){

        if(log)print(paste0("tstep:",t))
        ##move the agents 
        #pop[,"x"] = pop[,"x"]+runif(N,0,2*speed)-speed
        #pop[,"y"] = pop[,"y"]+runif(N,0,2*speed)-speed
        dir=runif(N)*2*pi
        pop[,"x"] = pop[,"x"]+pop[,"speed"] * cos(dir)
        pop[,"y"] = pop[,"y"]+pop[,"speed"] * sin(dir)
        pop[,"y"][pop[,"y"]>ysize]=ysize
        pop[,"x"][pop[,"x"]>xsize]=xsize
        pop[,"y"][pop[,"y"]<0]=0
        pop[,"x"][pop[,"x"]<0]=0

        pop[pop[,"health"] == I ,"recovery"]=pop[pop[,"health"] == I ,"recovery"]-1
        pop[pop[,"recovery"] < 1,"health"]=R

        #count effected by agents
        infected=table(factor(pop[,"health"],levels=1:3),pop[,"ages"])
        #behavior=table(factor(pop[,"behavior"],levels=1:2),pop[,"ages"])

        if(checkcountact){
            contacts=contacts+sapply(1:nrow(pop),function(i)sum(sqrt(abs(pop[-i,"x"]-pop[i,"x"])^2+abs(pop[-i,"y"]-pop[i,"y"])^2)<di))
            meancontact=c(meancontact,mean(contacts))
        }

        for(i in sample(N)){#for each individual S we check if the are close enought to an I individual
            ind=pop[i,]

            ### Policies and Behavioral changes 


            if(runif(1)<p_i){ #probability for individual learning (p_i)
                group_infection=infected[2,ind["ages"]]/sum(infected[,ind["ages"]]) #compute the percentage of infected people from the same group 

                if(ind["behavior"] == B){
                    proba_switch=sig(group_infection,a=sat,b=inf)
                    if(runif(1)<proba_switch)
                        pop[i,"behavior"]=G
                }
                else{
                    proba_switch=sig(1-group_infection,a=sat_r,b=inf_r)
                    if(runif(1)<proba_switch)
                        pop[i,"behavior"]=B
                }

            }
            else{ #probability of social learning (1-p_i)
                pop[i,"behavior"]=sample(pop[pop[,"ages"]==ind["ages"],"behavior"],1)
            }

            p_ind = p[ind["behavior"]]

            ### disease spread
            dist=sqrt(abs(pop[,"x"]-ind["x"])^2+abs(pop[,"y"]-ind["y"])^2) #check the distance of the all other agents
            ni=pop[,"health"][dist<di]==I #find infected neighbours

            if(any(ni)){  #if some agent are close enough 
                if(any(runif(sum(ni))<p_ind))pop[,"health"][i]=I #if one of the neighbours transmit the virus, the agent becomes Infected
            }
        }
        if(ts)timeseries=rbind(timeseries,c(table(factor(pop[,"health"],levels=1:3)),table(factor(pop[,"behavior"],levels=1:2))))#store the ratio S vs I
        if(ap)allpop[[t]]=pop
		if(visu){
            if(ap)
                visualize(allpop,timeseries,xsize=xsize,ysize=ysize,file=file)
            else 
                visualize(pop,timeseries,xsize=xsize,ysize=ysize,file=file)
        }
	}
    if(ts)output$timeseries=timeseries
    if(ap)output$allpop=allpop
    if(checkcountact){
        output$meancontact=meancontact
        output$contacts=contacts
    }
    return(output)
}


#' @param N number of agents
#' @param agedistrib a distribution defining the percentage of population represented in different age class
#' @param behavior a distribution of behaviors
#' @param xsize spatial limits (to keep in the model?)
#' @param ysize spatial limits 

generatePopulation <- function(N,agedistrib=NULL,behavior=NULL,xsize=100,ysize=100,recovery=NULL,speed=NULL){
    if(is.null(agedistrib)){
       agedistrib=c(.24,.09,.12,.26,.13,.16) #source: https://www.kff.org/other/state-indicator/distribution-by-age/
       names(agedistrib)=letters[1:length(agedistrib)]
    }
	if(is.null(rep))stop()
    ages=rep(names(agedistrib),agedistrib*N)
    if(sum(prop.table(table(ages))-agedistrib))stop()

    if(is.null(behavior))behavior=rep(B,N) #by default start with everyone as Bad

    pop=cbind(x=runif(N,0,xsize),y=runif(N,0,ysize)) #generate a population 
    health=rep(S,N) #set all population as Susceptible to be infected
    if(length(recovery)==1)recovery=rep(recovery,N) #set different recovery time for different individuals
    if(length(recovery)==2)recovery=runif(N,recovery[1],recovery[2]) #
    if(length(speed)==1)speed=rep(speed,N) #set different speed time for different individuals
    if(length(speed)==2)speed=abs(rnorm(N,speed[1],speed[2])) 
    pop=cbind(pop,health=health,behavior=behavior,ages=as.factor(ages),recovery=recovery,speed=speed)
    return(pop)
}


#'@param
#'@param a parameter to define the steepness of the sigmoid slope  
#'@param b parameter to define when the slope starts
sig<-function(x,a=10,b=.5)1/(1+exp(-a*(x-b)))

library(RColorBrewer)

visualize <- function(allpop,timeseries,xsize,ysize,file=F){
    if(file)png(sprintf("frame_%04d.png",nrow(timeseries)),width=1200,height=900)
    pop=c()
    if(is.null(dim(allpop)))pop=allpop[[length(allpop)]]
    else pop=allpop

    if(is.null(dim(allpop)))layout(matrix(c(1,1,1,2,3,4),nrow=3,ncol=2),width=c(2,1))  
    else par(mfrow=c(1,2))
    N=nrow(pop)
    plot(pop[,"x"],pop[,"y"],pch=20+pop[,"behavior"],bg=pop[,"health"]-1,ylim=c(0,ysize),lwd=.2,xlim=c(0,xsize),xlab="",ylab="",cex=2)
    plot(1:nrow(timeseries),timeseries[,2],col="red",type="l",ylim=c(0,N),xlab="time",ylab="# infected")
    lines(1:nrow(timeseries),timeseries[,1],col="green",type="l",ylim=c(0,N),xlab="time",ylab="# GOOD")
    lines(1:nrow(timeseries),timeseries[,3],col="blue",type="l",ylim=c(0,N),xlab="time",ylab="# GOOD")


    if(is.null(dim(allpop))){
        allheal=sapply(allpop,function(i){table(factor(i[,"health"],levels=1:3),i[,"ages"])[2,]})
        allbeh=sapply(allpop,function(i){table(factor(i[,"behavior"],levels=1:2),i[,"ages"])[2,]})
        agecol=brewer.pal(6,"Pastel1")
        barplot(allheal,border=NA,space=0,col=agecol)
        barplot(allbeh,border=NA,space=0,col=agecol)
        legend("toplef",legend=c("0-18","18-25","26-34","35-54","55-64","65+"),fill=agecol,title="age category",cex=1)
    }
    if(file)dev.off()
}
