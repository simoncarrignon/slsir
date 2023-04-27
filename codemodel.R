library(sf)
library(igraph)

#' Create n points within a square positioned using h, v, and n.
#' 
#' @param h A numeric value indicating the horizontal distance of the square from the origin.
#' @param v A numeric value indicating the vertical distance of the square from the origin.
#' @param n An integer indicating the number of points to generate within the square.
#' 
#' @return An object of class "sf" containing n points within a square.
#' 
#' @importFrom sf st_polygon st_sample
#' @export
createPolyPoints <- function(h,v=0,n){
    square <- st_polygon(list(rbind(c(-h,v), c(v,h), c(h,v), c(v,-h), c(-h,v))))
    st_sample(square,size = n,type="hexagonal")
}

#' Create n points within two circles with a distance of m.
#' 
#' @param m A numeric value indicating the distance between the two circles.
#' @param n An integer indicating the number of points to generate within the two circles.
#' 
#' @return An object of class "sf" containing n points within two circles with a distance of m.
#' 
#' @importFrom sf st_polygon st_sample st_difference
#' @export
createFringePoints <- function(m,n){
    # create circle polygon inside the square polygon
    center <- c(0, 0)
    radius <- 0.92
    theta <- seq(0, 2*pi, length.out = 50+1)[-1] # create n points evenly spaced around the circle
    circle_poly_coords <- center + radius * cbind(cos(theta), sin(theta))
    circle_poly1 <- st_polygon(list(rbind(circle_poly_coords, circle_poly_coords[1,])))
    radius <- radius-m 
    circle_poly_coords <- center + radius * cbind(cos(theta), sin(theta))
    circle_poly2 <- st_polygon(list(rbind(circle_poly_coords, circle_poly_coords[1,])))
    fringe=st_difference(circle_poly1,circle_poly2)
    st_sample(fringe,size = n,type="hexagonal")
}


#' Convert spatial points to an igraph graph.
#' 
#' @param points An object of class "sf" containing the spatial points.
#' @param contagious_period An integer indicating the duration (in time steps) that an individual is contagious.
#' 
#' @return An object of class "igraph" representing the graph.
#' 
#' @importFrom igraph graph.adjacency E V neighbors
#' @importFrom sf st_distance
#' @export
toGraph <- function(points,contagious_period=10){
    distances <- as.matrix(st_distance(points))
    g <- graph.adjacency(distances, mode = "undirected", weighted = TRUE)
    E(g)$width=log(1/E(g)$weight) #when drowing the graph we draw only the shorter edge ie the link between very close nodes
    V(g)$state=0
    V(g)$state[1]=contagious_period # the node "1" become contagious for `contagious_period` time step.
    g
}


#' the main function
#'
#' @param g A graph object representing the network
#' @param prob_diffuse The probability of the culture diffusing to a visited individual
#' @param nvisits The number of neighbors visited by each contagious individual in each time step
#' @param contagious_period The length of time that an individual remains contagious
#'
#' @return A graph object representing the network, with updated states and colors
#' @export
diffuse_culture <- function(g, prob_diffuse,nvisits=1,contagious_period=10) {

    #for all node in the graph of state > 0, which correspond to "contagious" individual we will see if they visit other indivdiual
    for (v in V(g)[V(g)$state>0]) {
        neighbors <- neighbors(g, v) #get all the neighbors of the node v
        probas=1/((E(g)[.from(v)]$weight+1)^4)    #this is where the probablilty of visiting someone is defined, very important. Here it is defined as depending on the invert of the squared distance. I added + 1 on the distance because the spatial geography of the model is very small and distance are ofen <1 so adding one simplify 
        visit <- sample(neighbors, 1,prob=probas)  #we sample the `nvisits` neighbors given these probability
        if (V(g)$state[visit]==0 && runif(1) < prob_diffuse) V(g)$state[visit] <- contagious_period
    }
  V(g)$state[V(g)$state==1]=-1 #here I decided that once individual are going to be not contagious again I put them in a "resistant" stat, so they won't be contaminated again. We could imagine this resistant state to be also -n and then at each time state we can increase this resistance ; when it reaches 0 they can be contaminated again 
  V(g)$state[V(g)$state>0]=V(g)$state[V(g)$state>0]-1 #decrease the time left for contagious people
  V(g)$color[V(g)$state==0]=1 #just some estetic
  V(g)$color[V(g)$state>0]=3 #just some estetic
  V(g)$color[V(g)$state<0]=2 #just some estetic
  return(g)
}

n=2000

#let's draw 3 maps

scenarios <- list( createFringePoints(.3,n), createPolyPoints(h=1.5,n=n), createPolyPoints(h=.5,n=n))

par(mfrow=c(1,3))
for(i in 1:3){
    plot(0,0,xlim=c(-1,1),ylim=c(-1,1),type="n",ann=F,axes=F)
    plot(scenarios[[i]],add=T,pch=21,bg="red")
}


# convert these maps as graph where we will run the simulation:
graphs=lapply(scenarios,toGraph)

# Run diffusion simulation for 25 time steps on each graph
prob_diffuse <- .5
par(mfrow=c(1,3))
par(mar=c(0,0,0,0))
for (t in 1:500) {
    st=Sys.time()
    #update 3 scenarios
    graphs=lapply(graphs, diffuse_culture,prob_diffuse=prob_diffuse)
    print(Sys.time()-st)
    #plots 3 scenarios
    lapply(1:length(graphs),function(i){
               plotgraph(graphs[[i]],scenarios[[i]]))
})
    #Sys.sleep(1)
}

plotgraph <- function(network,layout,...){
    plot(0,0,xlim=c(-1.5,1.5),ylim=c(-1.5,1.5),type="n",ann=F,axes=F)
    st=V(network)$state
    cl=st
    cl[st==0]=0
    cl[st>0]=2
    cl[st<0]="blue"
    plot(layout,bg=cl,add=T,pch=21,...)
}

#Final result after 25 time step:
pdf("network.pdf",width=24,height=8)

par(mfrow=c(1,3))
par(mar=c(0,0,0,0))
for(i in 1:3){
    plot(0,0,xlim=c(-1.5,1.5),ylim=c(-1.5,1.5),type="n",ann=F,axes=F)
    plot(graphs[[i]],layout=st_coordinates(scenarios[[i]]),add=T,rescale=F)
}
dev.off()



#reset the state of the graph:


# below example to run in parallel multiple runs


#run multiple simulation
library(parallel)


for(n in c(50,100)){
    for(prob_diffuse in c(.2,.3)){
        scenarios <- list( createFringePoints(.3,n), createPolyPoints(h=1.5,n=n), createPolyPoints(h=.5,n=n))
        graphs=lapply(scenarios,toGraph)

        st=Sys.time()
        #below the while loop run until only 2 individuals are still contagious and return the number of time step neede to do so.
        cl <- makeCluster(4)
        bigg=lapply(seq_along(graphs),function(g)
                    { 
                        parLapply(cl,1:10,function(i,g,graphs,n,prob_diffuse,diffuse_culture){

                                      library(igraph)
                                      #lapply(1:10,function(i){
                                      print(i);
                                      gt=graphs[[g]]
                                      t=0;
                                      crve=c();
                                      while( sum(V(gt)$state<0)<.99*n && t < 1500 && sum(V(gt)$state>0)>0 ){
                                          #if(i%%50==1){
                                          #    png(sprintf("N%d_pd%d_r%d_layout%d_t%03d_b.png",n,prob_diffuse*10,i,g,t),width=700,height=700,pointsize=10)
                                          #    par(mar=c(0,0,0,0))

                                          # plotgraph(gt,scenarios[[g]],lwd=.6)
                                          #}
                                          gt=diffuse_culture(gt,prob_diffuse);
                                          t=t+1;
                                          ni=sum(V(gt)$state>0); 
                                          crve=c(crve,ni)
                                          #if(i%%50==1){
                                          #    dev.off()
                                          #}
                                      };
                                      print(crve)
                                      return(list(t,crve)) 
               },graphs=graphs,g=g,n=n,prob_diffuse=prob_diffuse,diffuse_culture=diffuse_culture)
                    })
        stopCluster(cl)
        saveRDS(file=paste0("result_gis_n",n,"_p",prob_diffuse,".rds"),bigg)
        print(Sys.time()-st)
    }
}


boxplot(t(bigg))
