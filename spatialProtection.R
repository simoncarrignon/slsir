library(sf)
library(igraph)

#create points within a square
createPolyPoints <- function(h,v=0,n){
    square <- st_polygon(list(rbind(c(-h,v), c(v,h), c(h,v), c(v,-h), c(-h,v))))
    st_sample(square,size = n,type="hexagonal")
}

#create points within two circle
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
    st_sample(fringe,size = 50,type="hexagonal")
}


toGraph <- function(points){
    distances <- as.matrix(st_distance(points))
    g <- graph.adjacency(distances, mode = "undirected", weighted = TRUE)
    E(g)$width=log(1/E(g)$weight)
    E(g)$width=log(1/E(g)$weight)
    V(g)$state=0
    V(g)$state[1]=10
    g
}
diffuse_culture <- function(g, prob_diffuse) {
    for (v in V(g)[V(g)$state>0]) {
        neighbors <- neighbors(g, v)
        probas=1/(E(g)[.from(v)]$weight+1)^2
        visit <- sample(neighbors, 1,prob=probas)
        #message(paste("probas:",length(probas),"neighbors:", length(neighbors)))
        #message(paste("visit:",visit))
        if (V(g)$state[visit]==0 && runif(1) < prob_diffuse) V(g)$state[visit] <- 10
        }
  V(g)$state[V(g)$state==1]=-1
  V(g)$state[V(g)$state>0]=V(g)$state[V(g)$state>0]-1
  V(g)$color=ifelse(V(g)$state==0,1,2)
  return(g)
}

n=30
plot(0,0,xlim=c(-1,1),ylim=c(-1,1),type="n",ann=F,axes=F)
plot(createFringePoints(.3,n),add=T,pch=21,bg="red")
plot(createPolyPoints(h=1,n=n),add=T,pch=21,bg="red")
plot(createPolyPoints(h=.5,n=n),add=T,pch=21,bg="red")

scenarios <- list( createFringePoints(.3,n), createPolyPoints(h=1.5,n=n), createPolyPoints(h=.5,n=n))

n=30
plot(0,0,xlim=c(-1,1),ylim=c(-1,1),type="n",ann=F,axes=F)
plot(createFringePoints(.3,n),add=T,pch=21,bg="red")
plot(createPolyPoints(h=1,n=n),add=T,pch=21,bg="red")
plot(createPolyPoints(h=.5,n=n),add=T,pch=21,bg="red")

scenarios <- list( createFringePoints(.3,n), createPolyPoints(h=1.5,n=n), createPolyPoints(h=.5,n=n))


# Simulate cultural diffusion
graphs=lapply(scenarios,toGraph)
# Run diffusion simulation for 10 time steps
prob_diffuse <- .5
par(mfrow=c(1,3))
for (i in 1:25) {
    #update 3 scenarios
    graphs=lapply(graphs, diffuse_culture,prob_diffuse=prob_diffuse)
    #plots 3 scenarios
    lapply(1:length(graphs),function(i){
               plot(0,0,xlim=c(-1.5,1.5),ylim=c(-1.5,1.5),type="n",ann=F,axes=F)
               plot(graphs[[i]],layout=st_coordinates(scenarios[[i]]),add=T,rescale=F)
                })
               Sys.sleep(1)
}
plot(graphs[[1]],layout=st_coordinates(scenarios[[1]]))
plot(g1,layout=st_coordinates(scenarios[[1]]))

a=Sys.time()
graphs=lapply(scenarios,toGraph)
#run multiple simulation
testg=replicate(100,sapply(graphs,function(gt){ t=0; while(sum(V(gt)$state>2) ){gt=diffuse_culture(gt,prob_diffuse);t=t+1}; return(t) }))
Sys.time()-a
library(parallel)
cl <- makeCluster(10,type="FORK")
testgd=parSapply(cl,1:100,function(i){print(i);sapply(graphs,function(gt){ t=0; while(sum(V(gt)$state>2) ){gt=diffuse_culture(gt,prob_diffuse);t=t+1}; return(t) })})
stopCluster(cl)

cl <- makeCluster(10,type="FORK")
bigg=parSapply(cl,1:3000,function(i){print(i);sapply(graphs,function(gt){ t=0; while(sum(V(gt)$state>2) ){gt=diffuse_culture(gt,prob_diffuse);t=t+1}; return(t) })})
stopCluster(cl)

testgraphspeadd=cbind(testgraphspeadd,testgd3,testgd4)



 j  graphs=lapply(graphs, diffuse_culture,prob_diffuse=prob_diffuse)
