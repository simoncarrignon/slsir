#a bunch of tools
library(RColorBrewer)

bdistance <- function(x)(x-min(x))/(max(x)-min(x))
bscore <- function(x)(x-min(x))/(max(x)-min(x))

color_class=rev(brewer.pal(5,"PuBu"))
colorbest=color_class[1]
myred=rgb(r=213,g=94,b=0,maxColorValue=255)
myred="red"
mygreen=rgb(r=0,g=158,b=115,maxColorValue=255)
mygreen=colorbest


marginAndJoin  <- function(
                           distribA=list(dim1=posterior$worst$inf,dim2=posterior$worst$sat),
                           distribB=list(dim1=posterior$best$inf,dim2=posterior$best$sat),
                           expnames=c(distribA="Worst",distribB="Best"),
                           dimlab=list(dim1=expression(nu),dim2=expression(kappa)),
                           cols,
                           probas=c(75,90),
                           log="",
                           main="",
                           range=list(dim1=c(0,1),dim2=c(-1,3))
                           )
{

    #dimmnames=c(dim1="inf",dim2="sat")
    names(cols)=c("distribA","distribB")

    xs=list(distribA=distribA$dim1,distribB=distribB$dim1)
    ys=list(distribA=distribA$dim2,distribB=distribB$dim2)

    ## logarithmique transformation if need (if log != "") 
    if(length(grep("x",log))>0){
        xs=lapply(xs,log10)
    }
    if(length(grep("y",log))>0)
        ys=lapply(ys,log10)

    par(oma=c(3,3,0,0))

    layout(matrix(nrow=2,ncol=2,c(1,3,4,2)),widths=c(.7,.3),heights=c(.3,.7))

    par(mar=rep(1,4))
    plot.new()
    ds=lapply(xs,function(u)density(u,from=range$dim1[1],to=range$dim1[2]))
    plot.window(range$dim1,range(sapply(ds,"[[","y")), "", xaxs = "i",asp=NA)
    for(d in names(ds)){
        polygon(c(range$dim1[1],ds[[d]]$x,range$dim1[2]),c(0,ds[[d]]$y,0),col=cols[d],lwd=2)
    }
    #mtext(lab$dim1,3,1)

    plot.new()
    ds=lapply(ys,function(u)density(u,from=range$dim2[1],to=range$dim2[2]))
    plot.window(range(sapply(ds,"[[","y")),range$dim2, "", yaxs = "i",asp=NA)
    for(d in names(ds)){
        polygon(c(0,ds[[d]]$y,0),c(range$dim2[1],ds[[d]]$x,range$dim2[2]),col=cols[d],lwd=2)
    }
    #mtext(lab$dim2,4,1)


    for(s in names(cols)){
        x=xs[[s]]
        y=ys[[s]]
        hdr.boxplot.2d(x=x,y=y,prob=probas,shadecols=cols[s],xlim=range$dim1,ylim=range$dim2,xlab=dimlab$dim1,ylab=dimlab$dim2,yaxt="n",outside.points=F,yaxt="n",axes=F,frame.plot=F)
        par(new=T)
    }
    box()

    if(length(grep("x",log))>0){
        xaxp=par()$xaxp
        axis(1,at=seq(xaxp[1],xaxp[2],length.out=xaxp[3]+1),label=10^seq(xaxp[1],xaxp[2],length.out=xaxp[3]+1))
    }else{axis(1)}
    mtext(dimlab$dim2,2,3)
    mtext(dimlab$dim1,1,3)
    if(length(grep("y",log))>0){
        yaxp=par()$yaxp
        axis(2,at=seq(yaxp[1],yaxp[2],length.out=yaxp[3]+1),label=10^seq(yaxp[1],yaxp[2],length.out=yaxp[3]+1))
    }else{axis(2)}

    plot.new()
    par(xpd=NA)
    legend("center",legend=expnames,fill=cols,bty="n")
    mtext(main,3,0,out=T)
    par(xpd=F)

}

#' 
#' 
#' create a gradient of colors givent a vector of value
#' 
#' @param x a vector of value to be matched with colors
#' @param colors a vector of three colors
#' @param colsteps the number of cilors groups
color.gradient <- function(x, colors=c(myred,"yellow",mygreen), colsteps=100) {
    return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
}

#' plot posteriors distribution against priors
#'@param A : a vector with posterior
#'@param B : a vector with posterior 
#'@param prior : a vector with posterior 
#' @export plot2dens 
plot2dens <- function(A=NULL,B=NULL,C=NULL,from=NULL,to=NULL,prior=NULL,cols=c(adjustcolor("red",.8),adjustcolor("blue",.8),adjustcolor("yellow",.8)),hdr=F,yaxt=NULL,log=F,...){

    denseP=NULL
    denseA=NULL
    denseB=NULL
    denseC=NULL
    if(!is.null(prior))prior=prior[!is.na(prior)]
    if(is.null(yaxt))yaxt="n"
    if(is.null(from))from=min(A,B,prior)
    if(is.null(to))to=max(A,B,prior)
    if(!is.null(A))denseA=density(A,from=from,to=to)
    if(!is.null(B))denseB=density(B,from=from,to=to)
    if(!is.null(C))denseC=density(C,from=from,to=to)
    if(length(prior)==2)denseP=density(runif(100000,prior[1],prior[2]),from=from,to=to)
    else if(!is.null(prior))denseP=density(prior,from=from,to=to)

    if(is.null(names(cols)))names(cols)=c("P","A","B")
    rangex=range(denseB$x,denseA$x,denseP$x,denseC$x)
    rangey=range(0,denseB$y,denseA$y,denseP$y,denseC$y)
    stepy=max(rangey)*0.2
    if(hdr)miny=-1*stepy else miny=0
    plot(denseA,ylim=c(miny,max(rangey)),xlim=rangex,type="n",xaxt="n",yaxt=yaxt,ylab="",...)
    if(log){
        xaxp=par()$xaxp
        axis(1,at=seq(xaxp[1],xaxp[2],length.out=xaxp[3]+1),label=10^seq(xaxp[1],xaxp[2],length.out=xaxp[3]+1))
    }
    else axis(1)
    mtext("Density",2,1)
    if(!is.null(C))
        polygon(c(from,denseC$x,to),c(0,denseC$y,0),col=cols["C"],lwd=2)
    if(!is.null(A)){
        #polygon(c(from,denseA$x,to),c(0,denseA$y,0),col="white",border=NA,lwd=2)
        if(hdr){
            hdstaA=hdr(A,prob=c(75,95),lambda=0.9)
            hdrA=hdstaA$hdr
            print(hdrA)
            polygon(c(hdrA[1,1],hdrA[1,1],hdrA[1,2],hdrA[1,2]),c(-.2*stepy,-.9*stepy,-.9*stepy,-.2*stepy),col=cols["A"],lwd=1)

            polygon(c(hdrA[2,1],hdrA[2,1],hdrA[2,2],hdrA[2,2]),c(-.2*stepy,-.9*stepy,-.9*stepy,-.2*stepy),col=cols["A"],lwd=1)
            polygon(c(hdrA[2,1],hdrA[2,1],hdrA[2,2],hdrA[2,2]),c(-.2*stepy,-.9*stepy,-.9*stepy,-.2*stepy),col=cols["A"],lwd=1)
            segments(hdstaA$mode,-.3*stepy,hdstaA$mode,-.8*stepy)
            middle=(-.2*stepy+-.9*stepy)/2
            segments(min(A),middle,hdrA[1,1],middle)
            segments(hdrA[1,2],middle,max(A),middle)
            polygon(c(from,denseA$x,to),c(0,denseA$y,0),col=cols["A"],lwd=2)
        }
        else
            polygon(c(from,denseA$x,to),c(0,denseA$y,0),col=cols["A"],lwd=2)
    }
    if(!is.null(B)){
        polygon(c(from,denseB$x,to),c(0,denseB$y,0),col=adjustcolor("white",.6),lwd=2,density=NA,border=0)
        polygon(c(from,denseB$x,to),c(0,denseB$y,0),col=cols["B"],lwd=2,density=5,border=1)
    }
    if(!is.null(prior))
        polygon(c(from,denseP$x,to),c(0,denseP$y,0),col=cols["P"],lwd=2)

    #if(!is.null(A))
    #        polygon(c(from,denseA$x,to),c(0,denseA$y,0),col=NA,lwd=2)

}

