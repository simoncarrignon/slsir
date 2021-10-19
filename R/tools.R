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
