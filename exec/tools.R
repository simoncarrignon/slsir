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
        print(cols[s])
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

duocolrs=alpha(c(myred,colorbest),.8)
currange=range=list(dim1=c(0,1),dim2=c(-1,3))
dimlab=list(dim1=expression(nu),dim2=expression(kappa))

### Cmpound marginal with 2d  learning

pdf("5e8779c56517890001536101/figures/2d_marginal_learning_best_vs_worst.pdf",pointsize=30,  width = 10, height = 10)
expnames=c(distribA="Worst",distribB="Best")
marginAndJoin(
              distribA=list(dim1=posterior$worst$inf,dim2=posterior$worst$sat),
              distribB=list(dim1=posterior$best$inf,dim2=posterior$best$sat),
              dimlab=list(dim1=expression(nu),dim2=expression(kappa)),
              expnames=expnames, range=currange,cols=duocolrs,log="y",
              main=""
              )
dev.off()


pdf("5e8779c56517890001536101/figures/2d_marginal_learning_best_vs_worst_150.pdf",pointsize=30,  width = 10, height = 10)
expnames=c(distribA="Worst 150",distribB="Best 150")
marginAndJoin(
              distribA=list(dim1=posterior$wmax_infect150$inf,dim2=posterior$wmax_infect150$sat),
              distribB=list(dim1=posterior$max_infect150$inf,dim2=posterior$max_infect150$sat),
              dimlab=list(dim1=expression(nu),dim2=expression(kappa)),
              expnames=expnames, range=currange,cols=duocolrs,log="y",
              main=""
              )
dev.off()


### COmpound marginal with 2d revert learning

pdf("5e8779c56517890001536101/figures/2d_marginal_revert_learning_best_vs_worst.pdf",pointsize=30,  width = 10, height = 10)
expnames=c(distribA="Worst",distribB="Best")
marginAndJoin(
              distribA=list(dim1=posterior$worst$inf_r,dim2=posterior$worst$sat_r),
              distribB=list(dim1=posterior$best$inf_r,dim2=posterior$best$sat_r),
              dimlab=list(dim1=expression(nu[r]),dim2=expression(kappa[r])),
              expnames=expnames, range=currange,cols=duocolrs,log="y",
              main=""
              )
dev.off()


pdf("5e8779c56517890001536101/figures/2d_marginal_revert_learning_best_vs_worst_150.pdf",pointsize=30,  width = 10, height = 10)
expnames=c(distribA="Worst 150",distribB="Best 150")
marginAndJoin(
              distribA=list(dim1=posterior$wmax_infect150$inf_r,dim2=posterior$wmax_infect150$sat_r),
              distribB=list(dim1=posterior$max_infect150$inf_r,dim2=posterior$max_infect150$sat_r),
              dimlab=list(dim1=expression(nu[r]),dim2=expression(kappa[r])),
              expnames=expnames, range=currange,cols=duocolrs,log="y",
              main=""
              )
dev.off()


################ best vs prior

### Cmpound marginal with 2d  learning

pdf("5e8779c56517890001536101/figures/2d_marginal_learning_best_vs_prior.pdf",pointsize=30,  width = 10, height = 10)
expnames=c(distribA="Prior",distribB="Best")
marginAndJoin(
              distribA=list(dim1=allresults$inf,dim2=allresults$sat),
              distribB=list(dim1=posterior$best$inf,dim2=posterior$best$sat),
              dimlab=list(dim1=expression(nu),dim2=expression(kappa)),
              expnames=expnames, range=currange,cols=duocolrs,log="y",
              main=""
              )
dev.off()


pdf("5e8779c56517890001536101/figures/2d_marginal_learning_best_vs_prior_150.pdf",pointsize=30,  width = 10, height = 10)
expnames=c(distribA="Prior",distribB="Best 150")
marginAndJoin(
              distribA=list(dim1=allresults$inf,dim2=allresults$sat),
              distribB=list(dim1=posterior$max_infect150$inf,dim2=posterior$max_infect150$sat),
              dimlab=list(dim1=expression(nu),dim2=expression(kappa)),
              expnames=expnames, range=currange,cols=duocolrs,log="y",
              main=""
              )
dev.off()


### COmpound marginal with 2d revert learning

pdf("5e8779c56517890001536101/figures/2d_marginal_revert_learning_best_vs_prior.pdf",pointsize=30,  width = 10, height = 10)
expnames=c(distribA="Prior",distribB="Best")
marginAndJoin(
              distribA=list(dim1=allresults$inf_r,dim2=allresults$sat_r),
              distribB=list(dim1=posterior$best$inf_r,dim2=posterior$best$sat_r),
              dimlab=list(dim1=expression(nu[r]),dim2=expression(kappa[r])),
              expnames=expnames, range=currange,cols=duocolrs,log="y",
              main=""
              )
dev.off()


pdf("5e8779c56517890001536101/figures/2d_marginal_revert_learning_best_vs_prior_150.pdf",pointsize=30,  width = 10, height = 10, fillOddEven = F)
expnames=c(distribA="Prior",distribB="Best 150")
marginAndJoin(
              distribA=list(dim1=allresults$inf_r,dim2=allresults$sat_r),
              distribB=list(dim1=posterior$max_infect150$inf_r,dim2=posterior$max_infect150$sat_r),
              dimlab=list(dim1=expression(nu[r]),dim2=expression(kappa[r])),
              expnames=expnames, range=currange,cols=duocolrs,log="y",
              main=""
              )
dev.off()


#### AGE FROM PAPER 
ages=c("< 24"=.19,
"25-44"=.45,
"45+"=.36)
