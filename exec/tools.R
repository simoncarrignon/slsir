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
