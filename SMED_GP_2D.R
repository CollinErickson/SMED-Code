SMED_GP_2D <- function(f,n0=10,n=10,nc=100,GP.package='',opt.method='genoud',contour.fit=0,continue.option=F) {
  # Function that implements SMED with GP in 2 dimensions. Using the GP means you 
  #  only get the function value from points you have selected, selection is 
  #  based on the GP model predictions
  # Input:
  #  n0: # of pts to start with, selected with maximin LHS
  #  n: # of pts to select -> total pts is n0+n
  #  nc: grid size of contour plots are nc by nc
  #  GP.package: which R package to fit GP with. Current options are 
  #    GPfit, laGP, and mlegp. Also exact will use exact function (for debugging)
  #  opt.method: optimization method. Current options are 'genoud', 'LHS'.
  #  contour.fit: # of iterations after which to plot model contour plot
  #  continue.option: If TRUE, you are given choice to continue sampling
  # Output: 
  #  X: design matrix
  #  Y: function values at design points
  
  p <- 2 # dimension
  k <- 4 * p # MED distance power, k=4p is recommended by Roshan
  
  # source('TestFunctions.R') # Not sure if I should source test functions here or outside
  source('myfilledcontour.R') # Source myfilledcontour.R to get contour plots
  require('lhs') # Use lhs to get initial points
  
  # Set GP functions for each package
  #  Need functions for init, update, predict, and delete
  if(GP.package=='GPfit') {
    require('GPfit')
    predict.GP.SMED <- function(mod,xx){max(1e-16,GPfit::predict.GP(mod,matrix(xx,1,2))$Y_hat)}
    init.GP.SMED <- GPfit::GP_fit
    update.GP.SMED <- function(mod,X,Y){GPfit::GP_fit(X,Y)}
    delete.GP.SMED <- function(mod){}
  } else if(GP.package=='laGP') {
    require('laGP')
    init.GP.SMED <- function(X,Y) {laGP::newGPsep(X=X,Z=Y,d=2,g=1e-8)}
    update.GP.SMED <- function(mod,X,Y) {laGP::updateGPsep(gpsepi=mod,X=X,Z=Y);return(mod)}
    predict.GP.SMED <- function(mod,xx){max(1e-16,laGP::predGPsep(mod,matrix(xx,1,2))$mean)}
      # Changed this predict, really bad now
    delete.GP.SMED <- laGP::deleteGPsep
  } else if(GP.package=='mlegp') {
    require('mlegp')
    predict.GP.SMED <- function(mod,xx) {max(1e-16,mlegp::predict.gp(object=mod,newData=xx))}
    init.GP.SMED <- function(X,Y) {mlegp::mlegp(X=X,Z=Y,verbose=0)}
    update.GP.SMED <- function(mod,X,Y) {mlegp::mlegp(X=X,Z=Y,verbose=0)}
    delete.GP.SMED <- function(mod){}
  } else if (GP.package=='exact') {
    predict.GP.SMED <- function(mod,xx) {f(xx)}
    init.GP.SMED <- function(X,Y) {}
    update.GP.SMED <- function(mod,X,Y) {}
    delete.GP.SMED <- function(mod){}    
  } else {
    stop('No GP.package specified Error # 579238572')
  }
  
  
  # Charge function qq, no longer is actual function, now is based on predictions
  qq <- function(xx,mod){
    yy <- (predict.GP.SMED(mod=mod,xx=xx))^-(1/(2*p));
    if(!is.nan(yy)){
      return(yy)
    } else { 
      browser();stop(paste('NaN qq value Error 357923',predict.GP.SMED(mod=mod,xx=xx)))
    }
  }
  # Function we will optimize
  f_min <- function(xnew,xall,kk,mod,qq.scale=1) {
    (qq(xnew,mod=mod)*qq.scale)^kk*sum(apply(xall,1,function(xx){((qq(xx,mod=mod)*qq.scale)/(sqrt(sum((xx-xnew)^2))))^kk}))
  }
  # Log version
  log.f_min <- function(xnew,xall,kk,mod,qq.scale=1) {
    kk * log(qq(xnew,mod=mod)*qq.scale) + log(sum(apply(xall,1,function(xx){((qq(xx,mod=mod)*qq.scale)/(sqrt(sum((xx-xnew)^2))))^kk})))
  }
  
  
  # Get contour plot
  my.filled.contour.func(fn=f,n=nc)
  
  # Initialize with LHS
  X <- lhs::maximinLHS(n=n0,k=2)
  Y <- apply(X,1,f)
  Delta <- .01 * max(Y)
  keep.Delta <- (Y > Delta)
  # Plot where points are
  text(X[,1],X[,2],col=ifelse(keep.Delta,'magenta','orange'))
  
  # Initialize model
  mod <- init.GP.SMED(X,Y)
  
  # Get rest of points, while loop lets you add points in when prompted
  i <- 1
  while(i <= n) {
    # Plot contour of fit if iteration is multiple of contour.fit
    if(contour.fit>0 & i%%contour.fit==0) {
      my.filled.contour.func(function(xx)predict.GP.SMED(mod,xx))
      text(x=X[,1],y=X[,2],col=ifelse(keep.Delta,'magenta','springgreen3'))
    }
    # Perform optimzation to select next point, use log scale
    # Different options for optimization, genoud is default
    if (opt.method=='genoud') { # First optimization option is genoud
      opt.out <- rgenoud::genoud(fn=function(xx){points(xx[1],xx[2],pch=19,cex=.1,col=3);log(f_min(xx,X[keep.Delta,],kk=k,mod=mod))}
                      ,nvars=2,max.generations=3,hard.generation.limit=T
                      ,Domains=matrix(c(0,0,1,1),2,2,byrow=F),boundary.enforcement=T,pop.size=50
                      ,print.level=0
                      )
    } else if (opt.method == 'LHS') { # Optimization method: check a bunch of LHS points, NOT GOOD
      #browser()
      XX.LHS <- lhs::maximinLHS(n=1000,k=2)
      points(XX.LHS[,1],XX.LHS[,2],pch=19,cex=.1,col='orange')
      log.f_min.LHS <- apply(XX.LHS,1,function(xx){log(f_min(xx,X[keep.Delta,],kk=k,mod=mod))})
      points(XX.LHS[,1],XX.LHS[,2],pch=19,cex=.1,col=rgb(1,(1:100)/150,(1:100)/150)[floor(1+99*((log.f_min.LHS-min(log.f_min.LHS))/(max(log.f_min.LHS)-min(log.f_min.LHS))))])
      min.ind.LHS <- which.min(log.f_min.LHS)[1]
      opt.out <- list(par=XX.LHS[min.ind.LHS,],value=log.f_min.LHS[min.ind.LHS])
    } else {
      stop('opt.method not specified. Error #32572938')
    }
    # Add new point
    xnew <- opt.out$par
    ynew <- f(xnew)
    cat(n0+i,xnew,ynew,opt.out$val,log(f_min(xnew,X[keep.Delta,],kk=k,mod=mod)),log(f_min(runif(2),X[keep.Delta,],kk=k,mod=mod)),'\n')
    X <- rbind(X,unname(xnew))
    Y <- c(Y,ynew)
    text(x=xnew[1],y=xnew[2],labels=n0+i,col=1) # Plot it on contour
    
    # Delta is used to 'hide' points from the potential function that are not useful
    # Without this the points selected were terrible
    Delta <- .01 * (n0 / (n + i)) ^ (1 / p) * max(Y)
    keep.Delta <- (Y > Delta)
    
    # Update model
    mod <- update.GP.SMED(mod,X,Y)

    # If told to it will ask how many more points you want, negative integers are for debugging
    if(continue.option==T & i==n) {
      n.more <- as.integer(readline(prompt = 'How many more? '))
      if(n.more==-1) {browser();n.more <- as.integer(readline(prompt = 'How many more? '))}
      if(n.more==-2) {browser();my.filled.contour.func(function(xx){log(f_min(xx,X[keep.Delta,],kk=k,mod=mod))},n=30);n.more <- as.integer(readline(prompt = 'How many more? '))}
      if(is.integer(n.more) & length(n.more)==1 & !is.na(n.more)) {n <- n + n.more}
    }
    i <- i + 1
  }
  # Print final contour of estimated surface
  my.filled.contour.func(function(xx)predict.GP.SMED(mod,xx))
  text(x=X[,1],y=X[,2],col='magenta')
  # Delete model if needed (only laGP)
  delete.GP.SMED(mod)
  # Return design matrix and outputs
  return(cbind(X,Y))
}
if (F) {
  setwd("C:/Users/cbe117/School/DOE/SMED/SMED-Code")
  source('TestFunctions.R')
  SMED_GP_2D(f=banana,n0=10,n=30,contour.fit=1,GP.package='mlegp')
  SMED_GP_2D(f=banana,n0=50,n=5,contour.fit=1,GP.package='laGP',continue.option=T)
  SMED_GP_2D(f=function(x){exp(-(sum((x-.5)^2))/.01)},n0=10,n=5,contour.fit=1,GP.package='laGP',continue.option=T)
  SMED_GP_2D(f=banana,n0=10,n=5,contour.fit=1,GP.package='exact',continue.option=T)
  
  # Comparing with seed, pick one point
  set.seed(0);SMED_GP_2D(f=banana,n0=30,n=1,contour.fit=1,GP.package='GPfit',continue.option=T,opt.method='LHS')
  # 31 0.995832 0.9755436 9.061522e-41 18.48422 18.48422 22.16811 
  set.seed(0);SMED_GP_2D(f=banana,n0=30,n=1,contour.fit=1,GP.package='mlegp',continue.option=T,opt.method='LHS')
  # 31 0.9891865 0.986211 7.488378e-40 18.44477 18.44477 32.09332 
  set.seed(0);SMED_GP_2D(f=banana,n0=30,n=1,contour.fit=1,GP.package='laGP',continue.option=T,opt.method='LHS')
  # 31 0.6320691 0.9985909 0.01649461 17.28631 17.28631 89.38975  Much worse fit though
  set.seed(0);SMED_GP_2D(f=banana,n0=30,n=1,contour.fit=1,GP.package='exact',continue.option=T,opt.method='LHS')
  # 31 0.5693355 0.913252 0.6247042 20.49745 20.49745 118.301 
}