if (F) {
  #f <- function(xx){apply(xx,1,function(yy){exp(-.005*yy[1]-.5*(yy[2]+.03*yy[1]^2-3)^2)})}
  f <- function(xx){
    if(is.matrix(xx)) {
      return(apply(xx,1,function(yy){exp(-.005*(yy[1]*40-20)^2-.5*(yy[2]*15-10+.03*(yy[1]*40-20)^2-3)^2)}))
    } else {
      yy <- xx
      exp(-.005*(yy[1]*40-20)^2-.5*(yy[2]*15-10+.03*(yy[1]*40-20)^2-3)^2)
    }
  }
  
  # Get contour plot
  nc <- 100 # number contour each dimension
  fx <- fy <- seq(0,1,length.out = nc)
  fz <- matrix(0,nc,nc)
  for (xi in 1:nc) for(yi in 1:nc) fz[xi,yi] <- f(matrix(c(fx[xi],fy[yi]),1,2))
  contour(fx,fy,fz,nlevels = 5)
  my.filled.contour(fx,fy,fz,nlevels = 5)
  
  p <- 2 # dimension
  k <- 4*p # MED distance thing
  
  qq <- function(xx){f(xx)^-(1/(2*p))}
  f_min <- function(xnew,xall,kk) {
    #print(c())
    #print(c(qq(xnew)^kk,-12,(apply(xall,1,function(xx){(qq(xx)/(sqrt(sum((xx-xnew)^2))))^kk}))))
    #qq(xnew)^k*sum(sapply(xall,function(xx){(qq(xx)/abs(xx-xnew))^k}))
    qq(xnew)^kk*sum(apply(xall,1,function(xx){(qq(xx)/(sqrt(sum((xx-xnew)^2))))^kk}))
  }
  
  
  # Initialize with mode
  gsa.out <- GenSA::GenSA(par=NULL,fn=function(xx)-f(xx),lower=c(0,0),upper=c(1,1),control = list(trace.mat=F))
  X <- matrix(gsa.out$par,1,2)
  text(X[1],X[2],labels=1,col=1,pch=1)
  for(i in 2:30) {
    # GenSA is slow, maybe change parameters to limit time
    #gsa.out <- GenSA::GenSA(par=NULL,fn=function(xx)log(f_min(xx,X,kk=k)),lower=c(0,0),upper=c(1,1),control = list(maxit=10,trace.mat=F,max.time=1,max.call=10))
    
    # Or use default
    gsa.out <- GenSA::GenSA(par=NULL,fn=function(xx)log(f_min(xx,X,kk=k)),lower=c(0,0),upper=c(1,1),control = list(trace.mat=F,max.time=1))
    
    # Add new point
    xnew <- gsa.out$par
    X <- rbind(X,xnew)
    text(x=xnew[1],y=xnew[2],labels=i,col=1)
  }
}

SMED_GP_2D <- function(f,n0=10,n=10,nc=100,max.time=NULL,GP.package='GPfit') {
  p <- 2 # dimension
  k <- 4*p # MED distance thing
  GenSA.controls <- list(trace.mat=F,max.call=5) 
  if(!is.null(max.time)) GenSA.controls[['max.time']] <- max.time
  
  
  # source('TestFunctions.R')
  source('myfilledcontour.R')
  require('lhs')
  
  if(GP.package=='GPfit') {
    require('GPfit')
    predict.GP.SMED <- function(mod,xx){
      GPfit::predict.GP(mod,matrix(xx,1,2))$Y_hat
    }
    #fit.GP.SMED <- GPfit::GP_fit
    init.GP.SMED <- GPfit::GP_fit
    update.GP.SMED <- function(mod,X,Y){GPfit::GP_fit(X,Y)}
  } else if(GP.package=='laGP') {
    init.GP.SMED <- function(X,Y) {laGP::newGPsep(X=X,Z=Y,d=2,g=1e-6)}
    #fit.GP.SMED <- laGP::upGPsep
    update.GP.SMED <- function(mod,X,Y) {laGP::updateGPsep(gpsepi=mod,X=X,Z=Y);return(mod)}
    predict.GP.SMED <- function(mod,xx){yy <- laGP::predGPsep(mod,matrix(xx,1,2))$mean;print(yy);return(yy) }
    delete.GP.SMED <- laGP::deleteGPsep
  } else if(GP.package=='mlegp') {
    require('mlegp')
    predict.GP.SMED <- function(mod,xx) {mlegp::predict.gp(object=mod,newData=xx)}
    fit.GP.SMED <- function(X,Y) {mlegp::mlegp(X=X,Z=Y,verbose=0)}
    init.GP.SMED <- function(X,Y) {mlegp::mlegp(X=X,Z=Y,verbose=0)}
    update.GP.SMED <- function(mod,X,Y) {mlegp::mlegp(X=X,Z=Y,verbose=0)}
  } else {
    stop('No GP.package specified 579238572')
  }
  
  
  # Charge function qq
  #qq <- function(xx){f(xx)^-(1/(2*p))}
  qq <- function(xx,mod){(predict.GP.SMED(mod=mod,xx=xx))^-(1/(2*p))}
  # Function we will optimize
  f_min <- function(xnew,xall,kk,mod) {
    qq(xnew,mod=mod)^kk*sum(apply(xall,1,function(xx){(qq(xx,mod=mod)/(sqrt(sum((xx-xnew)^2))))^kk}))
  }
  
  # Get contour plot
  fx <- fy <- seq(0,1,length.out = nc)
  fz <- matrix(0,nc,nc)
  for (xi in 1:nc) for(yi in 1:nc) fz[xi,yi] <- f(matrix(c(fx[xi],fy[yi]),1,2))
  #contour(fx,fy,fz,nlevels = 5) # This had legend, don't want it
  my.filled.contour(fx,fy,fz,nlevels = 5)
  
  # Initialize with LHS
  X <- lhs::maximinLHS(n=n0,k=2)
  Y <- apply(X,1,f)
  text(X[,1],X[,2],labels=1:n0,col=2,pch=1)
  
  browser()
  # Get rest of points
  for(i in 1:n) {
    print(paste('Starting',i))
    
    # Create GP model
    #mod <- fit.GP.SMED(X=X, Y=Y)
    if(i==1) mod <- init.GP.SMED(X,Y)
    else mod <- update.GP.SMED(mod,X,Y)
    my.filled.contour.func(function(xx)predict.GP.SMED(mod,xx))
    text(x=X[,1],y=X[,2],col='magenta')
    print('   fit')
    # Use log scale for optimization, had trouble before when numbers were 1e88
    #gsa.out <- GenSA::GenSA(par=NULL,fn=function(xx)log(f_min(xx,X,kk=k,mod=mod)),lower=c(0,0),upper=c(1,1),control = GenSA.controls)
    gsa.fn.call.count <- 0
    #gsa.out <- GenSA::GenSA(par=NULL,
    #                        fn=function(xx){ gsa.fn.call.count <<- gsa.fn.call.count + 1
    #                                         print(gsa.fn.call.count)
    #                          points(xx[1],xx[2],pch=19,cex=.1);
    #                          if(gsa.fn.call.count>100) return(Inf)#stop('GenSA no good 54829357')
    #                          log(f_min(xx,X,kk=k,mod=mod))
    #                          },
    #                        lower=c(0,0),upper=c(1,1),control = GenSA.controls)
    #gsa.out <- optim(par=runif(2),fn=function(xx)log(f_min(xx,X,kk=k,mod=mod)),lower=c(0,0),upper=c(1,1),
    #                 method='L-BFGS-B')
    gsa.out <- rgenoud::genoud(fn=function(xx){points(xx[1],xx[2],pch=19,cex=.1);log(f_min(xx,X,kk=k,mod=mod))}
                    ,nvars=2,max.generations=1,hard.generation.limit=T
                    ,Domains=matrix(c(0,0,1,1),2,2,byrow=F),boundary.enforcement=T
                    )
    # Add new point
    xnew <- gsa.out$par
    print(xnew)
    X <- rbind(X,xnew)
    Y <- c(Y,f(xnew))
    text(x=xnew[1],y=xnew[2],labels=i,col=1)
    browser()
  }
  # Return design matrix
  rownames(X) <- NULL
  delete.GP.SMED(mod)
  return(X)
}
if (F) {
  setwd("C:/Users/cbe117/School/DOE/SMED/SMED-Code")
  source('TestFunctions.R')
  SMED_GP_2D(f=banana,n0=10,n=10,max.time=.2)
}