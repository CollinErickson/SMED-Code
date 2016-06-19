SMED <- function(f,p,n=10,nc=100,max.time=NULL) {
  # Function for SMED in 2D
  # Input:
  #  f: function
  #  p: # of dimensions
  #  n: # of pts to select
  #  nc: # of pts in contour plot
  #  max.time: max.time for GenSA optimization for each point
  
  # source('TestFunctions.R')
  # source('myfilledcontour.R')
  
  #p <- d # dimension
  k <- 4*p # MED distance thing
  GenSA.controls <- list(trace.mat=F) # Optimization parameters
  if(!is.null(max.time)) GenSA.controls[['max.time']] <- max.time
  
  # Charge function qq
  qq <- function(xx){f(xx)^-(1/(2*p))}
  # Function we will optimize
  f_min <- function(xnew,xall,kk) {
    qq(xnew)^kk*sum(apply(xall,1,function(xx){(qq(xx)/(sqrt(sum((xx-xnew)^2))))^kk}))
  }
  
  if (p==2) {
    # Get contour plot
    #my.filled.contour.func(f,nlevels=5)
    contourfilled.func(f,nlevels=5)
  }
  
  # Initialize with mode
  gsa.out <- GenSA::GenSA(par=NULL,fn=function(xx)-f(xx),lower=rep(0,p),upper=rep(1,p),control = list(trace.mat=F))
  X <- matrix(gsa.out$par,1,p)
  if (p==2) {
    text(X[1],X[2],labels=1,col=1,pch=1)
  }
  
  # Get rest of points
  for(i in 2:n) {
    # Use log scale for optimization, had trouble before when numbers were 1e88
    gsa.out <- GenSA::GenSA(par=NULL,fn=function(xx)log(f_min(xx,X,kk=k)),lower=rep(0,p),upper=rep(1,p),control = GenSA.controls)
    # Add new point
    xnew <- gsa.out$par
    X <- rbind(X,unname(xnew))
    if (p==2) {
      text(x=xnew[1],y=xnew[2],labels=i,col=1)
    }
  }
  if (p>2) {
    pairs(X)
  }
  # Return design matrix
  return(X)
}
if (F) {
  #setwd("C:/Users/cbe117/School/DOE/SMED/SMED-Code")
  source('TestFunctions.R')
  source('C:/Users/cbe117/School/DOE/Codes/contour/contourfilled/R/contourfilled.R')
  SMED(banana,p=2,n=10,max.time=.2)
  SMED(function(xx){xx[1]+xx[2]^2-sin(2*pi*xx[3])},p=3,n=10,max.time=.2)
}