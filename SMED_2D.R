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

SMED_2D <- function(f,n=10,nc=100,max.time=NULL) {
  # SMED in 2D
  #  f: function
  #  n: # of pts to select
  #  nc: # of pts in contour plot
  #  max.time: max.time for GenSA optimization for each point
  
  # source('TestFunctions.R')
  # source('myfilledcontour.R')
  
  p <- 2 # dimension
  k <- 4*p # MED distance thing
  GenSA.controls <- list(trace.mat=F)
  if(!is.null(max.time)) GenSA.controls[['max.time']] <- max.time
  
  # Charge function qq
  qq <- function(xx){f(xx)^-(1/(2*p))}
  # Function we will optimize
  f_min <- function(xnew,xall,kk) {
    qq(xnew)^kk*sum(apply(xall,1,function(xx){(qq(xx)/(sqrt(sum((xx-xnew)^2))))^kk}))
  }
  
  # Get contour plot
  fx <- fy <- seq(0,1,length.out = nc)
  fz <- matrix(0,nc,nc)
  for (xi in 1:nc) for(yi in 1:nc) fz[xi,yi] <- f(matrix(c(fx[xi],fy[yi]),1,2))
  #contour(fx,fy,fz,nlevels = 5) # This had legend, don't want it
  my.filled.contour(fx,fy,fz,nlevels = 5)
 
  # Initialize with mode
  gsa.out <- GenSA::GenSA(par=NULL,fn=function(xx)-f(xx),lower=c(0,0),upper=c(1,1),control = list(trace.mat=F))
  X <- matrix(gsa.out$par,1,2)
  text(X[1],X[2],labels=1,col=1,pch=1)
  
  # Get rest of points
  for(i in 2:n) {
    # Use log scale for optimization, had trouble before when numbers were 1e88
    gsa.out <- GenSA::GenSA(par=NULL,fn=function(xx)log(f_min(xx,X,kk=k)),lower=c(0,0),upper=c(1,1),control = GenSA.controls)
    # Add new point
    xnew <- gsa.out$par
    X <- rbind(X,xnew)
    text(x=xnew[1],y=xnew[2],labels=i,col=1)
  }
  # Return design matrix
  rownames(X) <- NULL
  return(X)
}
if (F) {
  SMED_2D(banana,n=10,max.time=.2)
}