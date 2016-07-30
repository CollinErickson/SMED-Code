SMED_select <- function(f,p,n=10,nc=100,max.time=NULL, X0=NULL, Xopt=NULL) {
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
  X <- X0
  n0 <- if(is.null(X0)) 0 else nrow(X0)
  
  # Charge function qq
  qq <- function(xx){f(xx)^-(1/(2*p))}
  # Function we will optimize
  #f_min <- function(xnew,xall,kk) {
  #  qq(xnew)^kk*sum(apply(xall,1,function(xx){(qq(xx)/(sqrt(sum((xx-xnew)^2))))^kk}))
  #}
  f_min <- function(xnew,xall, qqall,kk) {
    #qq(xnew)^kk*sum(apply(xall,1,function(xx){(qq(xx)/(sqrt(sum((xx-xnew)^2))))^kk}))
    qq(xnew)^kk*sum(sapply(1:nrow(xall),function(ii){(qqall[ii]/(sqrt(sum((xall[ii,]-xnew)^2))))^kk}))
  }
  
  if (p==2) {
    # Get contour plot
    #my.filled.contour.func(f,nlevels=5)
    #contourfilled.func(f,nlevels=5)
  }

  #points(X, pch=19)

  Y <- apply(X,1,f)
  qqX <- apply(X,1,qq)
  Delta <- .01 * max(Y)
  keep.Delta <- (Y > Delta)
  #browser()
  Xopt.inds <- 1:nrow(Xopt)
  Xopt.selected <- c()
  
  # Get rest of points
  for(i in 1:n) {#print(keep.Delta)
    # Use log scale for optimization, had trouble before when numbers were 1e88
    opt.func <- function(xx)log(f_min(xx,X[keep.Delta, , drop=F],qqX[keep.Delta],kk=k))
  
    func.vals <- apply(Xopt, 1, opt.func)
    best <- which.min(func.vals)
    xnew <- Xopt[best,]
    Xopt <- Xopt[-best, , drop=F]
    Xopt.selected <- c(Xopt.selected, Xopt.inds[best])
    Xopt.inds <- Xopt.inds[-best]

    X <- rbind(X,unname(xnew))
    Y <- apply(X,1,f) # could just append
    qqX <- apply(X,1,qq)
    # Delta is used to 'hide' points from the potential function that are not useful
    # Without this the points selected were terrible
    Delta <- .01 * max(Y) #.01 * (n0 / (n + i)) ^ (1 / p) * max(Y)
    #keep.Delta <- (Y > Delta)
    #browser()
    keep.Delta <- ifelse(1:nrow(X) <= n0, Y > Delta, T)
    
    if (p==2) {
      #text(x=xnew[1],y=xnew[2],labels=i,col=1)
    }
  }
  if (p>2) {
    #pairs(X)
  }
  # Return design matrix
  #return(X)
  return(Xopt.selected)
}
if (F) {
  #setwd("C:/Users/cbe117/School/DOE/SMED/SMED-Code")
  source('TestFunctions.R')
  source('C:/Users/cbe117/School/DOE/Codes/contour/contourfilled/R/contourfilled.R')
  SMED(banana,p=2,n=10,max.time=.2)
  SMED(function(xx){xx[1]+xx[2]^2-sin(2*pi*xx[3])},p=3,n=10,max.time=.2)
}