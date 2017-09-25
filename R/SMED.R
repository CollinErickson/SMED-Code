#' SMED
#'
#' @param f Function
#' @param p Number of dimensions
#' @param n Number of points
#' @param nc Number of points for contour plot
#' @param max.time max.time for GenSA optimization for each point
#' @param X0 Matrix of initial points
#' @param Xopt Matrix of candidate points
#' @importFrom graphics curve hist pairs par plot points text
#'
#' @return Points selected
#' @export
#'
#' @examples
#' SMED(function(xx){xx[1]+xx[2]^2-sin(2*pi*xx[3])},p=3,n=10,max.time=.2)
SMED <- function(f,p,n=10,nc=100,max.time=NULL, X0=NULL, Xopt=NULL) {
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
    # contourfilled.func(f,nlevels=5)
    ContourFunctions::cf_func(f)
  }

  already_got_one <- is.null(X)
  if (is.null(X)) {
    # Initialize with mode
    gsa.out <- GenSA::GenSA(par=NULL,fn=function(xx)-f(xx),lower=rep(0,p),upper=rep(1,p),control = list(trace.mat=F))
    X <- matrix(gsa.out$par,1,p)
    if (p==2) {
      text(X[1],X[2],labels=1,col=1,pch=1)
    }
  } else {
    points(X, pch=19)
  }

  Y <- apply(X,1,f)
  qqX <- apply(X,1,qq)
  Delta <- .01 * max(Y)
  keep.Delta <- (Y > Delta)

  # Get rest of points
  for(i in (1 + already_got_one):n) {#print(keep.Delta)
    # Use log scale for optimization, had trouble before when numbers were 1e88
    opt.func <- function(xx)log(f_min(xx,X[keep.Delta, , drop=F],qqX[keep.Delta],kk=k))
    if (is.null(Xopt)) {
      gsa.out <- GenSA::GenSA(par=NULL,
                              fn=opt.func,
                              lower=rep(0,p),upper=rep(1,p),
                              control = GenSA.controls)
      # Add new point
      xnew <- gsa.out$par
    } else {
      func.vals <- apply(Xopt, 1, opt.func)
      best <- which.min(func.vals)
      xnew <- Xopt[best,]
      Xopt <- Xopt[-best, , drop=F]
    }
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
  # source('TestFunctions.R')
  # source('C:/Users/cbe117/School/DOE/Codes/contour/contourfilled/R/contourfilled.R')
  SMED(TestFunctions::banana,p=2,n=10,max.time=.2)
  SMED(function(xx){xx[1]+xx[2]^2-sin(2*pi*xx[3])},p=3,n=10,max.time=.2)
}
