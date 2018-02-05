#' SMED 2D
#'
#' @param f Function
#' @param n Number of points to select
#' @param max.time max.time for GenSA optimization for each point
#'
#' @return Matrix of points
#' @export
#'
#' @examples
#' SMED_2D(TestFunctions::banana, 10, .1)
SMED_2D <- function(f,n=10,max.time=NULL) {
  # Function for SMED in 2D
  # Input:
  #  f: function
  #  n: # of pts to select
  #  nc: # of pts in contour plot
  #  max.time: max.time for GenSA optimization for each point

  # source('TestFunctions.R')
  # source('myfilledcontour.R')

  p <- 2 # dimension
  k <- 4*p # MED distance thing
  GenSA.controls <- list(trace.mat=F) # Optimization parameters
  if(!is.null(max.time)) GenSA.controls[['max.time']] <- max.time

  # Charge function qq
  qq <- function(xx){f(xx)^-(1/(2*p))}
  # Function we will optimize
  f_min <- function(xnew,xall,kk) {
    qq(xnew)^kk*sum(apply(xall,1,function(xx){(qq(xx)/(sqrt(sum((xx-xnew)^2))))^kk}))
  }

  # Get contour plot
  #my.filled.contour.func(f,nlevels=5)
  # contourfilled.func(f,nlevels=5)
  reset.plot <- ContourFunctions::cf_func(f, reset.par=F)

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
    X <- rbind(X,unname(xnew))
    text(x=xnew[1],y=xnew[2],labels=i,col=1)
  }
  reset.plot() # Reset plot
  # Return design matrix
  return(X)
}
if (F) {
  #setwd("C:/Users/cbe117/School/DOE/SMED/SMED-Code")
  source('TestFunctions.R')
  source('C:/Users/cbe117/School/DOE/Codes/contour/contourfilled/R/contourfilled.R')
  SMED_2D(banana,n=10,max.time=.2)
  SMED_2D(banana,n=40,max.time=.2)
}
