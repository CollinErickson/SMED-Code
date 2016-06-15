SMED_1D <- function(f,n=10,navail=1000,kk=4,xmin=0,xmax=1,include.histogram=F) {
  # Function that implements SMED in 1D. Estimates probability distribution
  #  f: density function
  #  n: # of pts to select
  #  navail: # of pts to use as candidates
  #  kk: power to use for MED. k=4p is default
  #  xmin, xmax: limits on inputs
  #  include.histogram - whether to plot histogram at end
  
  # Charge function
  qq <- function(xx){f(xx)^-.5} 
  
  # Function we minimize to select next point
  f_min <- Vectorize(function(xnew,xall,k=kk) {
    qq(xnew)^k*sum(sapply(xall,function(xx){(qq(xx)/abs(xx-xnew))^k}))
  }, vectorize.args='xnew')
  
  # Initialize with mode
  xavail <- runif(navail,xmin,xmax)
  xmodeind <- which.max(f(xavail))
  X <- xavail[xmodeind]
  xavail <- xavail[-xmodeind] # Removes selected pts from available pts  
  
  # Plot density and points
  curve(f,from=xmin,to=xmax)
  text(X,f(X),1,col=4)
  points(X,0,col=2)
  
  # Sequentially pick rest
  for(i in 2:n) {
    # Find f_opt minimum
    f_opt <- which.min(f_min(xavail,X))
    xnew <- xavail[f_opt]
    # Update X and plot new point
    X <- c(X,xnew)
    xavail <- xavail[-f_opt]
    text(xnew,f(xnew),i,col=4)
    points(xnew,0,col=2)
  }
  if (include.histogram) hist(X,breaks=20)
  return(X)
}
if (F) {
  # Two examples from paper
  SMED_1D(function(xx)dbeta(xx,4,2),n=30)
  SMED_1D(dnorm,xmin=-3,xmax=3,n=30)
}