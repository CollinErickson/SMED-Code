if (F){
  # Scratch work, look below for function
  f <- Vectorize(function(xx){dbeta(xx,4,2)})
  curve(f)
  qq <- function(xx){f(xx)^-.5}
  f_min <- Vectorize(function(xnew,xall,k=4) {
    #min(
      qq(xnew)^k*sum(sapply(xall,function(xx){(qq(xx)/abs(xx-xnew))^k}))
        #,1e100) # Thought Inf was a problem, still might be
  }, vectorize.args='xnew')
  curve(f_min(x,c(.75001)))
  
  # Initialize with mode
  curve(f)
  #X <- c(.75)
  xavail <- runif(1000)
  xmodeind <- which.max(f(xavail))
  X <- xavail[xmodeind]
  text(X,f(X),1,col='red')
  xavail <- xavail[-xmodeind]
  for(i in 2:30) {
    #f_opt <- optimize(function(xx){f_min(xx,X)},interval = c(0,1),temp='SANN')
    #xnew <- f_opt$min
    f_opt <- which.min(f_min(xavail,X))
    xnew <- xavail[f_opt]
    X <- c(X,xnew)
    xavail <- xavail[-f_opt]
    text(xnew,f(xnew),i,col='red')
    points(xnew,0,col=2)
  }
  #curve(f_min(x,X),ylim=c(1e6,1e11))
  hist(X,freq=F)
  curve(f,add=T,col='red')
}

# Here's the function
SMED_1D <- function(f,n=10,navail=1000,kk=4,xmin=0,xmax=1,include.histogram=F) {
  # Function that implements SMED in 1D. Estimates probability distribution
  #  f: density function
  #  n: # of pts to select
  #  navail: # of pts to use as candidates
  #  kk: power to use for MED. k=4p is default
  #  xmin, xmax: limits on inputs
  
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
  return(X)
}
if (F) {
  # Two examples from paper
  SMED_1D(function(xx)dbeta(xx,4,2),n=30)
  SMED_1D(dnorm,xmin=-3,xmax=3,n=30)
}