#' Select SMED points from a set of options
#'
#' SMED_select: R implementation
#'
#' @useDynLib SMED
#' @importFrom Rcpp evalCpp
#' @param f Function
#' @param n Number of points to select
#' @param X0 Design matrix of already selected points (points in rows)
#' @param Xopt Matrix where each row is a candidate point
#' @param theta Correlation parameters in each dimension for distance scaling
#' @param useC Should SMED_select call SMED_selectC?
#'
#' @return Vector with indices of selected points
#' @export
#'
#' @examples
#' X0 <- matrix(c(.45,.1,.55,.9),ncol=2,byrow=TRUE)
#' Xopt <- matrix(c(.45,.5,.5,.2),ncol=2,byrow=TRUE)
#' plot(NULL, xlim=c(0,1), ylim=c(0,1))
#' points(X0)
#' text(Xopt, col=2)
#' SMED_select(function(a) {1}, 1, X0, Xopt)
#'
#' # Checking if using theta works
#' X0 <- matrix(c(.1,.1,.5,.5,.9,.9), ncol=2, byrow=TRUE)
#' Xopt <- matrix(c(.1,.9,.3,.5,.5,.3), ncol=2, byrow=TRUE)
#' plot(X0,xlim=0:1,ylim=0:1, pch=19)
#' text(Xopt, col=2, pch=19)
#' f <- function(x){1} # all have same charge, so only theta matters
#' # Each dimension equal, pick 1 since it is most spread out
#' SMED_select(f=f, n=3, X0=X0, Xopt=Xopt)
#' # y dimension has more variability, pick 3 since it is furthest from X0 in y direction
#' SMED_select(f=f, n=3, X0=X0, Xopt=Xopt, theta=c(1,1e3))
#' # x dimension has more variability, pick 2 since it is furthest from X0 in x direction
#' SMED_select(f=f, n=3, X0=X0, Xopt=Xopt, theta=c(1,1e-3))
SMED_select <- function(f,n=1, X0=NULL, Xopt=NULL, theta=rep(1,ncol(X0)), useC=F) {
  if (useC) {
    return(SMED_selectC(f=f, n=n, X0=X0, Xopt=Xopt, theta=theta))
  }

  #p <- d # dimension
  p <- ncol(X0)
  k <- 4*p # MED distance thing
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
    #qq(xnew)^kk*sum(sapply(1:nrow(xall),function(ii){(qqall[ii]/(sqrt(sum((xall[ii,]-xnew)^2))))^kk}))
    qq(xnew)^kk*sum(sapply(1:nrow(xall),function(ii){(qqall[ii]/(sqrt(sum((xall[ii,]-xnew)^2 * theta))))^kk}))
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
  #Xopt.consider <- rep(T, nrow(X.opt)) # to only check good points when taking more than 1

  # Get rest of points
  for(i in 1:n) {#print(keep.Delta)
    # Use log scale for optimization, had trouble before when numbers were 1e88
    opt.func <- function(xx)log(f_min(xx,X[keep.Delta, , drop=F],
                                      qqX[keep.Delta],kk=k))

    func.vals <- apply(Xopt, 1, opt.func)
    best <- which.min(func.vals)
    xnew <- Xopt[best,]
    Xopt <- Xopt[-best, , drop=F]
    Xopt.selected <- c(Xopt.selected, Xopt.inds[best])
    Xopt.inds <- Xopt.inds[-best]
    #Xopt.consider <- Xopt.consider[-best]

    X <- rbind(X,unname(xnew))
    #Y <- apply(X,1,f) # could just append
    Y <- c(Y, f(xnew))
    #qqX <- apply(X,1,qq)
    qqX <- c(qqX, qq(xnew))
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
  # source('TestFunctions.R')
  # source('C:/Users/cbe117/School/DOE/Codes/contour/contourfilled/R/contourfilled.R')
  SMED(TestFunctions::banana,p=2,n=10,max.time=.2)
  SMED(function(xx){xx[1]+xx[2]^2-sin(2*pi*xx[3])},p=3,n=10,max.time=.2)
}
