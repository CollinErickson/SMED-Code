Rcpp::cppFunction('
double log_f_minC(NumericVector xnew, double qqxnew, NumericMatrix xall, NumericVector qqall, double kk) {
  int n = xall.nrow();
  int m = xall.ncol();
  double total = 0;
  double dist = 0;
  for(int i = 0; i < n; ++i) {
    dist = 0;
    for(int j = 0; j < m; ++j) {
      dist += pow(xall(i,j) - xnew[j], 2.0);
    }
    total += pow(qqall[i] / sqrt(dist), kk);
  }

  total *= pow(qqxnew, kk);
  total = log(total);
  return total;
}')


Rcpp::cppFunction('
                  IntegerVector SMED_selectC(Function f, int n, NumericMatrix X0, NumericMatrix Xopt) {  
                  int p = X0.ncol();
                  int k = 4 * p;
                  
                  // initiate values for X0
                  NumericVector Y(X0.nrow());
                  for (int i=0; i < Y.size(); ++i) {
                  Y[i] = as<double>(f(X0(i, _)));
                  }
                  NumericVector qqX(X0.nrow());
                  for (int i=0; i < qqX.size(); ++i) {
                  qqX[i] = pow(Y[i], -1.0 / (2 * p));
                  }
                  
                  // initiate values for Xopt
                  NumericVector Yopt(Xopt.nrow());
                  for (int i=0; i < Yopt.size(); ++i) {
                  Yopt[i] = as<double>(f(Xopt(i, _)));
                  }
                  NumericVector qqXopt(Xopt.nrow());
                  for (int i=0; i < qqXopt.size(); ++i) {
                  qqXopt[i] = pow(Yopt[i], -1.0 / (2 * p));
                  }
                  double Delta = .01 * max(Y);
                  LogicalVector keepDelta = (Y > Delta);
                  IntegerVector XoptSelectedIndsOrder(n);
                  LogicalVector XoptSelected(Xopt.nrow(), false);
                  
                  
                  //double total = 0;
                  double dist = 0;
                  double funcValMin;
                  double funcValMinInd = -1;
                  double funcVal;
                  NumericVector funcVals(Xopt.nrow());
                  
                  // Pick n next with SMED
                  for(int i = 0; i < n; ++i) {
                  //NumericVector funcVals(Xopt.nrow());
                  for(int ii=0; ii < funcVals.size(); ++ii){
                  funcVals(ii) = 0;
                  }
                  // Loop over points
                  for(int j = 0; j < Xopt.nrow(); ++j) {
                  if (!XoptSelected[j]) {
                  funcVals[j] = 0;
                  funcVal = 0;
                  // Loop over X0 (keptDelta) and selected Xopt to get funcVal
                  for(int l = 0; l < X0.nrow(); ++l) {
                  dist = sum(pow(Xopt(j, _) - X0(l, _), 2.0));
                  funcVals[j] += pow(qqX[l] / sqrt(dist), k);
                  funcVal += pow(qqX[l] / sqrt(dist), k);
                  }
                  // Loop over points already selected
                  for(int l = 0; l < Xopt.nrow(); ++l) {
                  if (XoptSelected[l]) {
                  dist = sum(pow(Xopt(j, _) - Xopt(l, _), 2.0));
                  funcVals[j] += pow(qqXopt[l] / sqrt(dist), k);
                  funcVal += pow(qqXopt[l] / sqrt(dist), k);
                  }
                  }
                  
                  funcVal *= pow(qqXopt[j], k);
                  funcVals[j] *= pow(qqXopt[j], k);
                  
                  // Check if it is the best
                  if ((funcValMinInd < 0) | (funcVal < funcValMin)) {
                  funcValMin = funcVal;
                  funcValMinInd = j;
                  } 
                  }
                  
                  } // end loop over Xopt points
                  XoptSelectedIndsOrder[i] = funcValMinInd;
                  XoptSelected[funcValMinInd] = true;
                  funcValMinInd = -1;
                  } // end loop to select n
                  
                  /*# Get rest of points
                  
                  Delta <- .01 * max(Y) #.01 * (n0 / (n + i)) ^ (1 / p) * max(Y)
                  #keep.Delta <- (Y > Delta)
                  keep.Delta <- ifelse(1:nrow(X) <= n0, Y > Delta, T)
                  
                  }
                  */
                  
                  //return funcVals;
                  //return funcValMin;
                  //return qqXopt;//as<double>(f(X0(0,_)));
                  return XoptSelectedIndsOrder + 1;
                  }')
print(
  #SMEDC(function(x)sum(abs(sin(x))), 2, matrix(runif(16),8,2), matrix(runif(8),4,2))
  #SMEDC(function(x)(sin(x*2*pi)^2), 3, matrix(c(.2,.3,.4,.25,.14,.8,.75,.93),8,1), matrix(c(.91,.21,.9,.77,.85),ncol=1))
  SMEDC(function(x)(sin(x*2*pi)^2), 1, matrix(c(.2,.3,.7),ncol=1), matrix(c(.301,.91,.21,.9,.77,.85,.8,.99),ncol=1))
  )
SMED_select(function(x)(sin(x*2*pi)^2), 7, matrix(c(.2,.3,.7),ncol=1), matrix(c(.301,.91,.21,.9,.77,.85,.8,.99),ncol=1))

#SMEDC(function(x)sum(abs(sin(x))), 2, matrix(runif(16),8,2), matrix(runif(8),4,2))


#f_min <- function(xnew,xall, qqall,kk) {
#  qq(xnew)^kk*sum(sapply(1:nrow(xall),function(ii){(qqall[ii]/(sqrt(sum((xall[ii,]-xnew)^2))))^kk}))
#}
#qq <- function(xx){f(xx)^-(1/(2*p))}

SMED_selectSLOW <- function(f, n, X0=NULL, Xopt=NULL) {#browser()
  # Function for SMED in 2D
  # Input:
  #  f: function
  #  p: # of dimensions
  #  n: # of pts to select
  #  X0: current design points
  #  Xopt: points to consider
  
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
  f_min <- function(xnew,qqnew,xall, qqall,kk) {
    qq(xnew)^kk*sum(apply(xall,1,function(xx){(qq(xx)/(sqrt(sum((xx-xnew)^2))))^kk}))
    #qq(xnew)^kk*sum(sapply(1:nrow(xall),function(ii){(qqall[ii]/(sqrt(sum((xall[ii,]-xnew)^2))))^kk}))
  }
  

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
    opt.func <- function(xx){#if(nrow(X[keep.Delta,])>50)browser()
      #log_f_minC(xx,qq(xx), X[keep.Delta, , drop=F],qqX[keep.Delta],kk=k)
      log(f_min(xx, qq(xx), X[keep.Delta, , drop=F],qqX[keep.Delta],kk=k))
    }
  
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
    
  }
  # Return indices of selected points
  return(Xopt.selected)
}
if (F) {
  #setwd("C:/Users/cbe117/School/DOE/SMED/SMED-Code")
  source('TestFunctions.R')
  source('C:/Users/cbe117/School/DOE/Codes/contour/contourfilled/R/contourfilled.R')
  SMED(banana,p=2,n=10,max.time=.2)
  SMED(function(xx){xx[1]+xx[2]^2-sin(2*pi*xx[3])},p=3,n=10,max.time=.2)
}