#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

//' Select points using SMED
//'
//' SMED_selectYC takes in Y values instead of a function. It
//' uses C code and should run faster than SMED_select
//'
//' @param Y0 Y values for X0
//' @param Yopt Y values for Xopt
//' @examples
//' n = 1
//' X0 <- matrix(c(.45,.1,.55,.9),ncol=2,byrow=TRUE)
//' Xopt <- matrix(c(.45,.5,.5,.2),ncol=2,byrow=TRUE)
//' Y0 <- c(.3,.5)
//' Yopt <- c(.2,.6)
//' SMED_selectYC(n, X0, Xopt, Y0, Yopt)
//'
//' @export
//' @rdname SMED_select
// [[Rcpp::export]]
IntegerVector SMED_selectYC(int n, NumericMatrix X0, NumericMatrix Xopt, NumericVector Y0, NumericVector Yopt, NumericVector theta = NumericVector::create()) {
  int p = X0.ncol();
  int k = 4 * p;

  // If theta wasn't given, set theta to vector of 1's, ie rep(1, Xopt.nrow()) but not sure how to do it simply
  if (theta.length() == 0) {
    theta = NumericVector(Xopt.ncol());
    for (int i = 0; i < Xopt.ncol(); i++) {
      theta(i) = 1;
    }
  }
  //Rcout << theta << "\n";

  // initiate values for X0
  //NumericVector Y(X0.nrow());
  //for (int i=0; i < Y.size(); ++i) {
  //  Y[i] = as<double>(f(X0(i, _)));
  //}
  NumericVector qqX(X0.nrow());
  for (int i=0; i < qqX.size(); ++i) {
    qqX[i] = pow(Y0[i], -1.0 / (2 * p));
  }
  //Rcout << qqX << "\n";

  // initiate values for Xopt
  //NumericVector Yopt(Xopt.nrow());
  //for (int i=0; i < Yopt.size(); ++i) {
  //  Yopt[i] = as<double>(f(Xopt(i, _)));
  //}
  NumericVector qqXopt(Xopt.nrow());
  for (int i=0; i < qqXopt.size(); ++i) {
    qqXopt[i] = pow(Yopt[i], -1.0 / (2 * p));
  }
  //Rcout << qqXopt << "\n";

  double Delta = .01 * max(Y0);
  //LogicalVector keepDelta = (Y > Delta);
  IntegerVector XoptSelectedIndsOrder(n);
  LogicalVector XoptSelected(Xopt.nrow(), false);


  //double total = 0;
  double dist2 = 0; // distance squared
  double funcValMin = 0; // Shouldnt need to be initialized
  double funcValMinInd = -1;
  double funcVal;
  //NumericVector funcVals(Xopt.nrow());

  // Pick n next with SMED
  for(int i = 0; i < n; ++i) {
    //NumericVector funcVals(Xopt.nrow());
    //for(int ii=0; ii < funcVals.size(); ++ii){
    //  funcVals(ii) = 0;
    //}
    // Loop over points still available
    for(int j = 0; j < Xopt.nrow(); ++j) {
      if (!XoptSelected[j]) {
        //funcVals[j] = 0;
        funcVal = 0;
        // Loop over X0 (keptDelta) and selected Xopt to get funcVal
        for(int l = 0; l < X0.nrow(); ++l) {
          if (Y0[l] >= Delta) {
            dist2 = sum(pow(Xopt(j, _) - X0(l, _), 2) * theta); // Adding theta to scale distance by dimension
            //funcVals[j] += pow(qqX[l] / sqrt(dist2), k);
            funcVal += pow(qqX[l] / sqrt(dist2), k);
          }
        }
        // Loop over points already selected
        for(int l = 0; l < Xopt.nrow(); ++l) {
          if (XoptSelected[l]) {
            dist2 = sum(pow(Xopt(j, _) - Xopt(l, _), 2) * theta); // Add theta here too
            //funcVals[j] += pow(qqXopt[l] / sqrt(dist2), k);
            funcVal += pow(qqXopt[l] / sqrt(dist2), k);
          }
        }

        funcVal *= pow(qqXopt[j], k);
        //Rcout << j << " " << funcVal << "\n";
        //funcVals[j] *= pow(qqXopt[j], k);

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

  return XoptSelectedIndsOrder + 1;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
#timesTwo(42)
*/
