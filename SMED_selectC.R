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
  //LogicalVector keepDelta = (Y > Delta);
  IntegerVector XoptSelectedIndsOrder(n);
  LogicalVector XoptSelected(Xopt.nrow(), false);
  
  
  //double total = 0;
  double dist = 0;
  double funcValMin;
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
          if (Y[l] >= Delta) {
            dist = sum(pow(Xopt(j, _) - X0(l, _), 2.0));
            //funcVals[j] += pow(qqX[l] / sqrt(dist), k);
            funcVal += pow(qqX[l] / sqrt(dist), k);
          }
        }
        // Loop over points already selected
        for(int l = 0; l < Xopt.nrow(); ++l) {
          if (XoptSelected[l]) {
            dist = sum(pow(Xopt(j, _) - Xopt(l, _), 2.0));
            //funcVals[j] += pow(qqXopt[l] / sqrt(dist), k);
            funcVal += pow(qqXopt[l] / sqrt(dist), k);
          }
        }
      
        funcVal *= pow(qqXopt[j], k);
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
}')
if (F) {
  #SMEDC(function(x)sum(abs(sin(x))), 2, matrix(runif(16),8,2), matrix(runif(8),4,2))
  #SMEDC(function(x)(sin(x*2*pi)^2), 3, matrix(c(.2,.3,.4,.25,.14,.8,.75,.93),8,1), matrix(c(.91,.21,.9,.77,.85),ncol=1))
  SMEDC(function(x)(sin(x*2*pi)^2), 1, matrix(c(.2,.3,.7),ncol=1), matrix(c(.301,.91,.21,.9,.77,.85,.8,.99),ncol=1))
  
  SMED_select(function(x)(sin(x*2*pi)^2), 7, matrix(c(.2,.3,.7),ncol=1), matrix(c(.301,.91,.21,.9,.77,.85,.8,.99),ncol=1))

}