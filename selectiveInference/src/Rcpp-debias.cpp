#include <Rcpp.h>      // need to include the main Rcpp header file 
#include <debias.h>    // where find_one_row is defined

// [[Rcpp::export]]
Rcpp::NumericVector find_one_row(Rcpp::NumericMatrix Sigma,
				 int row, // 0-based 
				 double bound,
				 int maxiter) {

  int nrow = Sigma.nrow(); // number of features
  int nactive = 1;
  Rcpp::IntegerVector ever_active(1);
  Rcpp::NumericVector Sigma_diag(nrow);
  Rcpp::NumericVector Sigma_theta(nrow);

}
 Rcpp::NumericVector xUL,
                       int maxEval, double absErr, double tol, int vectorInterface, unsigned norm) {

  count = 0; /* Zero count */
  fun = f;

  Rcpp::NumericVector integral(fDim);
  Rcpp::NumericVector errVals(fDim);
  int retCode;

  // Rcpp::Rcout<<"Call Integrator" <<std::endl;
  if (vectorInterface) {
    retCode = hcubature_v(fDim, fWrapper_v, NULL,
			  xLL.size(), xLL.begin(), xUL.begin(),
			  maxEval, absErr, tol, (error_norm) norm,
			  integral.begin(), errVals.begin());
  } else {
    retCode = hcubature(fDim, fWrapper, NULL,
			xLL.size(), xLL.begin(), xUL.begin(),
			maxEval, absErr, tol, (error_norm) norm,
			integral.begin(), errVals.begin());
  }
  return Rcpp::List::create(
			    Rcpp::_["integral"] = integral,
			    Rcpp::_["error"] = errVals,
			    Rcpp::_["functionEvaluations"] = count,
			    Rcpp::_["returnCode"] = retCode);
}

