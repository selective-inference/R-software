#include <Rcpp.h>      // need to include the main Rcpp header file 
#include <debias.h>    // where find_one_row_void is defined

// [[Rcpp::export]]
Rcpp::NumericVector find_one_row(Rcpp::NumericMatrix Sigma,
				 int row, // 0-based 
				 double bound,
				 int maxiter) {

  int nrow = Sigma.nrow(); // number of features

  // Active set

  int irow;
  Rcpp::IntegerVector nactive(1); // An array so we can easily modify it
  Rcpp::IntegerVector ever_active(1);
  int *ever_active_p = ever_active.begin();
  *ever_active_p = row;

  // Extract the diagonal
  Rcpp::NumericVector Sigma_diag(nrow);
  double *sigma_p = Sigma_diag.begin();

  for (irow=0; irow<nrow; irow++) {
    sigma_p[irow] = Sigma(irow, irow);
  }
  
  // The solution and its product with Sigma

  Rcpp::NumericVector theta(nrow);
  Rcpp::NumericVector Sigma_theta(nrow);
  
  // Now call our C function

  find_one_row_void((double *) Sigma.begin(),
		    (double *) Sigma_diag.begin(),
		    (double *) Sigma_theta.begin(),
		    (int *) ever_active.begin(),
		    (int *) nactive.begin(),
		    (int *) &nrow,
		    (double *) &bound,
		    (double *) theta.begin(),
		    (int *) &maxiter,
		    (int *) &row);

  return theta;
}
