#include <Rcpp.h>      // need to include the main Rcpp header file 
#include <debias.h>    // where find_one_row_void is defined

// [[Rcpp::export]]
Rcpp::List find_one_row_debiasingM(Rcpp::NumericMatrix Sigma,
				   int row, // 0-based 
				   double bound,
				   int maxiter,
				   Rcpp::NumericVector theta,
				   Rcpp::NumericVector Sigma_theta) {

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
  
  // Now call our C function

  int iter = find_one_row_((double *) Sigma.begin(),
			   (double *) Sigma_diag.begin(),
			   (double *) Sigma_theta.begin(),
			   (int *) ever_active.begin(),
			   (int *) nactive.begin(),
			   nrow,
			   bound,
			   (double *) theta.begin(),
			   maxiter,
			   row);
  
  // Check whether feasible

  int kkt_check = check_KKT(theta.begin(),
			    Sigma_theta.begin(),
			    nrow,
			    row,
			    bound);

  return(Rcpp::List::create(Rcpp::Named("soln") = theta,
			    Rcpp::Named("Sigma_soln") = Sigma_theta,
			    Rcpp::Named("iter") = iter,
			    Rcpp::Named("kkt_check") = kkt_check));

}
