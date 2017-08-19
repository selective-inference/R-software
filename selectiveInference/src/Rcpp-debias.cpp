#include <Rcpp.h>      // need to include the main Rcpp header file 
#include <debias.h>    // where find_one_row_void is defined

// [[Rcpp::export]]
Rcpp::List find_one_row_debiasingM(Rcpp::NumericMatrix Sigma,
				   double bound,
				   int maxiter,
				   Rcpp::NumericVector theta,
				   Rcpp::NumericVector linear_func,
				   Rcpp::NumericVector gradient,
				   Rcpp::IntegerVector ever_active,
				   Rcpp::IntegerVector nactive
				   ) {

  int nrow = Sigma.nrow(); // number of features

  // Active set

  int irow;

  // Extract the diagonal
  Rcpp::NumericVector Sigma_diag(nrow);
  double *sigma_diag_p = Sigma_diag.begin();

  for (irow=0; irow<nrow; irow++) {
    sigma_diag_p[irow] = Sigma(irow, irow);
  }
  
  // Now call our C function

  int iter = solve_qp((double *) Sigma.begin(),
		      (double *) linear_func.begin(),
		      (double *) Sigma_diag.begin(),
		      (double *) gradient.begin(),
		      (int *) ever_active.begin(),
		      (int *) nactive.begin(),
		      nrow,
		      bound,
		      (double *) theta.begin(),
		      maxiter);
  
  // Check whether feasible

  int kkt_check = check_KKT(theta.begin(),
			    gradient.begin(),
			    nrow,
			    bound);

  return(Rcpp::List::create(Rcpp::Named("soln") = theta,
			    Rcpp::Named("gradient") = gradient,
			    Rcpp::Named("linear_func") = linear_func,
			    Rcpp::Named("iter") = iter,
			    Rcpp::Named("kkt_check") = kkt_check,
			    Rcpp::Named("ever_active") = ever_active,
			    Rcpp::Named("nactive") = nactive));

}
