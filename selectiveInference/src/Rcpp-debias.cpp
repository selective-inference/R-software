#include <Rcpp.h>      // need to include the main Rcpp header file 
#include <debias.h>    // where find_one_row_void is defined

// Below, the gradient should be equal to Sigma * theta + linear_func!!
// No check is done on this.

// [[Rcpp::export]]
Rcpp::List solve_QP(Rcpp::NumericMatrix Sigma,
		    double bound,
		    int maxiter,
		    Rcpp::NumericVector theta,
		    Rcpp::NumericVector linear_func,
		    Rcpp::NumericVector gradient,
		    Rcpp::IntegerVector ever_active,
		    Rcpp::IntegerVector nactive,
		    double kkt_tol,
		    double objective_tol,
		    int max_active
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
		      maxiter,
		      kkt_tol,
		      objective_tol,
		      max_active);
  
  // Check whether feasible

  int kkt_check = check_KKT_qp(theta.begin(),
			       gradient.begin(),
			       nrow,
			       bound,
			       kkt_tol);

  int max_active_check = (*(nactive.begin()) >= max_active);

  return(Rcpp::List::create(Rcpp::Named("soln") = theta,
			    Rcpp::Named("gradient") = gradient,
			    Rcpp::Named("linear_func") = linear_func,
			    Rcpp::Named("iter") = iter,
			    Rcpp::Named("kkt_check") = kkt_check,
			    Rcpp::Named("ever_active") = ever_active,
			    Rcpp::Named("nactive") = nactive,
			    Rcpp::Named("max_active_check") = max_active_check));

}
