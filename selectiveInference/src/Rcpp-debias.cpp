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


// [[Rcpp::export]]
Rcpp::List solve_QP_wide(Rcpp::NumericMatrix X,
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

  int nrow = X.nrow(); // number of cases
  int ncol = X.ncol(); // number of features

  // Active set

  int irow, icol;

  // Extract the diagonal
  Rcpp::NumericVector X_diag(ncol);
  double *X_diag_p = X_diag.begin();

  for (icol=0; icol<ncol; icol++) {
    X_diag_p[irow] = 0;
    for (irow=0; irow<nrow; irow++) {
      X_diag_p[irow] += X(irow, icol) * X(irow, icol);
    }
    X_diag_p[irow] = X_diag_p[irow] / nrow;
  }
  
  // Now call our C function

  int iter = solve_wide((double *) X.begin(),
			(double *) linear_func.begin(),
			(double *) X_diag.begin(),
			(double *) gradient.begin(),
			(int *) ever_active.begin(),
			(int *) nactive.begin(),
			nrow,
			ncol,
			bound,
			(double *) theta.begin(),
			maxiter,
			kkt_tol,
			objective_tol,
			max_active);
  
  // Check whether feasible

  int kkt_check = check_KKT_wide(theta.begin(), // This is the same function as check_KKT_qp at the moment!!
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
