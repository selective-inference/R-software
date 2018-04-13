#include <Rcpp.h>      // need to include the main Rcpp header file 
#include <debias.h>    // where solve_QP, solve_QP_wide are defined

// Below, the gradient should be equal to Sigma * theta + linear_func!!
// No check is done on this.

// [[Rcpp::export]]
Rcpp::List solve_QP(Rcpp::NumericMatrix Sigma,
		    double bound,
		    int max_iter,
		    Rcpp::NumericVector theta,
		    Rcpp::NumericVector linear_func,
		    Rcpp::NumericVector gradient,
		    Rcpp::IntegerVector ever_active,
		    Rcpp::IntegerVector nactive,
		    double kkt_tol,
		    double objective_tol,
		    double parameter_tol,
		    int max_active,
		    int kkt_stop,
		    int objective_stop,
		    int param_stop
		    ) {

  int nrow = Sigma.nrow(); // number of features

  // Active set

  int irow;

  // Extract the diagonal
  Rcpp::NumericVector Sigma_diag(nrow);
  double *sigma_diag_p = Sigma_diag.begin();

  Rcpp::NumericVector theta_old(nrow);

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
		      (double *) theta_old.begin(),
		      max_iter,
		      kkt_tol,
		      objective_tol,
		      parameter_tol,
		      max_active,
		      kkt_stop,
		      objective_stop,
		      param_stop);
  
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
			 Rcpp::NumericVector bound,
			 double ridge_term,
			 int max_iter,
			 Rcpp::NumericVector theta,
			 Rcpp::NumericVector linear_func,
			 Rcpp::NumericVector gradient,
			 Rcpp::NumericVector X_theta,
			 Rcpp::IntegerVector ever_active,
			 Rcpp::IntegerVector nactive,
			 double kkt_tol,
			 double objective_tol,
			 double parameter_tol,
			 int max_active,
			 int kkt_stop,
			 int objective_stop,
			 int param_stop
			 ) {

  int ncase = X.nrow(); // number of cases
  int nfeature = X.ncol(); // number of features

  // Active set

  int icase, ifeature;

  // A vector to keep track of gradient updates

  Rcpp::IntegerVector need_update(nfeature);

  Rcpp::NumericVector theta_old(nfeature);

  // Extract the diagonal -- divide by ncase

  Rcpp::NumericVector nndef_diag(nfeature);
  double *nndef_diag_p = nndef_diag.begin();

  for (ifeature=0; ifeature<nfeature; ifeature++) {
    nndef_diag_p[ifeature] = 0;
    for (icase=0; icase<ncase; icase++) {
      nndef_diag_p[ifeature] += X(icase, ifeature) * X(icase, ifeature);
    }
    nndef_diag_p[ifeature] = nndef_diag_p[ifeature] / ncase;
  }
  
  // Now call our C function

  int iter = solve_wide((double *) X.begin(),
			(double *) X_theta.begin(),
			(double *) linear_func.begin(),
			(double *) nndef_diag.begin(),
			(double *) gradient.begin(),
			(int *) need_update.begin(),
			(int *) ever_active.begin(),
			(int *) nactive.begin(),
			ncase,
			nfeature,
			(double *) bound.begin(),
			ridge_term,
			(double *) theta.begin(),
			(double *) theta_old.begin(),
			max_iter,
			kkt_tol,
			objective_tol,
			parameter_tol,
			max_active,
			kkt_stop,
			objective_stop,
			param_stop);
  
  // Check whether feasible

  int kkt_check = check_KKT_wide((double *) theta.begin(),
				 (double *) gradient.begin(),
				 (double *) X_theta.begin(),
				 (double *) X.begin(),
				 (double *) linear_func.begin(),
				 (int *) need_update.begin(),
				 ncase,
				 nfeature,
				 (double *) bound.begin(),
				 ridge_term,
				 kkt_tol);

  int max_active_check = (*(nactive.begin()) >= max_active);

  // Make sure gradient is updated -- essentially a matrix multiply

  update_gradient_wide((double *) gradient.begin(),
		       (double *) X_theta.begin(),
		       (double *) X.begin(),
		       (double *) linear_func.begin(),
		       (int *) need_update.begin(),
		       ncase,
		       nfeature);

  return(Rcpp::List::create(Rcpp::Named("soln") = theta,
			    Rcpp::Named("gradient") = gradient,
			    Rcpp::Named("X_theta") = X_theta,
			    Rcpp::Named("linear_func") = linear_func,
			    Rcpp::Named("iter") = iter,
			    Rcpp::Named("kkt_check") = kkt_check,
			    Rcpp::Named("ever_active") = ever_active,
			    Rcpp::Named("nactive") = nactive,
			    Rcpp::Named("max_active_check") = max_active_check));

}
