#include <Rcpp.h>                // need to include the main Rcpp header file 
#include <randomized_lasso.h>    // where densities are defined
#include <selective_mle.h>       // where barrier_solve is defined

// [[Rcpp::export]]
Rcpp::NumericVector log_density_gaussian_(double noise_scale,                         // Scale of randomization
					  Rcpp::NumericMatrix internal_linear,        // A_D -- linear part for data
					  Rcpp::NumericMatrix internal_state,         // D -- data state -- matrix of shape (nopt, npts)
					  Rcpp::NumericMatrix optimization_linear,    // A_O -- linear part for optimization variables
					  Rcpp::NumericMatrix optimization_state,     // O -- optimization state -- matrix of shape (ninternal, npts)
					  Rcpp::NumericVector offset) {               // h -- offset in affine transform -- "p" dimensional 

  int npt = internal_state.ncol();         // Function is vectorized
  if (optimization_state.ncol() != npt) {  // Assuming each column is an internal or opt state because arrays are column major
    Rcpp::stop("Number of optimization samples should equal the number of (internally represented) data.");
  }

  int ndim = optimization_linear.nrow();  
  if (internal_linear.nrow() != ndim) {  
    Rcpp::stop("Dimension of optimization range should be the same as the dimension of the data range.");
  } 
  int ninternal = internal_linear.ncol();
  int noptimization = optimization_linear.ncol();

  Rcpp::NumericVector result(npt);

  int ipt;
  for (ipt=0; ipt<npt; ipt++) {
    result[ipt] = log_density_gaussian(noise_scale,
				       ndim,
				       ninternal,
				       noptimization,
				       (double *) internal_linear.begin(),
				       ((double *) internal_state.begin() + ipt * ninternal),
				       (double *) optimization_linear.begin(),
				       ((double *) optimization_state.begin() + ipt * noptimization),
				       (double *) offset.begin());
  }

  return(result);
}

// [[Rcpp::export]]
Rcpp::NumericVector log_density_gaussian_conditional_(double noise_scale,                         // Scale of randomization
						      Rcpp::NumericMatrix optimization_linear,    // A_O -- linear part for optimization variables
						      Rcpp::NumericMatrix optimization_state,     // O -- optimization state -- matrix of shape (ninternal, npts)
						      Rcpp::NumericVector offset) {               // h -- offset in affine transform -- "p" dimensional 

  int npt = optimization_state.ncol();         // Function is vectorized
  int ndim = optimization_linear.nrow();  
  int noptimization = optimization_linear.ncol();

  Rcpp::NumericVector result(npt);

  int ipt;
  for (ipt=0; ipt<npt; ipt++) {
    result[ipt] = log_density_gaussian_conditional(noise_scale,
						   ndim,
						   noptimization,
						   (double *) optimization_linear.begin(),
						   ((double *) optimization_state.begin() + ipt * noptimization),
						   (double *) offset.begin());
  }

  return(result);
}

// [[Rcpp::export]]
Rcpp::NumericVector log_density_laplace_(double noise_scale,                         // Scale of randomization
					 Rcpp::NumericMatrix internal_linear,        // A_D -- linear part for data
					 Rcpp::NumericMatrix internal_state,         // D -- data state -- matrix of shape (nopt, npts)
					 Rcpp::NumericMatrix optimization_linear,    // A_O -- linear part for optimization variables
					 Rcpp::NumericMatrix optimization_state,     // O -- optimization state -- matrix of shape (ninternal, npts)
					 Rcpp::NumericVector offset) {               // h -- offset in affine transform -- "p" dimensional 

  int npt = internal_state.ncol();         // Function is vectorized
  if (optimization_state.ncol() != npt) {  // Assuming each column is an internal or opt state because arrays are column major
    Rcpp::stop("Number of optimization samples should equal the number of (internally represented) data.");
  }

  int ndim = optimization_linear.nrow();  
  if (internal_linear.nrow() != ndim) {  
    Rcpp::stop("Dimension of optimization range should be the same as the dimension of the data range.");
  } 
  int ninternal = internal_linear.ncol();
  int noptimization = optimization_linear.ncol();

  Rcpp::NumericVector result(npt);

  int ipt;
  for (ipt=0; ipt<npt; ipt++) {
    result[ipt] = log_density_laplace(noise_scale,
				      ndim,
				      ninternal,
				      noptimization,
				      (double *) internal_linear.begin(),
				      ((double *) internal_state.begin() + ipt * ninternal),
				      (double *) optimization_linear.begin(),
				      ((double *) optimization_state.begin() + ipt * noptimization),
				      (double *) offset.begin());
  }

  return(result);
}

// [[Rcpp::export]]
Rcpp::NumericVector log_density_laplace_conditional_(double noise_scale,                         // Scale of randomization
						     Rcpp::NumericMatrix optimization_linear,    // A_O -- linear part for optimization variables
						     Rcpp::NumericMatrix optimization_state,     // O -- optimization state -- matrix of shape (ninternal, npts)
						     Rcpp::NumericVector offset) {               // h -- offset in affine transform -- "p" dimensional 

  int npt = optimization_state.ncol();         // Function is vectorized
  int ndim = optimization_linear.nrow();  
  int noptimization = optimization_linear.ncol();

  Rcpp::NumericVector result(npt);

  int ipt;
  for (ipt=0; ipt<npt; ipt++) {
    result[ipt] = log_density_laplace_conditional(noise_scale,
						  ndim,
						  noptimization,
						  (double *) optimization_linear.begin(),
						  ((double *) optimization_state.begin() + ipt * noptimization),
						  (double *) offset.begin());
  }

  return(result);
}

// [[Rcpp::export]]
Rcpp::List solve_barrier_(Rcpp::NumericVector conjugate_arg,     // Argument to conjugate
			  Rcpp::NumericMatrix precision,         // Precision matrix in conjugate optimization problem
			  Rcpp::NumericVector feasible_point,    // Feasible point -- must be nonnegative
			  int max_iter,                          // How many iterations to run
			  int min_iter,                          // Minimum iterations to run
			  double value_tol,                      // Tolerance for convergence
			  double initial_step) {                 // Initial step size

  int ndim = precision.ncol();
  int idim;

  Rcpp::NumericVector gradient(ndim);
  Rcpp::NumericVector opt_variable(ndim);
  Rcpp::NumericVector opt_proposed(ndim);
  Rcpp::NumericVector scaling(ndim);

  double *scaling_ptr = scaling.begin();
  double *opt_ptr = opt_variable.begin();

  for (idim=0; idim<ndim; idim++) {
    *scaling_ptr = precision(idim, idim); scaling_ptr++;
    *opt_ptr = feasible_point(idim); opt_ptr++;
  }

  double value = barrier_solve((double *) gradient.begin(),                   
			       (double *) opt_variable.begin(),               
			       (double *) opt_proposed.begin(),               
			       (double *) conjugate_arg.begin(),              
			       (double *) precision.begin(),                  
			       (double *) scaling.begin(),                    
			       ndim,                                          
			       max_iter,                                      
			       min_iter,                                      
			       value_tol,                                     
			       initial_step);                                 

  return(Rcpp::List::create(Rcpp::Named("soln") = opt_variable,
			    Rcpp::Named("value") = value,
			    Rcpp::Named("gradient") = gradient));
}
