#include <math.h> // for fabs
#include <stdio.h>

// Augmented density for randomized LASSO after
// Gaussian randomization

// Described in https://arxiv.org/abs/1609.05609

// Gaussian is product of IID N(0, noise_scale^2) density
// Evaluated at A_D D + A_O O + h

// Laplace is product of IID Laplace with scale noise_scale
// Also evaluated at A_D D + A_O O + h

// Matrices are assumed in column major order! 

double log_density_gaussian(double noise_scale,             // Scale of randomization
			    int ndim,                       // Number of features -- "p"
			    int ninternal,                  // Dimension of internal data representation often 1
			    int noptimization,              // Dimension of optimization variables -- "p"
			    double *internal_linear,        // A_D -- linear part for data
			    double *internal_state,         // D -- data state
			    double *optimization_linear,    // A_O -- linear part for optimization variables
			    double *optimization_state,     // O -- optimization state
			    double *offset)                 // h -- offset in affine transform -- "p" dimensional 
{
  int irow, icol;
  double denom = 2 * noise_scale * noise_scale;
  double value = 0;
  double reconstruction = 0;
  double *offset_ptr;
  double *internal_linear_ptr;
  double *internal_state_ptr;
  double *optimization_linear_ptr;
  double *optimization_state_ptr;

  for (irow=0; irow<ndim; irow++) {

    // Compute the irow-th entry of the ndim vector

    offset_ptr = ((double *) offset + irow);
    reconstruction = *offset_ptr;

    // Internal (i.e. data) contribution
    for (icol=0; icol<ninternal; icol++) {
      
      internal_linear_ptr = ((double *) internal_linear + icol * ndim + irow);
      internal_state_ptr = ((double *) internal_state + icol);

      reconstruction += (*internal_linear_ptr) * (*internal_state_ptr);
    }

    // Optimization contribution
    for (icol=0; icol<noptimization; icol++) {
      
      optimization_linear_ptr = ((double *) optimization_linear + icol * ndim + irow);
      optimization_state_ptr = ((double *) optimization_state + icol);

      reconstruction += (*optimization_linear_ptr) * (*optimization_state_ptr);
    }

    value -= (reconstruction * reconstruction) / denom;
  }

  return(value);
}

double log_density_laplace(double noise_scale,             // Scale of randomization
			   int ndim,                       // Number of features -- "p"
			   int ninternal,                  // Dimension of internal data representation often 1
			   int noptimization,              // Dimension of optimization variables -- "p"
			   double *internal_linear,        // A_D -- linear part for data
			   double *internal_state,         // D -- data state
			   double *optimization_linear,    // A_O -- linear part for optimization variables
			   double *optimization_state,     // O -- optimization state
			   double *offset)                 // h -- offset in affine transform -- "p" dimensional 
{
  int irow, icol;
  double value = 0;
  double reconstruction = 0;
  double *offset_ptr;
  double *internal_linear_ptr;
  double *internal_state_ptr;
  double *optimization_linear_ptr;
  double *optimization_state_ptr;

  for (irow=0; irow<ndim; irow++) {

    // Compute the irow-th entry of the ndim vector

    offset_ptr = ((double *) offset + irow);
    reconstruction = *offset_ptr;

    // Internal (i.e. data) contribution
    for (icol=0; icol<ninternal; icol++) {
      
      internal_linear_ptr = ((double *) internal_linear + icol * ndim + irow);
      internal_state_ptr = ((double *) internal_state + icol);

      reconstruction += (*internal_linear_ptr) * (*internal_state_ptr);
    }

    // Optimization contribution
    for (icol=0; icol<noptimization; icol++) {
      
      optimization_linear_ptr = ((double *) optimization_linear + icol * ndim + irow);
      optimization_state_ptr = ((double *) optimization_state + icol);

      reconstruction += (*optimization_linear_ptr) * (*optimization_state_ptr);
    }

    value -= fabs(reconstruction) / noise_scale;
  }

  return(value);
}

// Keeping internal (data) state fixed

double log_density_gaussian_conditional(double noise_scale,             // Scale of randomization
					int ndim,                       // Number of features -- "p"
					int noptimization,              // Dimension of optimization variables -- "p"
					double *optimization_linear,    // A_O -- linear part for optimization variables
					double *optimization_state,     // O -- optimization state
					double *offset)                 // h -- offset in affine transform -- "p" dimensional 
{
  int irow, icol;
  double value = 0;
  double denom = 2 * noise_scale * noise_scale;
  double reconstruction = 0;
  double *offset_ptr;
  double *optimization_linear_ptr;
  double *optimization_state_ptr;

  for (irow=0; irow<ndim; irow++) {

    // Compute the irow-th entry of the ndim vector

    offset_ptr = ((double *) offset + irow);
    reconstruction = *offset_ptr;

    // Optimization contribution
    for (icol=0; icol<noptimization; icol++) {
      
      optimization_linear_ptr = ((double *) optimization_linear + icol * ndim + irow);
      optimization_state_ptr = ((double *) optimization_state + icol);

      reconstruction += (*optimization_linear_ptr) * (*optimization_state_ptr);
    }

    value -= reconstruction * reconstruction / denom;
  }

  return(value);
}

double log_density_laplace_conditional(double noise_scale,             // Scale of randomization
				       int ndim,                       // Number of features -- "p"
				       int noptimization,              // Dimension of optimization variables -- "p"
				       double *optimization_linear,    // A_O -- linear part for optimization variables
				       double *optimization_state,     // O -- optimization state
				       double *offset)                 // h -- offset in affine transform -- "p" dimensional 
{
  int irow, icol;
  double value = 0;
  double reconstruction = 0;
  double *offset_ptr;
  double *optimization_linear_ptr;
  double *optimization_state_ptr;

  for (irow=0; irow<ndim; irow++) {

    // Compute the irow-th entry of the ndim vector

    offset_ptr = ((double *) offset + irow);
    reconstruction = *offset_ptr;

    // Optimization contribution
    for (icol=0; icol<noptimization; icol++) {
      
      optimization_linear_ptr = ((double *) optimization_linear + icol * ndim + irow);
      optimization_state_ptr = ((double *) optimization_state + icol);

      reconstruction += (*optimization_linear_ptr) * (*optimization_state_ptr);
    }

    value -= fabs(reconstruction) / noise_scale;
  }

  return(value);
}

//    objective = lambda u: -u.T.dot(conjugate_arg) + u.T.dot(precision).dot(u)/2. + np.log(1.+ 1./(u / scaling)).sum()

double barrier_objective(double *opt_variable,               // Optimization variable
			 double *conjugate_arg,              // Argument to conjugate of Gaussian
			 double *precision,                  // Precision matrix of Gaussian
			 double *scaling,                    // Diagonal scaling matrix for log barrier
			 int ndim)                           // Dimension of conjugate_arg, precision
{
  int idim, jdim;
  double *conjugate_arg_ptr;
  double *opt_variable_ptr;
  double *precision_ptr;
  double *scaling_ptr;
  double product_entry, value;

  value = 0.;
  for (idim=0; idim<ndim; idim++) {
    
    // Compute an entry of precision.dot(conjugate_arg)

    product_entry = 0;
    for (jdim=0; jdim<ndim; jdim++) {
    
      precision_ptr = ((double *) precision + idim * ndim + jdim); // precision is a symmetric matrix
      opt_variable_ptr = ((double *) opt_variable + jdim);
      product_entry += (*precision_ptr) * (*opt_variable_ptr);
    }
      
    opt_variable_ptr = ((double *) opt_variable + idim);
    value += 0.5 * (*opt_variable_ptr) * product_entry;

    // now add linear term

    conjugate_arg_ptr = ((double *) conjugate_arg + idim);
    value -= (*conjugate_arg_ptr) * (*opt_variable_ptr);

    // now log term

    scaling_ptr = ((double *) scaling + idim);
    value += log(((*opt_variable_ptr) + (*scaling_ptr)) / (*opt_variable_ptr));

  }

  return(value);
}

//    grad = lambda u: -conjugate_arg + precision.dot(u) + (1./(scaling + u) - 1./u)

void barrier_gradient(double *gradient,                   // Gradient vector
		      double *opt_variable,               // Optimization variable
		      double *conjugate_arg,              // Argument to conjugate of Gaussian
		      double *precision,                  // Precision matrix of Gaussian
		      double *scaling,                    // Diagonal scaling matrix for log barrier
		      int ndim)                           // Dimension of conjugate_arg, precision
{
  int idim, jdim;
  double *gradient_ptr;
  double *conjugate_arg_ptr;
  double *opt_variable_ptr;
  double *precision_ptr;
  double *scaling_ptr;
  double product_entry;

  for (idim=0; idim<ndim; idim++) {
    
    gradient_ptr = ((double *) gradient + idim);

    // Compute an entry of precision.dot(conjugate_arg)

    product_entry = 0;
    for (jdim=0; jdim<ndim; jdim++) {
    
      precision_ptr = ((double *) precision + idim * ndim + jdim); // precision is a symmetric matrix
      opt_variable_ptr = ((double *) opt_variable + jdim);
      product_entry += (*precision_ptr) * (*opt_variable_ptr);
    }
      
    opt_variable_ptr = ((double *) opt_variable + idim);
    *gradient_ptr = product_entry;

    // now add linear term

    conjugate_arg_ptr = ((double *) conjugate_arg + idim);
    *gradient_ptr -= (*conjugate_arg_ptr);

    // now log term

    scaling_ptr = ((double *) scaling + idim);
    *gradient_ptr = *gradient_ptr + 1. / ((*opt_variable_ptr) + (*scaling_ptr)) - 1. / (*opt_variable_ptr);
  }

}

double barrier_gradient_step(double *gradient,                   // Gradient vector
			     double *opt_variable,               // Optimization variable
			     double *opt_proposed,               // Proposed value of optimization variable
			     double *conjugate_arg,              // Argument to conjugate of Gaussian
			     double *precision,                  // Precision matrix of Gaussian
			     double *scaling,                    // Diagonal scaling matrix for log barrier
			     double step,                        // Step size for gradient step
			     int ndim)                           // Dimension of conjugate_arg, precision
{
  int idim;
  double *gradient_ptr;
  double *opt_variable_ptr;
  double *opt_proposed_ptr;

  for (idim=0; idim<ndim; idim++) {
    opt_variable_ptr = ((double *) opt_variable + idim);
    opt_proposed_ptr = ((double *) opt_proposed + idim);
    gradient_ptr = ((double *) gradient + idim);
    *opt_proposed_ptr = (*opt_variable_ptr) - step * (*gradient_ptr);
  }

  return barrier_objective(opt_proposed,
			   conjugate_arg,         
			   precision,             
			   scaling,               
			   ndim);
}

double barrier_solve(double *gradient,                   // Gradient vector
		     double *opt_variable,               // Optimization variable
		     double *opt_proposed,               // New value of optimization variable
		     double *conjugate_arg,              // Argument to conjugate of Gaussian
		     double *precision,                  // Precision matrix of Gaussian
		     double *scaling,                    // Diagonal scaling matrix for log barrier
		     int ndim,                           // Dimension of conjugate_arg, precision
		     int max_iter,                       // Maximum number of iterations
		     double value_tol,                   // Tolerance for convergence based on value
		     double initial_step)                // Initial step size
{
  int iter, idim, istep;
  double *gradient_ptr;
  double *opt_variable_ptr;
  double *opt_proposed_ptr;
  double current_value = barrier_objective(opt_variable,
					   conjugate_arg,         
					   precision,             
					   scaling,               
					   ndim);
  double proposed_value;
  double step = initial_step;
  int any_negative;

  for (iter=0; iter<max_iter; iter++) {
    
    // Compute the gradient

    barrier_gradient(gradient,    
		     opt_variable,
		     conjugate_arg,
		     precision,    
		     scaling,
		     ndim);
    
    // Find a valid step size

    istep = 0;
    while (1) {
      any_negative = 0;
      for (idim=0; idim<ndim; idim++) {
	opt_variable_ptr = ((double *) opt_variable + idim);
	gradient_ptr = ((double *) gradient + idim);
	if ((*opt_variable) - step * (*gradient_ptr) < 0) {
	  any_negative += 1;
	}
      }
      if (any_negative == 0) {
	break;
      }
      step = step * 0.5;
      istep++;
      if (istep == 50) {
	break;            // Terminate eventually -- but this will be a failure.
	                  // Should mean that opt_variable is in feasible.
      }
    }

    // Find a step size that is a descent

    istep = 0;
    while (1) {
      proposed_value = barrier_gradient_step(gradient,
					     opt_variable,
					     opt_proposed,
					     conjugate_arg,         
					     precision,             
					     scaling,               
					     step,
					     ndim);
      if (proposed_value < current_value) {
	for (idim=0; idim<ndim; idim++) {
	  opt_variable_ptr = ((double *) opt_variable + idim);
	  opt_proposed_ptr = ((double *) opt_proposed + idim);
	  *opt_variable_ptr = *opt_proposed_ptr;
	}
	break;
      }
      step = step * 0.5;
      istep++;
      if (istep == 50) {
	break;            // Terminate eventually -- this will mean no descent.
	                  // We've solved the problem
      }
    }

    if (fabs(current_value - proposed_value) < value_tol * fabs(current_value)) {
      current_value = proposed_value;
      break;
    }
    current_value = proposed_value;
    
  }

  return(current_value);

}

