#include <math.h> // for fabs

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
