#include <math.h> // for fabs

// Find an approximate row of \hat{Sigma}^{-1}

// Solves a dual version of problem (4) of https://arxiv.org/pdf/1306.3171.pdf

// Dual problem: \text{min}_{\theta} 1/2 \|X\theta\|^2 - l^T\theta + \mu \|\theta\|_1
// where l is `linear_func` below

// This is the "negative" of the problem as in https://gist.github.com/jonathan-taylor/07774d209173f8bc4e42aa37712339bf
// Therefore we don't have to negate the answer to get theta.
// Update one coordinate 

// Throughout X is a design matrix

double objective_wide(double *X_ptr,           /* A design matrix */
		      double *X_theta_ptr,     /* Fitted values */
		      double *linear_func_ptr, /* Linear term in objective */
		      int *ever_active_ptr,    /* Ever active set: 0-based */ 
		      int *nactive_ptr,        /* Size of ever active set */
		      int nrow,                /* how many rows in X */
		      int ncol,                /* how many columns in X */
		      double bound,            /* Lagrange multipler for \ell_1 */
		      double *theta)           /* current value */
{
  int irow, icol;
  double value = 0;
  double *X_theta_ptr_tmp = X_theta_ptr;
  double *linear_func_ptr_tmp = linear_func_ptr;
  double *theta_row_ptr, *theta_col_ptr;
  int *active_row_ptr, *active_col_ptr;
  int active_row, active_col;
  int nactive = *nactive_ptr;

  theta_row_ptr = theta;
  theta_col_ptr = theta;

  double entry = 0; // An entry of X\theta

  // The term \|X\theta\|^2_2/nrow

  for (irow=0; irow<nrow; irow++) {

    X_theta_ptr = ((double *) X_theta_ptr + irow);
    value += (*(X_theta_ptr)) * (*(X_theta_ptr));

  }

  for (irow=0; irow<nactive; irow++) {

    // The linear term in the objective

    active_row_ptr = ((int *) ever_active_ptr + irow);
    active_row = *active_row_ptr - 1;          // Ever-active is 1-based
    theta_row_ptr = ((double *) theta + active_row);

    linear_func_ptr_tmp = ((double *) linear_func_ptr + active_row);
    value += (*linear_func_ptr_tmp) * (*theta_row_ptr); 

    // The \ell_1 term

    value += bound * fabs((*theta_row_ptr));

  }
  
  return(value);
}

// Ever-active is 1-based
// coord is 0-based
int update_ever_active_wide(int coord,
			  int *ever_active_ptr,
			  int *nactive_ptr) {
  int iactive;
  int active_var;
  int nactive = *nactive_ptr;
  int *ever_active_ptr_tmp = ever_active_ptr;

  for (iactive=0; iactive<nactive; iactive++) {
    ever_active_ptr_tmp = ((int *) ever_active_ptr + iactive);
    active_var = *ever_active_ptr_tmp;
    if (active_var - 1 == coord) {          // Ever-active is 1-based
      return(1);
    }
  }
  
  // If we haven't returned yet, this means the coord was not in 
  // ever_active.

  // Add it to the active set and increment the 
  // number of active variables

  ever_active_ptr_tmp = ((int *) ever_active_ptr + *nactive_ptr);
  *ever_active_ptr_tmp = coord + 1;         // Ever-active is 1-based
  *nactive_ptr += 1;

  return(0);
}

int check_KKT_wide(double *theta,        /* current theta */
		 double *gradient_ptr, /* X^TX/n times theta */
		 int nrow,             /* how many rows in X */
		 double bound,         /* Lagrange multipler for \ell_1 */
		 double tol)           /* precision for checking KKT conditions */        
{
  // First check inactive

  int irow;
  double *theta_ptr, *gradient_ptr_tmp;
  double gradient;

  for (irow=0; irow<nrow; irow++) {
    theta_ptr = ((double *) theta + irow);
    gradient_ptr_tmp = ((double *) gradient_ptr + irow);

    // Compute this coordinate of the gradient

    gradient = *gradient_ptr_tmp;

    if (*theta_ptr != 0) { // these coordinates of gradients should be equal to -bound
      if ((*theta_ptr > 0) &&  (fabs(gradient + bound) > tol * bound)) {
	return(0);
      }
      else if ((*theta_ptr < 0) && (fabs(gradient - bound) > tol * bound)) {
	return(0);
      }
    }
    else {
      if (fabs(gradient) > (1. + tol) * bound) {
	return(0);
      }
    }
  }

  return(1);
}

double update_one_coord_wide(double *X_ptr,               /* A design matrix*/
			     double *linear_func_ptr,     /* Linear term in objective */
			     double *X_diag_ptr,          /* Diagonal entries of Sigma */
			     double *gradient_ptr,        /* X^TX/n times theta */
			     int *ever_active_ptr,        /* Ever active set: 1-based */ 
			     int *nactive_ptr,            /* Size of ever active set */
			     int nrow,                    /* How many rows in X */
			     int ncol,                    /* How many rows in X */
			     double bound,                /* feasibility parameter */
			     double *theta,               /* current value */
			     int coord,                   /* which coordinate to update: 0-based */
			     int is_active)               /* Is this coord in ever_active */     
{

  double delta;
  double linear_term = 0;
  double value = 0;
  double old_value;
  double *X_ptr_tmp;
  double *gradient_ptr_tmp;
  double *theta_ptr;
  int icol = 0;

  double *quadratic_ptr = ((double *) X_diag_ptr + coord);
  double quadratic_term = *quadratic_ptr;

  gradient_ptr_tmp = ((double *) gradient_ptr + coord);
  linear_term = *gradient_ptr_tmp;

  theta_ptr = ((double *) theta + coord);
  old_value = *theta_ptr;

  // The coord entry of gradient_ptr term has a diagonal term in it:
  // X[coord, coord] * theta[coord]
  // This removes it. 

  linear_term -= quadratic_term * old_value;

  // Now soft-threshold the coord entry of theta 

  // Objective is t \mapsto q/2 * t^2 + l * t + bound |t|
  // with q=quadratic_term and l=linear_term

  // With a negative linear term, solution should be
  // positive

  if (linear_term < -bound) {
    value = (-linear_term - bound) / quadratic_term;
  }
  else if (linear_term > bound) {
    value = -(linear_term - bound) / quadratic_term;
  }

  // Add to active set if necessary

  if ((is_active == 0) && (value != 0)) {
    update_ever_active_wide(coord, ever_active_ptr, nactive_ptr);
  }

  // Update the linear term

  if (fabs(old_value - value) > 1.e-6 * (fabs(value) + fabs(old_value))) { 

    delta = value - old_value;
    X_ptr_tmp = ((double *) X_ptr + coord * nrow);
    gradient_ptr_tmp = ((double *) gradient_ptr);

    for (icol=0; icol<nrow; icol++) {
      *gradient_ptr_tmp = *gradient_ptr_tmp + delta * (*X_ptr_tmp);
      gradient_ptr_tmp += 1;
      X_ptr_tmp += 1;
    }

    theta_ptr = ((double *) theta + coord);
    *theta_ptr = value;

  }

  return(value);

}

int solve_wide(double *X_ptr,              /* A design matrix */
	       double *linear_func_ptr,    /* Linear term in objective */
	       double *X_diag_ptr,         /* Diagonal entry of covariance matrix */
	       double *gradient_ptr,       /* X times theta */
	       int *ever_active_ptr,       /* Ever active set: 1-based */ 
	       int *nactive_ptr,           /* Size of ever active set */
	       int nrow,                   /* How many rows in X */
	       int ncol,                   /* How many columns in X */
	       double bound,               /* feasibility parameter */
	       double *theta,              /* current value */
	       int maxiter,                /* max number of iterations */
	       double kkt_tol,             /* precision for checking KKT conditions */
	       double objective_tol,       /* precision for checking relative decrease in objective value */
	       int max_active)             /* Upper limit for size of active set -- otherwise break */ 
{

  int iter = 0;
  int icoord = 0;
  int iactive = 0;
  int *active_ptr;

  int check_objective = 1;

  double old_value, new_value; 

  if (check_objective) {

    old_value = objective_wide(X_ptr,
			       linear_func_ptr,
			       ever_active_ptr,
			       nactive_ptr,
			       nrow,
			       ncol,
			       bound,
			       theta);

  }


  for (iter=0; iter<maxiter; iter++) {

    // Update the active variables first

    active_ptr = (int *) ever_active_ptr;

    for (iactive=0; iactive < *nactive_ptr; iactive++) {
      update_one_coord_wide(X_ptr,
			    linear_func_ptr,
			    X_diag_ptr,
			    gradient_ptr,
			    ever_active_ptr,
			    nactive_ptr,
			    nrow,
			    ncol,
			    bound,
			    theta,
			    *active_ptr - 1,   // Ever-active is 1-based
			    1);
      active_ptr++;
    }

    // Check KKT

    if (check_KKT_wide(theta,
		     gradient_ptr,
		     nrow,
		     bound,
		     kkt_tol) == 1) {
      break;
    }
					  
    // Update all variables 

    for (icoord=0; icoord<nrow; icoord++) {

      update_one_coord_wide(X_ptr,
			    linear_func_ptr,
			    X_diag_ptr,
			    gradient_ptr,
			    ever_active_ptr,
			    nactive_ptr,
			    nrow,
			    ncol,
			    bound,
			    theta,
			    icoord,
			    0);
    }

    // Check KKT

    if (check_KKT_wide(theta,
		       gradient_ptr,
		       nrow,
		       bound,
		       kkt_tol) == 1) {
      break;
    }
					  
    // Check size of active set

    if (*nactive_ptr >= max_active) {
      break;
    }

    // Check relative decrease of objective

    if (check_objective) {
      new_value = objective_wide(X_ptr,
				 linear_func_ptr,
				 ever_active_ptr,
				 nactive_ptr,
				 nrow,
				 ncol,
				 bound,
				 theta);

      if ((fabs(old_value - new_value) < objective_tol * fabs(new_value)) && (iter > 0)) {
	break;
      }
      old_value = new_value;
    }
  }
  return(iter);
}

