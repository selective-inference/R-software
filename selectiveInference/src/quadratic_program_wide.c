#include <math.h> // for fabs
#include <stdio.h>

// Find an approximate row of \hat{Sigma}^{-1}

// Solves a dual version of problem (4) of https://arxiv.org/pdf/1306.3171.pdf

// Dual problem: \text{min}_{\theta} 1/2 \|X\theta\|^2 - l^T\theta + \mu \|\theta\|_1
// where l is `linear_func` below

// This is the "negative" of the problem as in https://gist.github.com/jonathan-taylor/07774d209173f8bc4e42aa37712339bf
// Therefore we don't have to negate the answer to get theta.
// Update one coordinate 

// Throughout X is a design matrix

double objective_wide(double *X_theta_ptr,     /* Fitted values */
		      double *linear_func_ptr, /* Linear term in objective */
		      int *ever_active_ptr,    /* Ever active set: 0-based */ 
		      int *nactive_ptr,        /* Size of ever active set */
		      int ncase,               /* how many rows in X */
		      int nfeature,            /* how many columns in X */
		      double bound,            /* Lagrange multipler for \ell_1 */
		      double *theta_ptr)       /* current value */
{
  int icase, iactive;
  double value = 0;
  double *X_theta_ptr_tmp = X_theta_ptr;
  double *linear_func_ptr_tmp = linear_func_ptr;
  double *theta_ptr_tmp;
  int *active_feature_ptr;
  int active_feature;
  int nactive = *nactive_ptr;

  // The term \|X\theta\|^2_2/n, with n=ncase

  for (icase=0; icase<ncase; icase++) {

    X_theta_ptr_tmp = ((double *) X_theta_ptr + icase);
    value += (*(X_theta_ptr_tmp)) * (*(X_theta_ptr_tmp));

  }

  value *= 0.5 / ncase;

  for (iactive=0; iactive<nactive; iactive++) {

    // The linear term in the objective

    active_feature_ptr = ((int *) ever_active_ptr + iactive);
    active_feature = *active_feature_ptr - 1;          // Ever-active is 1-based

    // fprintf(stderr, "active %d\n", active_feature);

    theta_ptr_tmp = ((double *) theta_ptr + active_feature);
    linear_func_ptr_tmp = ((double *) linear_func_ptr + active_feature);
    value += (*linear_func_ptr_tmp) * (*theta_ptr_tmp); 

    // The \ell_1 term

    value += bound * fabs((*theta_ptr_tmp));

  }
  
  return(value);
}

// Compute, update and return one coordinate of the gradient of \|X\theta\|^2_2/2n

double compute_gradient_coord(double *gradient_ptr,        /* Gradient -- one coordinate will be updated if needed */
			      double *X_theta_ptr,         /* Current fitted values */
			      double *X_ptr,               /* A design matrix */
			      double *linear_func_ptr,     /* Linear term in objective */   
			      int *need_update_ptr,        /* Which coordinates need to be updated? */
			      int coord,                   /* Coordinate we are trying to update */
			      int ncase)                   /* How many rows in X */
{

  double *gradient_ptr_tmp;
  double *X_theta_ptr_tmp;
  double *X_ptr_tmp;  
  double *linear_func_ptr_tmp;  
  int *need_update_ptr_tmp;
  int icase;
  double value = 0.;

  need_update_ptr_tmp = ((int *) need_update_ptr + coord);
  // Check if this coordinate needs updating
  if (*need_update_ptr_tmp == 1) {

    // fprintf(stderr, "compute grad %d %d\n", ncase, coord);

    for (icase=0; icase<ncase; icase++) {
      X_ptr_tmp = ((double *) X_ptr + coord * ncase + icase);
      X_theta_ptr_tmp = ((double *) X_theta_ptr + icase);
      value += (*X_ptr_tmp) * (*X_theta_ptr_tmp);
    }

    value /= ncase;

    linear_func_ptr_tmp = ((double *) linear_func_ptr + coord);
    value += *linear_func_ptr_tmp;

    gradient_ptr_tmp = ((double *) gradient_ptr + coord);
    *gradient_ptr_tmp = value;

    *need_update_ptr_tmp = 0;
  }
  else {
    gradient_ptr_tmp = ((double *) gradient_ptr + coord);
    value = *gradient_ptr_tmp;
  }

  // fprintf(stderr, "grad %d, %f\n", coord, value);
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
    // fprintf(stderr, "active check %d\n", active_var-1);
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

void update_gradient_wide(double *gradient_ptr,     /* X^TX/n times theta */
			  double *X_theta_ptr,      /* Current fitted values */
			  double *X_ptr,            /* A design matrix */
			  double *linear_func_ptr,     /* Linear term in objective */   
			  int *need_update_ptr,     /* Which coordinates need to be updated? */
			  int nfeature,             /* how many columns in X */
			  int ncase)                /* how many rows in X */
{
  int ifeature;

  for (ifeature=0; ifeature<nfeature; ifeature++) {
    compute_gradient_coord(gradient_ptr, X_theta_ptr, X_ptr, linear_func_ptr, need_update_ptr, ifeature, ncase);
  }

}

int check_KKT_wide(double *theta_ptr,        /* current theta */
		   double *gradient_ptr,     /* X^TX/n times theta */
		   double *X_theta_ptr,      /* Current fitted values */
		   double *X_ptr,            /* A design matrix */
		   double *linear_func_ptr,     /* Linear term in objective */   
		   int *need_update_ptr,     /* Which coordinates need to be updated? */
		   int nfeature,             /* how many columns in X */
		   int ncase,                /* how many rows in X */
   		   double bound,             /* Lagrange multipler for \ell_1 */
		   double tol)               /* precision for checking KKT conditions */        
{
  // First check inactive

  int ifeature;
  double *theta_ptr_tmp;
  double gradient;

  for (ifeature=0; ifeature<nfeature; ifeature++) {

    theta_ptr_tmp = ((double *) theta_ptr + ifeature);

    // Compute this coordinate of the gradient

    // fprintf(stderr, "%d \n", ifeature);
    gradient = compute_gradient_coord(gradient_ptr, X_theta_ptr, X_ptr, linear_func_ptr, need_update_ptr, ifeature, ncase);

    if (*theta_ptr_tmp != 0) { // these coordinates of gradients should be equal to -bound
      // fprintf(stderr, "WTF\n");
      if ((*theta_ptr_tmp > 0) &&  (fabs(gradient + bound) > tol * bound)) {
	// fprintf(stderr, "WTF2\n");
	return(0);
      }
      else if ((*theta_ptr_tmp < 0) && (fabs(gradient - bound) > tol * bound)) {
	// fprintf(stderr, "WTF3\n");
	return(0);
      }
    }
    else {
      if (fabs(gradient) > (1. + tol) * bound) {
	// fprintf(stderr, "inactive\n");
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
			     double *X_theta_ptr,         /* X\theta -- fitted values */
			     int *need_update_ptr,        /* Whether a gradient coordinate needs update or not */
			     int ncase,                    /* How many rows in X */
			     int nfeature,                    /* How many rows in X */
			     double bound,                /* feasibility parameter */
			     double *theta_ptr,               /* current value */
			     int coord,                   /* which coordinate to update: 0-based */
			     int is_active)               /* Is this coord in ever_active */     
{

  double delta;
  double linear_term = 0;
  double value = 0;
  double old_value;
  double *X_ptr_tmp;
  double *X_theta_ptr_tmp;
  int *need_update_ptr_tmp;
  double *theta_ptr_tmp;
  int ifeature, icase;

  double *quadratic_ptr = ((double *) X_diag_ptr + coord);
  double quadratic_term = *quadratic_ptr;

  linear_term = compute_gradient_coord(gradient_ptr, X_theta_ptr, X_ptr, linear_func_ptr, need_update_ptr, coord, ncase);

  theta_ptr_tmp = ((double *) theta_ptr + coord);
  old_value = *theta_ptr_tmp;

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

    // Set the update_gradient_ptr to 1

    need_update_ptr_tmp = need_update_ptr;
    for (ifeature=0; ifeature<nfeature; ifeature++) {
      *need_update_ptr_tmp = 1;
      need_update_ptr_tmp++;
    }

    // Update X\theta

    delta = value - old_value;

    for (icase=0; icase<ncase; icase++) {
      X_theta_ptr_tmp = ((double *) X_theta_ptr + icase);
      X_ptr_tmp = ((double *) X_ptr + coord * ncase + icase);
      *X_theta_ptr_tmp = (*X_theta_ptr_tmp) + delta * (*X_ptr_tmp);
    }

    theta_ptr_tmp = ((double *) theta_ptr + coord);
    *theta_ptr_tmp = value;

  }

  return(value);

}

int solve_wide(double *X_ptr,              /* A design matrix */
	       double *X_theta_ptr,        /* Fitted values   */
	       double *linear_func_ptr,    /* Linear term in objective */
	       double *X_diag_ptr,         /* Diagonal entry of covariance matrix */
	       double *gradient_ptr,       /* X times theta */
	       int *need_update_ptr,       /* Keeps track of updated gradient coords */
	       int *ever_active_ptr,       /* Ever active set: 1-based */ 
	       int *nactive_ptr,           /* Size of ever active set */
	       int ncase,                  /* How many rows in X */
	       int nfeature,               /* How many columns in X */
	       double bound,               /* feasibility parameter */
	       double *theta_ptr,          /* current value */
	       int maxiter,                /* max number of iterations */
	       double kkt_tol,             /* precision for checking KKT conditions */
	       double objective_tol,       /* precision for checking relative decrease in objective value */
	       int max_active)             /* Upper limit for size of active set -- otherwise break */ 
{

  int iter = 0;
  int ifeature = 0;
  int iactive = 0;
  int *active_ptr;

  int check_objective = 1;

  double old_value, new_value; 

  // fprintf(stderr, "here1\n");

  if (check_objective) {
    // fprintf(stderr, "here9 %d %d\n", ncase, nfeature);

    old_value = objective_wide(X_theta_ptr,
			       linear_func_ptr,
			       ever_active_ptr,
			       nactive_ptr,
			       ncase,
			       nfeature,
			       bound,
			       theta_ptr);

  }


  // fprintf(stderr, "here2\n");

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
			    X_theta_ptr,
			    need_update_ptr,
			    ncase,
			    nfeature,
			    bound,
			    theta_ptr,
			    *active_ptr - 1,   // Ever-active is 1-based
			    1);
      active_ptr++;
    }

  // fprintf(stderr, "here3\n");

    // Check KKT

    if (check_KKT_wide(theta_ptr,
		       gradient_ptr,
		       X_theta_ptr,
		       X_ptr,
		       linear_func_ptr,
		       need_update_ptr,
		       nfeature,
		       ncase,
		       bound,
		       kkt_tol) == 1) {
      fprintf(stderr, "break1\n");
      break;
    }

    // Update all variables 

    for (ifeature=0; ifeature<nfeature; ifeature++) {

      // fprintf(stderr, "updating %d\n", ifeature);
      update_one_coord_wide(X_ptr,
			    linear_func_ptr,
			    X_diag_ptr,
			    gradient_ptr,
			    ever_active_ptr,
			    nactive_ptr,
			    X_theta_ptr,
			    need_update_ptr,
			    ncase,
			    nfeature,
			    bound,
			    theta_ptr,
			    ifeature,
			    0);
    }

    // fprintf(stderr, "here6\n");

    // Check KKT

    if (check_KKT_wide(theta_ptr,
		       gradient_ptr,
		       X_theta_ptr,
		       X_ptr,
		       linear_func_ptr,
		       need_update_ptr,
		       nfeature,
		       ncase,
		       bound,
		       kkt_tol) == 1) {
      fprintf(stderr, "break2\n");
      break;
    }
					  
    // Check size of active set

    // fprintf(stderr, "here18 %d\n", *nactive_ptr);

    if (*nactive_ptr >= max_active) {
      fprintf(stderr, "break3\n");
      break;
    }

    // Check relative decrease of objective

    // fprintf(stderr, "here7\n");

    if (check_objective) {
      new_value = objective_wide(X_theta_ptr,
				 linear_func_ptr,
				 ever_active_ptr,
				 nactive_ptr,
				 ncase,
				 nfeature,
				 bound,
				 theta_ptr);

      // fprintf(stderr, "here8\n");

      if ((fabs(old_value - new_value) < objective_tol * fabs(new_value)) && (iter > 0)) {
	fprintf(stderr, "break5 %f %f %f %d\n", old_value, new_value, objective_tol, iter);
	break;
      }
      old_value = new_value;
    }

      // fprintf(stderr, "here10\n");

  }
  return(iter);
}

