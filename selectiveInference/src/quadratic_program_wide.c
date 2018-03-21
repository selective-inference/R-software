#include <math.h> // for fabs

// Find an approximate row of \hat{nndef}^{-1}

// Solves a dual version of problem (4) of https://arxiv.org/pdf/1306.3171.pdf

// Dual problem: \text{min}_{\theta} 1/2 \|X\theta\|^2/n + l^T\theta + \mu \|\theta\|_1 + \frac{\epsilon}{2} \|\theta\|^2_2
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
		      double *bound_ptr,       /* Lagrange multiplers for \ell_1 */
		      double ridge_term,       /* Ridge / ENet term */
		      double *theta_ptr)       /* current value */
{
  int icase, iactive;
  double value = 0;
  double *bound_ptr_tmp;
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

    theta_ptr_tmp = ((double *) theta_ptr + active_feature);
    linear_func_ptr_tmp = ((double *) linear_func_ptr + active_feature);
    value += (*linear_func_ptr_tmp) * (*theta_ptr_tmp); 

    // The \ell_1 term

    bound_ptr_tmp = ((double *) bound_ptr + active_feature);
    value += (*bound_ptr_tmp) * fabs((*theta_ptr_tmp));

    // The ridge term

    value += 0.5 * ridge_term * (*theta_ptr_tmp) * (*theta_ptr_tmp);
  }
  
  return(value);
}

// Compute, update and return one coordinate of the gradient of \|X\theta\|^2_2/2n

double compute_gradient_coord(double *gradient_ptr,        /* Gradient -- one coordinate will be updated if needed */
			      double *X_theta_ptr,         /* Current fitted values */
			      double *X_ptr,               /* Sqrt of non-neg def matrix -- X^TX/ncase = nndef */
			      double *linear_func_ptr,     /* Linear term in objective */   
			      int *need_update_ptr,        /* Which coordinates need to be updated? */
			      int coord,                   /* Coordinate we are trying to update */
			      int ncase,                   /* How many rows in X */
			      int nfeature)                /* How many columns in X */
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

    for (icase=0; icase<ncase; icase++) {
#ifdef COLUMN_MAJOR_ORDER
      X_ptr_tmp = ((double *) X_ptr + coord * ncase + icase);
#else
      X_ptr_tmp = ((double *) X_ptr + icase * nfeature + coord);
#endif
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

void update_gradient_wide(double *gradient_ptr,     /* X^TX/ncase times theta + linear_func */
			  double *X_theta_ptr,      /* Current fitted values */
			  double *X_ptr,            /* Sqrt of non-neg def matrix -- X^TX/ncase = nndef */
			  double *linear_func_ptr,  /* Linear term in objective */   
			  int *need_update_ptr,     /* Which coordinates need to be updated? */
			  int ncase,             /* how many rows in X */
			  int nfeature)                /* how many columns in X */
{
  int ifeature;

  for (ifeature=0; ifeature<nfeature; ifeature++) {
    compute_gradient_coord(gradient_ptr, X_theta_ptr, X_ptr, linear_func_ptr, need_update_ptr, ifeature, ncase, nfeature);
  }

}

int check_KKT_wide(double *theta_ptr,        /* current theta */
		   double *gradient_ptr,     /* X^TX/ncase times theta + linear_func*/
		   double *X_theta_ptr,      /* Current fitted values */
		   double *X_ptr,            /* Sqrt of non-neg def matrix -- X^TX/ncase = nndef */
		   double *linear_func_ptr,  /* Linear term in objective */   
		   int *need_update_ptr,     /* Which coordinates need to be updated? */
		   int ncase,                /* how many rows in X */
		   int nfeature,             /* how many columns in X */
		   double *bound_ptr,        /* Lagrange multiplers for \ell_1 */
		   double ridge_term,        /* Ridge / ENet term */
		   double tol)               /* precision for checking KKT conditions */        
{
  // First check inactive

  int ifeature;
  double *theta_ptr_tmp;
  double *bound_ptr_tmp;
  double bound;
  double gradient;

  for (ifeature=0; ifeature<nfeature; ifeature++) {

    theta_ptr_tmp = ((double *) theta_ptr + ifeature);
    bound_ptr_tmp = ((double *) bound_ptr + ifeature);
    bound = *bound_ptr_tmp;

    // Compute this coordinate of the gradient

    gradient = compute_gradient_coord(gradient_ptr, X_theta_ptr, X_ptr, linear_func_ptr, need_update_ptr, ifeature, ncase, nfeature);

    if ((*theta_ptr_tmp != 0) && (bound != 0)) { // these coordinates of gradients should be equal to -bound
      
      if ((*theta_ptr_tmp > 0) &&  (fabs(gradient + ridge_term * (*theta_ptr_tmp) + bound) > tol * bound)) {
	return(0);
      }
      else if ((*theta_ptr_tmp < 0) && (fabs(gradient + ridge_term * (*theta_ptr_tmp) - bound) > tol * bound)) {
	return(0);
      }
      
    }
    else if (bound != 0) {
      if (fabs(gradient) > (1. + tol) * bound) {
	return(0);
      }
    }
  }

  return(1);
}

int check_KKT_wide_active(int *ever_active_ptr,           /* Ever active set: 0-based */ 
			  int *nactive_ptr,               /* Size of ever active set */
			  double *theta_ptr,              /* current theta */
			  double *gradient_ptr,           /* X^TX/ncase times theta + linear_func*/
			  double *X_theta_ptr,            /* Current fitted values */
			  double *X_ptr,                  /* Sqrt of non-neg def matrix -- X^TX/ncase = nndef */
			  double *linear_func_ptr,        /* Linear term in objective */   
			  int *need_update_ptr,           /* Which coordinates need to be updated? */
			  int ncase,                      /* how many rows in X */
			  int nfeature,                   /* how many columns in X */
			  double *bound_ptr,              /* Lagrange multipliers for \ell_1 */
			  double ridge_term,              /* Ridge / ENet term */
			  double tol)                     /* precision for checking KKT conditions */        
{
  // First check inactive

  int iactive;
  double *theta_ptr_tmp;
  double *bound_ptr_tmp;
  double bound;
  double gradient;
  int nactive = *nactive_ptr;
  int active_feature;
  int *active_feature_ptr;

  for (iactive=0; iactive<nactive; iactive++) {

    active_feature_ptr = ((int *) ever_active_ptr + iactive);
    active_feature = *active_feature_ptr - 1;          // Ever-active is 1-based

    theta_ptr_tmp = ((double *) theta_ptr + active_feature);
    bound_ptr_tmp = ((double *) bound_ptr + active_feature);
    bound = *bound_ptr_tmp;

    // Compute this coordinate of the gradient

    gradient = compute_gradient_coord(gradient_ptr, X_theta_ptr, X_ptr, linear_func_ptr, need_update_ptr, active_feature, ncase, nfeature);

    if ((*theta_ptr_tmp != 0) && (bound != 0)) { // these coordinates of gradients should be equal to -bound

      if ((*theta_ptr_tmp > 0) &&  (fabs(gradient + ridge_term * (*theta_ptr_tmp) + bound) > tol * bound)) {
	return(0);
      }
      else if ((*theta_ptr_tmp < 0) && (fabs(gradient + ridge_term * (*theta_ptr_tmp) - bound) > tol * bound)) {
	return(0);
      }

    }
    else if (bound != 0) {
      if (fabs(gradient) > (1. + tol) * bound) {
	return(0);
      }
    }
  }

  return(1);
}

double update_one_coord_wide(double *X_ptr,               /* A design matrix*/
			     double *linear_func_ptr,     /* Linear term in objective */
			     double *nndef_diag_ptr,      /* Diagonal of nndef */
			     double *gradient_ptr,        /* X^TX/ncase times theta + linear_func*/
			     int *ever_active_ptr,        /* Ever active set: 1-based */ 
			     int *nactive_ptr,            /* Size of ever active set */
			     double *X_theta_ptr,         /* X\theta -- fitted values */
			     int *need_update_ptr,        /* Whether a gradient coordinate needs update or not */
			     int ncase,                   /* How many rows in X */
			     int nfeature,                /* How many rows in X */
			     double *bound_ptr,           /* Lagrange multipliers */
			     double ridge_term,           /* Ridge / ENet term */
			     double *theta_ptr,           /* current value */
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
  double *bound_ptr_tmp;
  double bound;
  int ifeature, icase;

  double *diagonal_ptr = ((double *) nndef_diag_ptr + coord);
  double diagonal_entry = *diagonal_ptr;

  linear_term = compute_gradient_coord(gradient_ptr, X_theta_ptr, X_ptr, linear_func_ptr, need_update_ptr, coord, ncase, nfeature);

  theta_ptr_tmp = ((double *) theta_ptr + coord);
  old_value = *theta_ptr_tmp;

  bound_ptr_tmp = ((double *) bound_ptr + coord);
  bound = *bound_ptr_tmp;

  // The coord entry of gradient_ptr term has a diagonal term in it:
  // (X^TX)[coord, coord] * theta[coord] / ncase
  // This removes it. 

  linear_term -= diagonal_entry * old_value;

  // Now soft-threshold the coord entry of theta 

  // Objective is t \mapsto (q+eps)/2 * t^2 + l * t + bound |t| + 
  // with q=diagonal_entry and l=linear_term and eps=ridge_Term

  // With a negative linear term, solution should be
  // positive

  if (linear_term < -bound) {
    value = (-linear_term - bound) / (diagonal_entry + ridge_term);
  }
  else if (linear_term > bound) {
    value = -(linear_term - bound) / (diagonal_entry + ridge_term);
  }

  // Add to active set if necessary

  if ((is_active == 0) && (value != 0)) {
    update_ever_active_wide(coord, ever_active_ptr, nactive_ptr);
  }

  // Update X\theta if changed

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
#ifdef COLUMN_MAJOR_ORDER
      X_ptr_tmp = ((double *) X_ptr + coord * ncase + icase);
#else
      X_ptr_tmp = ((double *) X_ptr + icase * nfeature + coord);
#endif
      *X_theta_ptr_tmp = (*X_theta_ptr_tmp) + delta * (*X_ptr_tmp);
    }

    theta_ptr_tmp = ((double *) theta_ptr + coord);
    *theta_ptr_tmp = value;

    theta_ptr_tmp = theta_ptr;
    for (ifeature=0; ifeature<nfeature; ifeature++) {
      theta_ptr_tmp++;
    }

  }

  return(value);

}

int solve_wide(double *X_ptr,              /* Sqrt of non-neg def matrix -- X^TX/ncase = nndef */
	       double *X_theta_ptr,        /* Fitted values   */
	       double *linear_func_ptr,    /* Linear term in objective */
	       double *nndef_diag_ptr,     /* Diagonal entries of non-neg def matrix */
	       double *gradient_ptr,       /* X^TX/ncase times theta + linear_func*/
	       int *need_update_ptr,       /* Keeps track of updated gradient coords */
	       int *ever_active_ptr,       /* Ever active set: 1-based */ 
	       int *nactive_ptr,           /* Size of ever active set */
	       int ncase,                  /* How many rows in X */
	       int nfeature,               /* How many columns in X */
	       double *bound_ptr,          /* Lagrange multipliers */
	       double ridge_term,          /* Ridge / ENet term */
	       double *theta_ptr,          /* current value */
	       double *theta_old_ptr,      /* previous value */
	       int maxiter,                /* max number of iterations */
	       double kkt_tol,             /* precision for checking KKT conditions */
	       double objective_tol,       /* precision for checking relative decrease in objective value */
	       double parameter_tol,       /* precision for checking relative convergence of parameter */
	       int max_active,             /* Upper limit for size of active set -- otherwise break */ 
	       int kkt_stop,               /* Break based on KKT? */
	       int objective_stop,         /* Break based on convergence of objective value? */
	       int param_stop)             /* Break based on parameter convergence? */
{

  int iter = 0;
  int iter_old = 1;
  int ifeature = 0;
  int iactive = 0;
  int *active_ptr;

  double old_value, new_value; 
  int niter_active = 5;
  int iter_active;

  double norm_diff = 1.;
  double norm_last = 1.;
  double delta;
  double *theta_ptr_tmp, *theta_old_ptr_tmp;

  new_value = objective_wide(X_theta_ptr,
			     linear_func_ptr,
			     ever_active_ptr,
			     nactive_ptr,
			     ncase,
			     nfeature,
			     bound_ptr,
			     ridge_term,
			     theta_ptr);

  old_value = new_value + 2000000000; // hack for a big number... should do better


  for (iter=0; iter<maxiter; iter++) {

    // Update the active variables first -- do this niter_active times

    for (iter_active=0; iter_active<niter_active; iter_active++) { 

        active_ptr = (int *) ever_active_ptr;
        for (iactive=0; iactive < *nactive_ptr; iactive++) {

          update_one_coord_wide(X_ptr,
    			        linear_func_ptr,
				nndef_diag_ptr,
				gradient_ptr,
				ever_active_ptr,
				nactive_ptr,
				X_theta_ptr,
				need_update_ptr,
				ncase,
				nfeature,
				bound_ptr,
				ridge_term,
				theta_ptr,
				*active_ptr - 1,   // Ever-active is 1-based
				1);
	  active_ptr++;
	}

	// Check KKT of active subproblem

	if (check_KKT_wide_active(ever_active_ptr,
				  nactive_ptr,
				  theta_ptr,
				  gradient_ptr,
				  X_theta_ptr,
				  X_ptr,
				  linear_func_ptr,
				  need_update_ptr,
				  ncase,
				  nfeature,
				  bound_ptr,
				  ridge_term,
				  kkt_tol) == 1) {
	  break;
	}

    }

    // Check KKT

    if (kkt_stop) { 
      if (check_KKT_wide(theta_ptr,
			 gradient_ptr,
			 X_theta_ptr,
			 X_ptr,
			 linear_func_ptr,
			 need_update_ptr,
			 ncase,
			 nfeature,
			 bound_ptr,
			 ridge_term,
			 kkt_tol) == 1) {
	break;
      }
    }

    // Update all variables 

    for (ifeature=0; ifeature<nfeature; ifeature++) {

      update_one_coord_wide(X_ptr,
			    linear_func_ptr,
			    nndef_diag_ptr,
			    gradient_ptr,
			    ever_active_ptr,
			    nactive_ptr,
			    X_theta_ptr,
			    need_update_ptr,
			    ncase,
			    nfeature,
			    bound_ptr,
			    ridge_term,
			    theta_ptr,
			    ifeature,
			    0);
    }

    // Check KKT

    if (kkt_stop) {
      if (check_KKT_wide(theta_ptr,
			 gradient_ptr,
			 X_theta_ptr,
			 X_ptr,
			 linear_func_ptr,
			 need_update_ptr,
			 ncase,
			 nfeature,
			 bound_ptr,
			 ridge_term,
			 kkt_tol) == 1) {
	break;
      }
    }

    if (iter == 2 * iter_old) { // Geometric iterations from Adel's code

      // Check based on norm 

      if (param_stop) {
	iter_old = iter;
	norm_diff = 0;
	norm_last = 0;
	for (ifeature=0; ifeature<nfeature; ifeature++) {
	  theta_old_ptr_tmp = ((double *) theta_old_ptr + ifeature);
	  theta_ptr_tmp = ((double *) theta_ptr + ifeature);
	  delta = (*theta_ptr_tmp - *theta_old_ptr_tmp);
	  norm_diff += delta * delta;
	  norm_last += (*theta_ptr_tmp) * (*theta_ptr_tmp);
	  *theta_old_ptr = (*theta_ptr_tmp);
	}
	norm_diff = sqrt(norm_diff);
	norm_last = sqrt(norm_last);
	
	if (norm_diff < parameter_tol * norm_last) {
	  break;
	}
      }

      // Check relative decrease of objective

      if (objective_stop) {
	new_value = objective_wide(X_theta_ptr,
				   linear_func_ptr,
				   ever_active_ptr,
				   nactive_ptr,
				   ncase,
				   nfeature,
				   bound_ptr,
				   ridge_term,
				   theta_ptr);

	if ((fabs(old_value - new_value) < objective_tol * fabs(new_value)) && (iter > 0)) {
	  break;
	}
	old_value = new_value;
      }
    }

    // Check size of active set

    if (*nactive_ptr >= max_active) {
      break;
    }

  }
  return(iter);
}

