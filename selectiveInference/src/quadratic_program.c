#include <math.h> // for fabs

// Find an approximate row of \hat{nndef}^{-1}

// Solves a dual version of problem (4) of https://arxiv.org/pdf/1306.3171.pdf

// Dual problem: \text{min}_{\theta} 1/2 \theta^T \Sigma \theta - l^T\theta + \mu \|\theta\|_1
// where l is `linear_func` below

// This is the "negative" of the problem as in https://gist.github.com/jonathan-taylor/07774d209173f8bc4e42aa37712339bf
// Therefore we don't have to negate the answer to get theta.
// Update one coordinate 

double objective_qp(double *nndef_ptr,       /* A non-negative definite matrix */
		    double *linear_func_ptr, /* Linear term in objective */
		    int *ever_active_ptr,    /* Ever active set: 0-based */ 
		    int *nactive_ptr,        /* Size of ever active set */
		    int nrow,                /* how many rows in nndef */
		    double bound,            /* Lagrange multipler for \ell_1 */
		    double *theta_ptr)           /* current value */
{
  int irow, icol;
  double value = 0;
  double *nndef_ptr_tmp = nndef_ptr;
  double *linear_func_ptr_tmp = linear_func_ptr;
  double *theta_row_ptr, *theta_col_ptr;
  int *active_row_ptr, *active_col_ptr;
  int active_row, active_col;
  int nactive = *nactive_ptr;

  theta_row_ptr = theta_ptr;
  theta_col_ptr = theta_ptr;

  for (irow=0; irow<nactive; irow++) {

    active_row_ptr = ((int *) ever_active_ptr + irow);
    active_row = *active_row_ptr - 1;          // Ever-active is 1-based
    theta_row_ptr = ((double *) theta_ptr + active_row);

    for (icol=0; icol<nactive; icol++) {
      
      active_col_ptr = ((int *) ever_active_ptr + icol);
      active_col = *active_col_ptr - 1;          // Ever-active is 1-based
      theta_col_ptr = ((double *) theta_ptr + active_col);

      nndef_ptr_tmp = ((double *) nndef_ptr + nrow * active_col + active_row); // Matrices are column-major order

      value += 0.5 * (*nndef_ptr_tmp) * (*theta_row_ptr) * (*theta_col_ptr);
    }
    value += bound * fabs((*theta_row_ptr)); // the \ell_1 term

    // The linear term in the objective

    linear_func_ptr_tmp = ((double *) linear_func_ptr + active_row);
    value += (*linear_func_ptr_tmp) * (*theta_row_ptr); 

  }
  
  return(value);
}

// Ever-active is 1-based
// coord is 0-based
int update_ever_active_qp(int coord,
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

int check_KKT_qp(double *theta_ptr,        /* current theta */
		 double *gradient_ptr, /* nndef times theta + linear_func */
		 int nfeature,             /* how many features in nndef */
		 double bound,         /* Lagrange multipler for \ell_1 */
		 double tol)           /* precision for checking KKT conditions */        
{
  // First check inactive

  int ifeature;
  double *theta_ptr_tmp, *gradient_ptr_tmp;
  double gradient;

  for (ifeature=0; ifeature<nfeature; ifeature++) {
    theta_ptr_tmp = ((double *) theta_ptr + ifeature);
    gradient_ptr_tmp = ((double *) gradient_ptr + ifeature);

    // Compute this coordinate of the gradient

    gradient = *gradient_ptr_tmp;

    if (*theta_ptr_tmp != 0) { // these coordinates of gradients should be equal to -bound
      if ((*theta_ptr_tmp > 0) &&  (fabs(gradient + bound) > tol * bound)) {
	return(0);
      }
      else if ((*theta_ptr_tmp < 0) && (fabs(gradient - bound) > tol * bound)) {
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

int check_KKT_qp_active(int *ever_active_ptr,           /* Ever active set: 0-based */ 
		        int *nactive_ptr,               /* Size of ever active set */
			double *theta_ptr,        /* current theta */
			double *gradient_ptr, /* nndef times theta + linear_func */
			int nfeature,             /* how many features in nndef */
			double bound,         /* Lagrange multipler for \ell_1 */
			double tol)           /* precision for checking KKT conditions */        
{
  // First check inactive

  int iactive;
  double *theta_ptr_tmp;
  double gradient;
  double *gradient_ptr_tmp;
  int nactive = *nactive_ptr;
  int active_feature;
  int *active_feature_ptr;

  for (iactive=0; iactive<nactive; iactive++) {

    active_feature_ptr = ((int *) ever_active_ptr + iactive);
    active_feature = *active_feature_ptr - 1;          // Ever-active is 1-based
    theta_ptr_tmp = ((double *) theta_ptr + active_feature);

    gradient_ptr_tmp = ((double *) gradient_ptr + active_feature);

    // Compute this coordinate of the gradient

    gradient = *gradient_ptr_tmp;

    if (*theta_ptr_tmp != 0) { // these coordinates of gradients should be equal to -bound

      if ((*theta_ptr_tmp > 0) &&  (fabs(gradient + bound) > tol * bound)) {
	return(0);
      }
      else if ((*theta_ptr_tmp < 0) && (fabs(gradient - bound) > tol * bound)) {
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


double update_one_coord_qp(double *nndef_ptr,           /* A non-negative definite matrix */
			   double *linear_func_ptr,     /* Linear term in objective */
			   double *nndef_diag_ptr,      /* Diagonal of nndef */
			   double *gradient_ptr,        /* nndef times theta + linear_func */
			   int *ever_active_ptr,        /* Ever active set: 1-based */ 
			   int *nactive_ptr,            /* Size of ever active set */
			   int nfeature,                    /* How many features in nndef */
			   double bound,                /* feasibility parameter */
			   double *theta_ptr,               /* current value */
			   int coord,                   /* which coordinate to update: 0-based */
			   int is_active)               /* Is this coord in ever_active */     
{

  double delta;
  double linear_term = 0;
  double value = 0;
  double old_value;
  double *nndef_ptr_tmp;
  double *gradient_ptr_tmp;
  double *theta_ptr_tmp;
  int icol = 0;

  double *quadratic_ptr = ((double *) nndef_diag_ptr + coord);
  double quadratic_term = *quadratic_ptr;

  gradient_ptr_tmp = ((double *) gradient_ptr + coord);
  linear_term = *gradient_ptr_tmp;

  theta_ptr_tmp = ((double *) theta_ptr + coord);
  old_value = *theta_ptr_tmp;

  // The coord entry of gradient_ptr term has a diagonal term in it:
  // nndef[coord, coord] * theta[coord]
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
    update_ever_active_qp(coord, ever_active_ptr, nactive_ptr);
  }

  // Update the linear term

  if (fabs(old_value - value) > 1.e-6 * (fabs(value) + fabs(old_value))) { 

    delta = value - old_value;
    nndef_ptr_tmp = ((double *) nndef_ptr + coord * nfeature);
    gradient_ptr_tmp = ((double *) gradient_ptr);

    for (icol=0; icol<nfeature; icol++) {
      *gradient_ptr_tmp = *gradient_ptr_tmp + delta * (*nndef_ptr_tmp);
      gradient_ptr_tmp += 1;
      nndef_ptr_tmp += 1;
    }

    theta_ptr_tmp = ((double *) theta_ptr + coord);
    *theta_ptr_tmp = value;

  }

  return(value);

}

int solve_qp(double *nndef_ptr,          /* A non-negative definite matrix */
	     double *linear_func_ptr,    /* Linear term in objective */
	     double *nndef_diag_ptr,     /* Diagonal of nndef */
	     double *gradient_ptr,       /* nndef times theta */
	     int *ever_active_ptr,       /* Ever active set: 1-based */ 
	     int *nactive_ptr,           /* Size of ever active set */
	     int nfeature,               /* How many features in nndef */
	     double bound,               /* feasibility parameter */
	     double *theta,              /* current value */
	     double *theta_old,          /* previous value */
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

  int iter_active;
  int niter_active=5;

  double old_value, new_value; 
  double norm_diff = 1.;
  double norm_last = 1.;
  double delta;
  double *theta_ptr, *theta_old_ptr;

  if (objective_stop) {

    old_value = objective_qp(nndef_ptr,
			     linear_func_ptr,
			     ever_active_ptr,
			     nactive_ptr,
			     nfeature,
			     bound,
			     theta);

  }


  for (iter=0; iter<maxiter; iter++) {

    // Update the active variables first -- do this niter_active times

    for (iter_active=0; iter_active<niter_active; iter_active++) { 

        active_ptr = (int *) ever_active_ptr;
        for (iactive=0; iactive < *nactive_ptr; iactive++) {

	  update_one_coord_qp(nndef_ptr,
			      linear_func_ptr,
			      nndef_diag_ptr,
			      gradient_ptr,
			      ever_active_ptr,
			      nactive_ptr,
			      nfeature,
			      bound,
			      theta,
			      *active_ptr - 1,   // Ever-active is 1-based
			      1);
	  active_ptr++;
	}

	// Check KKT of active subproblem

	if (check_KKT_qp_active(ever_active_ptr,
				nactive_ptr,
				theta,
				gradient_ptr,
				nfeature,
				bound,
				kkt_tol) == 1) {
	  break;
	}
    }

    // Check KKT

    if (kkt_stop) {
      if (check_KKT_qp(theta, 
		       gradient_ptr,
		       nfeature,
		       bound,
		       kkt_tol) == 1) {
	break;
      }
    }

    // Update all variables 

    for (ifeature=0; ifeature<nfeature; ifeature++) {

      update_one_coord_qp(nndef_ptr,
			  linear_func_ptr,
			  nndef_diag_ptr,
			  gradient_ptr,
			  ever_active_ptr,
			  nactive_ptr,
			  nfeature,
			  bound,
			  theta,
			  ifeature,
			  0);
    }

    // Check KKT

    if (kkt_stop) {
      if (check_KKT_qp(theta, 
		       gradient_ptr,
		       nfeature,
		       bound,
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
	  theta_old_ptr = ((double *) theta_old + ifeature);
	  theta_ptr = ((double *) theta + ifeature);
	  delta = (*theta_ptr - *theta_old_ptr);
	  norm_diff += delta * delta;
	  norm_last += (*theta_ptr) * (*theta_ptr);
	  *theta_old_ptr = (*theta_ptr);
	}
	norm_diff = sqrt(norm_diff);
	norm_last = sqrt(norm_last);
	
	if (norm_diff < parameter_tol * norm_last) {
	  break;
	}
      }

      // Check relative decrease of objective

      if (objective_stop) {
	new_value = objective_qp(nndef_ptr,
				 linear_func_ptr,
				 ever_active_ptr,
				 nactive_ptr,
				 nfeature,
				 bound,
				 theta);

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

