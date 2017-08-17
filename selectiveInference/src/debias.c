#include <stdio.h>
#include <math.h> // for fabs

// Find an approximate row of \hat{Sigma}^{-1}

// Solves a dual version of problem (4) of https://arxiv.org/pdf/1306.3171.pdf

// Dual problem: \text{min}_{\theta} 1/2 \theta^T \Sigma \theta - e_i^T\theta + \mu \|\theta\|_1

// This is the "negative" of the problem as in https://gist.github.com/jonathan-taylor/07774d209173f8bc4e42aa37712339bf
// Therefore we don't have to negate the answer to get theta.
// Update one coordinate 

double objective(double *Sigma,       /* A covariance matrix: X^TX/n */
		 int *ever_active,    /* Ever active set: 0-based */ 
		 int *nactive_ptr,    /* Size of ever active set */
		 int nrow,            /* how many rows in Sigma */
		 int row,             /* which row: 0-based */
		 double bound,        /* Lagrange multipler for \ell_1 */
		 double *theta)       /* current value */
{
  int irow, icol;
  double value = 0;
  double *Sigma_ptr = Sigma;
  double *theta_row_ptr, *theta_col_ptr;
  int *active_row_ptr, *active_col_ptr;
  int active_row, active_col;
  int nactive = *nactive_ptr;

  theta_row_ptr = theta;
  theta_col_ptr = theta;

  for (irow=0; irow<nactive; irow++) {

    active_row_ptr = ((int *) ever_active + irow);
    active_row = *active_row_ptr;
    theta_row_ptr = ((double *) theta + active_row);

    for (icol=0; icol<nactive; icol++) {
      
      active_col_ptr = ((int *) ever_active + icol);
      active_col = *active_col_ptr;
      theta_col_ptr = ((double *) theta + active_col);
      
      fprintf(stderr, "%d %d \n", active_row, active_col);

      Sigma_ptr = ((double *) Sigma + nrow * active_row + active_col);

      value += 0.5 * (*Sigma_ptr) * (*theta_row_ptr) * (*theta_col_ptr);
    }
    value = value + bound * fabs((*theta_row_ptr)); // the \ell_1 term
  }

  theta_row_ptr = ((double *) theta + row);
  value -= (*theta_row_ptr); // the elementary basis vector term

  return(value);
}

int is_active(int coord,
	      int *nactive_ptr,
	      int *ever_active) {
  int iactive;
  int active_var;
  int nactive = *nactive_ptr;
  int *ever_active_ptr = ever_active;

  for (iactive=0; iactive<nactive; iactive++) {
    active_var = (*ever_active_ptr);
    if (active_var == coord) {
      return(1);
    }
  }
  return(0);
}

int check_KKT(double *theta,       /* current theta */
	      double *Sigma_theta, /* Sigma times theta */
	      int nrow,            /* how many rows in Sigma */
	      int row,             /* which row: 0-based */
	      double bound)        /* Lagrange multipler for \ell_1 */
{
  // First check inactive

  int irow;
  int fail = 0;
  double tol = 1.e-4;
  double *theta_ptr, *Sigma_theta_ptr;
  double gradient;

  for (irow=0; irow<nrow; irow++) {
    theta_ptr = ((double *) theta + irow);
    Sigma_theta_ptr = ((double *) Sigma_theta + irow);

    // Compute this coordinate of the gradient

    gradient = *Sigma_theta_ptr;
    if (row == irow) {
      gradient -= 1;
    }

    if (*theta_ptr != 0) { // these coordinates of gradients should be equal to -bound
      if ((*theta_ptr > 0) &&  (fabs(gradient + bound) > (1. + tol) * bound)) {
	fail += 1;
      }
      else if ((*theta_ptr < 0) && (fabs(gradient - bound) > (1. + tol) * bound)) {
	fail += 1;
      }
    }
    else {
      if (fabs(gradient) > (1. + tol) * bound) {
	fail += 1;
      }
    }
  }

  return(fail == 0);
}


double update_one_coord(double *Sigma,           /* A covariance matrix: X^TX/n */
                        double *Sigma_diag,      /* Diagonal entries of Sigma */
                        double *Sigma_theta,     /* Sigma times theta */
			int *ever_active,        /* Ever active set: 0-based */ 
			int *nactive_ptr,        /* Size of ever active set */
			int nrow,                /* How many rows in Sigma */
			double bound,            /* feasibility parameter */
			double *theta,           /* current value */
			int row,                 /* which row: 0-based */
			int coord)               /* which coordinate to update: 0-based */
{

  double delta;
  double linear_term = 0;
  double value = 0;
  double old_value;
  double *Sigma_ptr;
  double *Sigma_theta_ptr;
  double *theta_ptr;
  int icol = 0;

  double *quadratic_ptr = ((double *) Sigma_diag + coord);
  double quadratic_term = *quadratic_ptr;

  int *ever_active_ptr;

  Sigma_theta_ptr = ((double *) Sigma_theta + coord);
  linear_term = *Sigma_theta_ptr;

  theta_ptr = ((double *) theta + coord);
  old_value = *theta_ptr;

  // The coord entry of Sigma_theta term has a diagonal term in it:
  // Sigma[coord, coord] * theta[coord]
  // This removes it. 
  linear_term -= quadratic_term * old_value;

  if (row == coord) {
    linear_term -= 1;
  }

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

  if ((value != 0) && (is_active(coord, ever_active, nactive_ptr) == 0)) {
    ever_active_ptr = ((int *) ever_active + *nactive_ptr);
    *ever_active_ptr = coord;
    *nactive_ptr += 1;
  }

  // Update the linear term

  if (fabs(old_value - value) > 1.e-6 * (fabs(value) + fabs(old_value))) { 

    delta = value - old_value;
    Sigma_ptr = ((double *) Sigma + coord * nrow);
    Sigma_theta_ptr = ((double *) Sigma_theta);

    for (icol=0; icol<nrow; icol++) {
      *Sigma_theta_ptr = *Sigma_theta_ptr + delta * (*Sigma_ptr);
      Sigma_theta_ptr += 1;
      Sigma_ptr += 1;
    }

    theta_ptr = ((double *) theta + coord);
    *theta_ptr = value;

  }

  return(value);

}

void find_one_row_void(double *Sigma,          /* A covariance matrix: X^TX/n */
		       double *Sigma_diag,     /* Diagonal entry of covariance matrix */
		       double *Sigma_theta,    /* Sigma times theta */
		       int *ever_active,       /* Ever active set: 0-based */ 
		       int *nactive_ptr,       /* Size of ever active set */
		       int *nrow_ptr,          /* How many rows in Sigma */
		       double *bound_ptr,      /* feasibility parameter */
		       double *theta,          /* current value */
		       int *maxiter_ptr,       /* how many iterations */
		       int *row_ptr)           /* which coordinate to update: 0-based */
{

  int maxiter = *maxiter_ptr;
  int iter = 0;
  int icoord = 0;
  int row = *row_ptr;
  double bound = *bound_ptr;
  int nrow = *nrow_ptr;

  fprintf(stderr, "starting now\n");

  double old_value = objective(Sigma,
			       ever_active,
			       nactive_ptr,
			       nrow,
			       row,
			       bound,
			       theta);
  double new_value; 
  double tol=1.e-5;

  for (iter=0; iter<maxiter; iter++) {

    // Update the diagonal first

    update_one_coord(Sigma,
		     Sigma_diag,
		     Sigma_theta,
		     ever_active,
		     nactive_ptr,
		     nrow,
		     bound,
		     theta,
		     row,
		     row);

    if (check_KKT(theta,
		  Sigma_theta,
		  nrow,
		  row,
		  bound) == 1) {
      fprintf(stderr, "ending in first KKT check\n");
      break;
    }
					  
    for (icoord=0; icoord<nrow; icoord++) {

      update_one_coord(Sigma,
		       Sigma_diag,
		       Sigma_theta,
		       ever_active,
		       nactive_ptr,
		       nrow,
		       bound,
		       theta,
		       row,
		       icoord);
    }

    if (check_KKT(theta,
		  Sigma_theta,
		  nrow,
		  row,
		  bound) == 1) {
      fprintf(stderr, "ending in second KKT check\n");
      break;
    }
					  
    new_value = objective(Sigma,
			  ever_active,
			  nactive_ptr,
			  nrow,
			  row,
			  bound,
			  theta);

    if (((old_value - new_value) < tol * fabs(new_value)) && (iter > 0)) {
      fprintf(stderr, "ending in objective value check\n");
      break;
    }

    old_value = new_value;
  }

  *nrow_ptr = iter-1;
}

