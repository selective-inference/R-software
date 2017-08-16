#include <stdio.h>
#include <math.h> // for fabs

// Find an approximate row of \hat{Sigma}^{-1}

// Solves a dual version of problem (4) of https://arxiv.org/pdf/1306.3171.pdf

// Dual problem: \text{min}_{\theta} 1/2 \theta^T \Sigma \theta - e_i^T\theta + \mu \|\theta\|_1

// This is the "negative" of the problem as in https://gist.github.com/jonathan-taylor/07774d209173f8bc4e42aa37712339bf
// Therefore we don't have to negate the answer to get theta.
// Update one coordinate 

double objective(double *Sigma,       /* A covariance matrix: X^TX/n */
		 int nrow,            /* how many rows in Sigma */
		 int row,             /* which row: 0-based */
		 double bound,        /* Lagrange multipler for \ell_1 */
		 double *theta)       /* current value */
{
  int irow, icol;
  double value = 0;
  double *Sigma_ptr = Sigma;
  double *theta_row_ptr, *theta_col_ptr;

  theta_row_ptr = theta;
  theta_col_ptr = theta;

  for (irow=0; irow<nrow; irow++) {
    double *theta_col_ptr = theta;
    for (icol=0; icol<nrow; icol++) {
      value += 0.5 * (*Sigma_ptr) * (*theta_row_ptr) * (*theta_col_ptr);
      Sigma_ptr++;
      theta_col_ptr++;
    }
    if (irow == row) {
      value -= (*theta_row_ptr); // the elementary basis vector term
    }

    value = value + bound * fabs((*theta_row_ptr)); // the \ell_1 term
    theta_row_ptr++;
  }

  return(value);
}


double update_one_coord(double *Sigma,           /* A covariance matrix: X^TX/n */
                        double *Sigma_diag,      /* Diagonal entries of Sigma */
                        double *Sigma_theta,     /* Sigma times theta */
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

  if (fabs(old_value - value) > 1.e-6 * (fabs(value) + fabs(old_value))) { // Update the linear term

    delta = value - old_value;
    Sigma_ptr = ((double *) Sigma + coord * nrow);
    Sigma_theta_ptr = ((double *) Sigma_theta);

    for (icol=0; icol<nrow; icol++) {
      *Sigma_theta_ptr = *Sigma_theta_ptr + delta * (*Sigma_ptr);
      Sigma_theta_ptr += 1;
      Sigma_ptr += 1;
    }

    double before = objective(Sigma,
			      nrow,
			      row,
			      bound,
			      theta);
  fprintf(stderr, "before %f\n", before);

  theta_ptr = ((double *) theta + coord);
  *theta_ptr = value;

  double after = objective(Sigma,
			   nrow,
			   row,
			   bound,
			   theta);

  fprintf(stderr, "after %f\n", after);
  if (after > before) {
    fprintf(stderr, "not a descent step!!!!!!!!!!!!!!!!!!!!!\n");
  }


  }

    Sigma_ptr = ((double *) Sigma + coord * nrow);
    Sigma_theta_ptr = ((double *) Sigma_theta);
    for (icol=0; icol<nrow; icol++) {
      
      Sigma_theta_ptr += 1;
      Sigma_ptr += 1;
    }


  return(value);

}

void find_one_row(double *Sigma,          /* A covariance matrix: X^TX/n */
                  double *Sigma_diag,     /* Diagonal entry of covariance matrix */
                  double *Sigma_theta,    /* Sigma times theta */
                  int *nrow_ptr,          /* How many rows in Sigma */
		  double *bound_ptr,      /* feasibility parameter */
                  double *theta,          /* current value */
                  int *maxiter_ptr,       /* how many iterations */
                  int *row_ptr)         /* which coordinate to update: 0-based */
{

  int maxiter = *maxiter_ptr;
  int iter = 0;
  int icoord = 0;
  int row = *row_ptr;
  double bound = *bound_ptr;
  int nrow = *nrow_ptr;

  double old_value = objective(Sigma,
			       nrow,
			       row,
			       bound,
			       theta);
  double new_value; 
  double tol=1.e-10;

  for (iter=0; iter<maxiter; iter++) {

    // Update the diagonal first

    update_one_coord(Sigma,
		     Sigma_diag,
		     Sigma_theta,
		     nrow,
		     bound,
		     theta,
		     row,
		     row);

    for (icoord=0; icoord<nrow; icoord++) {

      update_one_coord(Sigma,
		       Sigma_diag,
		       Sigma_theta,
		       nrow,
		       bound,
		       theta,
		       row,
		       icoord);
    }

    new_value = objective(Sigma,
			  nrow,
			  row,
			  bound,
			  theta);

    if (((old_value - new_value) < tol * fabs(new_value)) && (iter > 5)) {
      break;
    }

    fprintf(stderr, "%f %f value\n", old_value, new_value);
    old_value = new_value;
  }

  *nrow_ptr = iter-1;
}

