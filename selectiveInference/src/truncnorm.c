#include <R.h>
#include <Rmath.h>

// Take a Gibbs hit and run step along a given direction

// Assumes the covariance is identity

void gibbs_step(double *state,     /* state has law N(0,I) constrained to polyhedral set \{y:Ay \leq b\}*/ 
		double *direction, /* direction we will take Gibbs step */
		double *U,         /* A %*% state - b */
		double *alpha,     /* A %*% direction */
		int nconstraint,   /* number of rows of A */
		int nstate)        /* dimension of state */
{
               
  int istate;
  double value = 0;

  /* Compute V=\eta^Ty */

  for (istate = 0; istate < nstate; istate++) {
    value += direction[istate] * state[istate];
  }

  /* Compute upper and lower bounds */

  double lower_bound = -1e12;
  double upper_bound = 1e12;
  double bound_val = 0;
  double tol=1.e-7;
  int iconstraint;

  for (iconstraint = 0; iconstraint < nconstraint; iconstraint++) {

    bound_val = -U[iconstraint] / alpha[iconstraint] + value;

    if ((alpha[iconstraint] > tol) &&
	(bound_val < upper_bound)) {
	upper_bound = bound_val;
    }
    else if ((alpha[iconstraint] < -tol) &&
	  (bound_val > lower_bound)) {
	lower_bound = bound_val;
    }

  }

  /* Ensure constraints are satisfied */

  if (lower_bound > value) {
    lower_bound = value - tol;
  }
  else if (upper_bound < value) {
    upper_bound = value + tol;
  }

  /* Check to see if constraints are satisfied */

  /* if (lower_bound > upper_bound) {
    
     }*/

  /* Now, take a step */

  double tnorm; /* the 1D gaussian variable */
  double cdfU, cdfL, unif; /* temp variables */

  if (upper_bound < -10) {

    /* use Exp approximation */
    /* the approximation is that */
    /* Z | lower_bound < Z < upper_bound */
    /* is fabs(upper_bound) * (upper_bound - Z) = E approx Exp(1) */
    /* so Z = upper_bound - E / fabs(upper_bound) */
    /* and the truncation of the exponential is */
    /* E < fabs(upper_bound - lower_bound) * fabs(upper_bound) = D */
    
    /* this has distribution function (1 - exp(-x)) / (1 - exp(-D)) */
    /* so to draw from this distribution */
    /* we set E = - log(1 - U * (1 - exp(-D))) where U is Unif(0,1) */
    /* and Z (= tnorm below) is as stated */

    unif = runif(0., 1.) * (1 - exp(-fabs((lower_bound - upper_bound) * upper_bound)));
    tnorm = (upper_bound + log(1 - unif) / fabs(upper_bound));
  }
  else if (lower_bound > 10) {

    /* here Z = lower_bound + E / fabs(lower_bound) (though lower_bound is positive) */
    /* and D = fabs((upper_bound - lower_bound) * lower_bound) */

    unif = runif(0., 1.) * (1 - exp(-fabs((upper_bound - lower_bound) * lower_bound)));
    tnorm = (lower_bound - log(1 - unif) / lower_bound);
  }
  else if (lower_bound < 0) {
    cdfL = pnorm(lower_bound, 0., 1., 1, 0); 
    cdfU = pnorm(upper_bound, 0., 1., 1, 0); 
    unif = runif(0., 1.) * (cdfU - cdfL) + cdfL; 
    if (unif < 0.5) {
      tnorm = qnorm(unif, 0., 1., 1, 0); 
    }
    else {
      tnorm = -qnorm(1-unif, 0., 1., 1, 0); 
    }
  }
  else {
    cdfL = pnorm(-lower_bound, 0., 1., 1, 0); 
    cdfU = pnorm(-upper_bound, 0., 1., 1, 0); 
    unif = runif(0., 1.) * (cdfL - cdfU) + cdfU;
    if (unif < 0.5) {
      tnorm = -qnorm(unif, 0., 1., 1, 0); 
    }
    else {
      tnorm = qnorm(1-unif, 0., 1., 1, 0);
    }
  }

  /* Now update the state and U */

  double delta = tnorm - value;

  for (istate = 0; istate < nstate; istate++) {
    state[istate] += delta * direction[istate];
  }
  for (iconstraint = 0; iconstraint < nconstraint; iconstraint++) {
    U[iconstraint] += delta * alpha[iconstraint] ;
  }

  /* End of gibbs_step */

}

void sample_truncnorm_white(double *state,      /* state has law N(0,I) constrained to polyhedral set \{y:Ay \leq b\}*/ 
			    double *U,          /* A %*% state - b */
			    double *directions, /* possible steps for sampler to take */
                                                /* assumed to be stored as list of columns of dimension nstate */
			                        /* has shape (nstate, ndirection) */
			    double *alphas,     /* The matrix A %*% directions */
      			                        /* has shape (nconstraint, ndirection) */
			    double *output,     /* array in which to store samples */
                                                /* assumed will stored as list of vectors of dimension nstate */
                                                /* has shape (nstate, ndraw) */  
			    int *pnconstraint,  /* number of rows of A */
			    int *pndirection,   /* the possible number of directions to choose from */
                                                /* `directions` should have size nstate*ndirection */
			    int *pnstate,       /* dimension of state */
			    int *pburnin,       /* number of burnin steps */
			    int *pndraw)        /* total number of samples to return */
{

  int iter_count;
  int which_direction;

  int nconstraint = *pnconstraint;
  int ndirection = *pndirection;
  int nstate = *pnstate;
  int burnin = *pburnin;
  int ndraw = *pndraw;

  double *direction, *alpha;

  for (iter_count = 0; iter_count < burnin + ndraw; iter_count++) {
  
    which_direction = (int) floor(runif(0., 1.) * ndirection); 
    direction = ((double *) directions) + nstate * which_direction; 
    alpha = ((double *) alphas) + nconstraint * which_direction; 

    /* take a step, which implicitly updates `state` and `U` */

    gibbs_step(state,
               direction,
	       U,
	       alpha,
	       nconstraint,
	       nstate);

    /* Store result if after burnin */

    int istate;
    if (iter_count >= burnin) {
      for (istate = 0; istate < nstate; istate++) {
	*output = state[istate];
	output++;
      }
    }
  }

}

