#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

void multiply_by_2(double *X, int nval);

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
	     int param_stop);            /* Break based on parameter convergence? */

int check_KKT_qp(double *theta,        /* current theta */
		 double *gradient_ptr, /* nndef times theta + linear_func */
		 int nrow,             /* how many rows in nndef */
		 double bound,         /* Lagrange multipler for \ell_1 */
		 double tol);          /* precision for checking KKT conditions */        

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
	       int param_stop);            /* Break based on parameter convergence? */

int check_KKT_wide(double *theta_ptr,        /* current theta */
		   double *gradient_ptr,     /* X^TX/ncase times theta + linear_func*/
		   double *X_theta_ptr,      /* Current fitted values */
		   double *X_ptr,            /* Sqrt of non-neg def matrix -- X^TX/ncase = nndef */
		   double *linear_func_ptr,  /* Linear term in objective */   
		   int *need_update_ptr,     /* Which coordinates need to be updated? */
		   int nfeature,             /* how many columns in X */
		   int ncase,                /* how many rows in X */
		   double *bound_ptr,        /* Lagrange multiplers for \ell_1 */
		   double ridge_term,        /* Ridge / ENet term */
		   double tol);              /* precision for checking KKT conditions */        
  
void update_gradient_wide(double *gradient_ptr,     /* X^TX/ncase times theta + linear_func */
			  double *X_theta_ptr,      /* Current fitted values */
			  double *X_ptr,            /* Sqrt of non-neg def matrix -- X^TX/ncase = nndef */
			  double *linear_func_ptr,  /* Linear term in objective */   
			  int *need_update_ptr,     /* Which coordinates need to be updated? */
			  int nfeature,             /* how many columns in X */
			  int ncase);               /* how many rows in X */


#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */
