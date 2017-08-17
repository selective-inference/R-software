#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

void find_one_row_void(double *Sigma,          /* A covariance matrix: X^TX/n */
		       double *Sigma_diag,     /* Diagonal entry of covariance matrix */
		       double *Sigma_theta,    /* Sigma times theta */
		       int *ever_active,       /* Ever active set: 0-based */ 
		       int *nactive_ptr,       /* Size of ever active set */
		       int nrow,          /* How many rows in Sigma */
		       double bound,      /* feasibility parameter */
		       double *theta,          /* current value */
		       int maxiter,       /* how many iterations */
		       int row);          /* which coordinate to update: 0-based */
#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */
