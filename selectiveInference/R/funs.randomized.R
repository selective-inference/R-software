# Functions to fit and "infer" about parameters in the
# randomized LASSO
#
# min 1/2 || y - \beta_0 - X \beta ||_2^2 + \lambda || \beta ||_1 - \omega^T\beta + \frac{\epsilon}{2} \|\beta\|^2_2

randomizedLASSO = function(X, 
                           y, 
                           lam, 
                           noise_scale, 
                           ridge_term, 
                           noise_type=c('gaussian', 'laplace'),
                           max_iter=100,        # how many iterations for each optimization problem
                           kkt_tol=1.e-4,       # tolerance for the KKT conditions
                           parameter_tol=1.e-8, # tolerance for relative convergence of parameter
                           objective_tol=1.e-8, # tolerance for relative decrease in objective
                           objective_stop=FALSE,
                           kkt_stop=TRUE,
                           param_stop=TRUE)
{

    n = nrow(X); p = ncol(X)
    			
    noise_type = match.arg(noise_type)

    if (noise_scale > 0) {
        if (noise_type == 'gaussian') {
            perturb_ = rnorm(p) * noise_scale
        }
        else if (noise_type == 'laplace') {
            perturb_ = rexp(p) * (2 * rbinom(p, 1, 0.5) - 1) * noise_scale
        }
    } else {
        perturb_ = rep(0, p)
    }

    lam = as.numeric(lam)
    if (length(lam) == 1) {
       lam = rep(lam, p)
    }

    if (length(lam) != p) {
       stop("Lagrange parameter should be single float or of length ncol(X)")
    }    

    soln = rep(0, p)
    Xsoln = rep(0, n)
    linear_func = (- t(X) %*% y - perturb_) / n
    gradient = 1. * linear_func
    ever_active = rep(0, p)
    nactive = as.integer(0)

    result = solve_QP_wide(X,                  # design matrix
    	                   lam / n,                # vector of Lagrange multipliers
		           ridge_term / n,     # ridge_term 
                           max_iter, 
                           soln, 
                           linear_func, 
                           gradient, 
                           Xsoln,
                           ever_active, 
                           nactive, 
                           kkt_tol, 
                           objective_tol, 
			   parameter_tol,
                           p,
		           objective_stop,     # objective_stop
			   kkt_stop,           # kkt_stop
			   param_stop)         # param_stop

    
    sign_soln = sign(result$soln)

    unpenalized = lam == 0
    active = !unpenalized * (sign_soln != 0)
    inactive = !unpenzlied * (sign_soln == 0)

    unpenalized_set = which(unpenalized)
    active_set = which(active)
    inactive_set = which(inactive)

    coef_term = t(X) %*% X[,c(unpenalized_set,  # the coefficients
                              active_set)]
    coef_term = coef_term %*% diag(c(rep(1, sum(unpenalized)), sign_soln[active]))  # coefficients are non-negative
    coef_term[active,] = coef_term[active,] + ridge_term * diag(rep(1, sum(active)))  # ridge term

    subgrad_term = cbind(matrix(0, sum(inactive), sum(active) + sum(unpenalized)),
                         diag(rep(1, sum(inactive))))
    linear_term = rbind(coef_term,
                        subgrad_term)

    offset_term = rep(0, p)
    offset_term[active] = lam[active] * sign_soln[active]

    

    list(active_set = active_set,
         inactive_set = inactive_set,
         unpenalized_set = unpenalized_set,
         sign_soln = sign_soln)
         
         
    return(result)
}
