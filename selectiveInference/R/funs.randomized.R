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
    	                   lam / n,            # vector of Lagrange multipliers
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
    active = (!unpenalized) & (sign_soln != 0)
    inactive = (!unpenalized) & (sign_soln == 0)

    unpenalized_set = which(unpenalized)
    active_set = which(active)
    inactive_set = which(inactive)

    # observed opt state

    observed_scalings = abs(result$soln)[active]
    observed_unpen = result$soln[unpenalized]
    observed_subgrad = result$gradient[inactive]

    observed_opt_state = c(observed_unpen, observed_scalings, observed_subgrad)

    # affine transform for optimization variables

    E = c(unpenalized_set, active_set)
    I = inactive_set
    X_E = X[,E]
    X_I = X[,I]
    L_E = t(X) %*% X[,E]

    coef_term = L_E
    coef_term = coef_term %*% diag(c(rep(1, sum(unpenalized)), sign_soln[active]))  # coefficients are non-negative
    coef_term[active,] = coef_term[active,] + ridge_term * diag(rep(1, sum(active)))  # ridge term

    subgrad_term = matrix(0, p, sum(inactive)) # for subgrad
    for (i in 1:sum(inactive)) {
        subgrad_term[inactive_set[i], i] = 1
    }

    linear_term = cbind(coef_term,
                        subgrad_term)

    offset_term = rep(0, p)
    offset_term[active] = lam[active] * sign_soln[active]

    opt_transform = list(linear_term=linear_term,
                         offset_term=offset_term)

    # affine transform for internal (data) variables
    # for now just use parametric in terms of
    # (\bar{\beta}_E, X_{-E}^T(y-X_E\bar{\beta}_E)
    # 
    # we have to reconstruct -X^TY from this pair
    #

    active_term = -L_E                           # for \bar{\beta}_E

    inactive_term = -subgrad_term
    linear_term = cbind(active_term,
                        inactive_term)
    offset_term = rep(0, p)
    internal_transform = list(linear_term = linear_term,
                              offset_term = offset_term)

    # density for sampling optimization variables

    observed_raw = -t(X) %*% y
    inactive_lam = lam[inactive_set]
    inactive_start = sum(unpenalized) + sum(active)
    active_start = sum(unpenalized)

    # XXX only for Gaussian so far

    log_optimization_density = function(opt_state 
                                        ) {


        if ((sum(abs(opt_state[(inactive_start + 1):p]) > inactive_lam) > 0) ||
            (sum(opt_state[(active_start+1):inactive_start] < 0) > 0)) {
            return(-Inf)
        }

        D = log_density_gaussian_conditional_(noise_scale,
                                              opt_transform$linear_term,
                                              as.matrix(opt_state),
                                              observed_raw)
        return(D)
    }

    return(list(active_set = active_set,
                inactive_set = inactive_set,
                unpenalized_set = unpenalized_set,
                sign_soln = sign_soln,
                optimization_transform = opt_transform,
                internal_transform = internal_transform,
                log_optimization_density = log_optimization_density,
		observed_opt_state = observed_opt_state,
                observed_raw = observed_raw
                ))

}

sample_opt_variables = function(randomizedLASSO_obj, jump_scale, nsample=10000) {
      return(MCMC(randomizedLASSO_obj$log_optimization_density, 
                  nsample,
		  randomizedLASSO_obj$observed_opt_state,
                  acc.rate=0.2,
                  scale=jump_scale))
}

# Carry out a linear decompositon of an internal
# representation with respect to a target

# Returns an affine transform into raw coordinates (i.e. \omega or randomization coordinates)

linear_decomposition = function(observed_target,
                                observed_internal,
                                var_target,
                                cov_target_internal,
                                internal_transform) {
    var_target = as.matrix(var_target) 
    if (nrow(var_target) == 1) {
        nuisance = observed_internal - cov_target_internal * observed_target / var_target
        target_linear = internal_transform$linear_term %*% cov_target_internal / var_target[1,1]
    } else {
        nuisance = observed_internal - cov_target_internal %*% solve(var_target) %*% observed_target 
        target_linear = internal_transform$linear_term %*% cov_target_internal %*% solve(var_target)
    }
    target_offset = internal_transform$linear_term %*% nuisance + internal_transform$offset_term
    return(list(linear_term=target_linear,
                offset_term=target_offset))
}

# XXX only for Gaussian so far

importance_weight = function(noise_scale,
                             target_sample,
                             opt_sample,
                             opt_transform,
                             target_transform,
                             observed_raw) {

    log_num = log_density_gaussian_(noise_scale,
                                    target_transform$linear_term,
                                    as.matrix(target_sample),
                                    opt_transform$linear_term,
                                    as.matrix(opt_sample),
                                    target_transform$offset_term + opt_transform$offset_term)

    log_den = log_density_gaussian_conditional_(noise_scale,
                                                opt_transform$linear_term,
                                                as.matrix(opt_sample),
                                                observed_raw+opt_transform$offset_term)
    W = log_num - log_den
    W = W - max(W)
    return(exp(W))
}
                             
