# Functions to fit and "infer" about parameters in the
# randomized LASSO
#
# min 1/2 || y - \beta_0 - X \beta ||_2^2 + \lambda || \beta ||_1 - \omega^T\beta + \frac{\epsilon}{2} \|\beta\|^2_2

randomizedLasso = function(X, 
                           y, 
                           lam, 
                           noise_scale=NULL, 
                           ridge_term=NULL, 
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
    			
    mean_diag = mean(apply(X^2, 2, sum))

    # default ridge term

    if (is.null(ridge_term)) {
        ridge_term = sqrt(mean_diag) * sd(y) / sqrt(n)
    }

    # default noise level

    if (is.null(noise_scale)) {
        noise_scale = 0.5 * sd(y) * sqrt(mean_diag)
    }

    print(c(noise_scale, ridge_term))
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

    print('here')
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

    print('now')
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
                                              observed_raw+opt_transform$offset_term)
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
                             
conditional_density = function(noise_scale, lasso_soln) {
  
  active_set = lasso_soln$active_set
  observed_raw = lasso_soln$observed_raw
  opt_linear = lasso_soln$optimization_transform$linear_term
  opt_offset =  lasso_soln$optimization_transform$offset_term
  observed_opt_state = lasso_soln$observed_opt_state
  
  nactive = length(active_set)
  B = opt_linear[,1:nactive]
  beta_offset = opt_offset
  p=length(observed_opt_state)
  if (nactive<p){
    beta_offset = beta_offset+(opt_linear[,(nactive+1):p] %*% observed_opt_state[(nactive+1):p])
  }
  opt_transform = list(linear_term=B, 
                       offset_term = beta_offset)
  reduced_B = chol(t(B) %*% B)
  beta_offset = beta_offset + observed_raw
  reduced_beta_offset = solve(t(reduced_B)) %*% (t(B) %*% beta_offset)
  
  log_condl_optimization_density = function(opt_state) {
    if  (sum(opt_state < 0) > 0) {
      return(-Inf)
    }
    D = log_density_gaussian_conditional_(noise_scale,
                                          reduced_B,
                                          as.matrix(opt_state),
                                          reduced_beta_offset)
    return(D)
  }
  lasso_soln$log_optimization_density = log_condl_optimization_density
  lasso_soln$observed_opt_state = observed_opt_state[1:nactive]
  lasso_soln$optimization_transform = opt_transform
  return(lasso_soln)
}

randomizedLassoInf = function(X, 
                              y, 
                              lam, 
                              sigma=NULL, 
                              noise_scale=NULL, 
                              ridge_term=NULL, 
                              condition_subgrad=TRUE, 
                              level=0.9) {

  n = nrow(X)
  p = ncol(X)
  lasso_soln = randomizedLasso(X, y, lam, noise_scale=noise_scale, ridge_term=ridge_term)
  active_set = lasso_soln$active_set
  inactive_set = lasso_soln$inactive_set
  nactive = length(active_set)
  
  if (condition_subgrad==TRUE){
    lasso_soln=conditional_density(noise_scale,lasso_soln)
  } 
    
  dim = length(lasso_soln$observed_opt_state)
  print(paste("chain dim", dim))
  S = sample_opt_variables(lasso_soln, jump_scale=rep(1/sqrt(n), dim), nsample=10000)
  opt_samples = S$samples[2001:10000,]
  print(paste("dim opt samples", toString(dim(opt_samples))))
  
  X_E = X[, active_set]
  X_minusE = X[, inactive_set]

  # if no sigma given, use OLS estimate

  if (is.null(sigma)) {
        lm_y = lm(y ~ X_E - 1)
        sigma = sqrt(sum(resid(lm_y)^2) / lm_y$df.resid)
  }        
  print(c(sigma, 'sigma'))
  target_cov = solve(t(X_E) %*% X_E)*sigma^2
  cov_target_internal = rbind(target_cov, matrix(0, nrow=p-nactive, ncol=nactive))
  observed_target = solve(t(X_E) %*% X_E) %*% t(X_E) %*% y
  observed_internal = c(observed_target, t(X_minusE) %*% (y-X_E%*% observed_target))
  internal_transform = lasso_soln$internal_transform
  opt_transform = lasso_soln$optimization_transform
  observed_raw = lasso_soln$observed_raw
  
  pvalues = rep(0, nactive)
  ci = matrix(0, nactive, 2)
  for (i in 1:nactive){
    target_transform = linear_decomposition(observed_target[i], 
                                            observed_internal, 
                                            target_cov[i,i], 
                                            cov_target_internal[,i],
                                            internal_transform)
    target_sample = rnorm(nrow(opt_samples)) * sqrt(target_cov[i,i])
    
    pivot = function(candidate){
      weights = importance_weight(noise_scale,
                                  t(as.matrix(target_sample)) + candidate,
                                  t(opt_samples),
                                  opt_transform,
                                  target_transform,
                                  observed_raw)
      return(mean((target_sample+candidate<observed_target[i])*weights)/mean(weights))
    }
    rootU = function(candidate){
      return (pivot(observed_target[i]+candidate)-(1-level)/2)
    }
    rootL = function(candidate){
      return (pivot(observed_target[i]+candidate)-(1+level)/2)
    }
    pvalues[i] = pivot(0)
    line_min = -20*sd(target_sample)
    line_max = 20*sd(target_sample)
    if (rootU(line_min)*rootU(line_max)<0){
      ci[i,2] = uniroot(rootU, c(line_min, line_max))$root+observed_target[i]
    } else{
      print("non inv u")
      ci[i,2]=line_max
    }
    if (rootL(line_min)*rootL(line_max)<0){
      ci[i,1] = uniroot(rootL, c(line_min, line_max))$root+observed_target[i]
    } else{
      print("non inv u")
      ci[i,1] = line_min
    }
  }
  return(list(active_set=active_set, pvalues=pvalues, ci=ci))
}
