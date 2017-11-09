# Functions to fit and "infer" about parameters in the
# randomized LASSO
#
# min 1/2 || y - \beta_0 - X \beta ||_2^2 + \lambda || \beta ||_1 - \omega^T\beta + \frac{\epsilon}{2} \|\beta\|^2_2

randomizedLasso = function(X, 
                           y, 
                           lam, 
                           family="gaussian",
                           noise_scale=NULL, 
                           ridge_term=NULL, 
                           noise_type=c('gaussian', 'laplace'),
                           max_iter=100,        # how many iterations for each optimization problem
                           kkt_tol=1.e-4,       # tolerance for the KKT conditions
                           parameter_tol=1.e-8, # tolerance for relative convergence of parameter
                           objective_tol=1.e-8, # tolerance for relative decrease in objective
                           objective_stop=FALSE,
                           kkt_stop=TRUE,
                           parameter_stop=TRUE)
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
    
    print(paste("ridge term", ridge_term))
    print(paste("noise scale", noise_scale))
    
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
			                     parameter_stop)         # param_stop
    
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
    observed_subgrad = -n*result$gradient[inactive]
    
    if (length(which(abs(observed_subgrad)>lam[1]))){
      print("subgradient eq not satisfied")
    }

    observed_opt_state = c(observed_unpen, observed_scalings, observed_subgrad)

    # affine transform for optimization variables

    E = c(unpenalized_set, active_set)
    I = inactive_set
    X_E = X[,E]
    X_I = X[,I]
    
    if (family=="binomial"){
      unpen_reg = glm(y~X_E-1, family="binomial")
      unpen_est = unpen_reg$coefficients
      pi_fn = function(beta){
        temp = X_E %*% as.matrix(beta)
        return(as.vector(exp(temp)/(1+exp(temp)))) # n-dimensional
      }
      pi_vec = pi_fn(unpen_est)
      W_E = diag(pi_vec*(1-pi_vec))
    } else if (family=="gaussian"){
      W_E = diag(rep(1,n))
    }
    L_E = t(X) %*% W_E %*% X[,E]
    
    coef_term = L_E

    signs_ = c(rep(1, sum(unpenalized)), sign_soln[active])
    
    coef_term[active,] = coef_term[active,] + ridge_term * diag(rep(1, sum(active)))  # ridge term
  
    if (length(signs_) == 1) {
      coef_term = coef_term * signs_
    } else {
      coef_term = coef_term %*% diag(signs_)  # scaligns are non-negative
    }
    
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
    if (family=="binomial"){
      beta_E = result$soln[active_set]
      observed_raw = observed_raw+t(X)%*%pi_fn(beta_E)-L_E %*% beta_E
    }
    inactive_lam = lam[inactive_set]
    inactive_start = sum(unpenalized) + sum(active)
    active_start = sum(unpenalized)
    
    
    # XXX only for Gaussian so far

    log_optimization_density = function(opt_state) {

        if ((sum(abs(opt_state[(inactive_start + 1):p]) > inactive_lam) > 0) ||
            (sum(opt_state[(active_start + 1):inactive_start] < 0) > 0)) {
            return(-Inf)
        }

        use_C_code = TRUE
        if (!use_C_code) {
            A = opt_transform$linear_term %*% opt_state + observed_raw + opt_transform$offset_term
            D = -apply(A^2, 2, sum) / noise_scale^2
        } else {
            D = log_density_gaussian_conditional_(noise_scale,
                                                  opt_transform$linear_term,
                                                  as.matrix(opt_state),
                                                  observed_raw + opt_transform$offset_term)
        }
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
                observed_raw = observed_raw,
		            noise_scale = noise_scale,
		            soln = result$soln,
		            perturb = perturb_
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

    use_C_code = TRUE
    if (!use_C_code) {
        A = (opt_transform$linear_term %*% opt_sample + 
             target_transform$linear_term %*% target_sample)
        A = apply(A, 2, function(x) {return(x + target_transform$offset_term + opt_transform$offset_term)})
        log_num = -apply(A^2, 2, sum) / noise_scale^2
    } else {
        log_num = log_density_gaussian_(noise_scale,
                                        target_transform$linear_term,
                                        as.matrix(target_sample),
                                        opt_transform$linear_term,
                                        as.matrix(opt_sample),
                                        target_transform$offset_term + opt_transform$offset_term)
    }

    if (!use_C_code) {
        A = opt_transform$linear_term %*% opt_sample 
        A = apply(A, 2, function(x) {return(x + observed_raw + opt_transform$offset_term)})
        log_den = -apply(A^2, 2, sum) / noise_scale^2
    } else {
       log_den = log_density_gaussian_conditional_(noise_scale,
                                                   opt_transform$linear_term,
                                                   as.matrix(opt_sample),
                                                   observed_raw+opt_transform$offset_term)
    }
    W = log_num - log_den
    W = W - max(W)
    return(exp(W))
}

get_mean_cov = function(noise_scale, linear_term, offset_term){
    temp = solve(t(linear_term) %*% linear_term)
    cov = noise_scale^2*temp
    mean = temp %*% t(linear_term) %*% offset_term
    return(list(mean=mean, cov=cov))
}


                             
conditional_density = function(noise_scale, lasso_soln) {
  
  active_set = lasso_soln$active_set
  observed_raw = lasso_soln$observed_raw
  opt_linear = lasso_soln$optimization_transform$linear_term
  opt_offset =  lasso_soln$optimization_transform$offset_term
  observed_opt_state = lasso_soln$observed_opt_state
  
  nactive = length(active_set)
  B = opt_linear[,1:nactive,drop=FALSE]
  beta_offset = opt_offset
  p = length(observed_opt_state)

  if (nactive < p) {
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

    use_C_code = TRUE
    if (!use_C_code) {
        A = reduced_B %*% as.matrix(opt_state) + reduced_beta_offset
        A = apply(A, 2, function(x) {x + reduced_beta_offset})
        log_den = -apply(A^2, 2, sum) / noise_scale^2
    } else {
        log_den = log_density_gaussian_conditional_(noise_scale,
                                                    reduced_B,
                                                    as.matrix(opt_state),
                                                    reduced_beta_offset)
    }
    return(log_den)
  }
  lasso_soln$log_optimization_density = log_condl_optimization_density
  lasso_soln$observed_opt_state = observed_opt_state[1:nactive]
  lasso_soln$optimization_transform = opt_transform
  reduced_opt_transform =list(linear_term = reduced_B, offset_term = reduced_beta_offset)
  return(list(lasso_soln=lasso_soln, 
              reduced_opt_transform = reduced_opt_transform))
}

randomizedLassoInf = function(X, 
                              y, 
                              lam, 
                              family="gaussian",
                              sampler="A",
                              sigma=NULL, 
                              noise_scale=NULL, 
                              ridge_term=NULL, 
                              condition_subgrad=TRUE, 
                              level=0.9,
			                        nsample=10000,
			                        burnin=2000,
                              max_iter=100,        # how many iterations for each optimization problem
                              kkt_tol=1.e-4,       # tolerance for the KKT conditions
                              parameter_tol=1.e-8, # tolerance for relative convergence of parameter
                              objective_tol=1.e-8, # tolerance for relative decrease in objective
                              objective_stop=FALSE,
                              kkt_stop=TRUE,
                              parameter_stop=TRUE)
 {

  n = nrow(X)
  p = ncol(X)
  
  lasso_soln = randomizedLasso(X, 
                               y, 
                               lam, 
                               family=family,
                               noise_scale=noise_scale, 
                               ridge_term=ridge_term,
                               max_iter=max_iter,
                               kkt_tol=kkt_tol,       
                               parameter_tol=parameter_tol,
                               objective_tol=objective_tol,
                               objective_stop=objective_stop,
                               kkt_stop=kkt_stop,
                               parameter_stop=parameter_stop)

  active_set = lasso_soln$active_set
  nactive = length(active_set)
  print(paste("nactive", nactive))
  if (nactive==0){
    return (list(active_set=active_set, pvalues=c(), ci=c()))
  }
  inactive_set = lasso_soln$inactive_set
  
  noise_scale = lasso_soln$noise_scale # set to default value in randomizedLasso
  
  constraints = matrix(0,nactive,2)
  constraints[,2] = Inf
  if (condition_subgrad==TRUE){
   condl_lasso=conditional_density(noise_scale, lasso_soln)
   lasso_soln = condl_lasso$lasso_soln
   cur_opt_transform = condl_lasso$reduced_opt_transform
   } else{
   if (nactive<p){
     subgrad_constraints = matrix(-lam, p-nactive, 2)
     subgrad_constraints[,2]=lam
     constraints = rbind(constraints, subgrad_constraints)
   }
    cur_opt_transform = list(linear_term = lasso_soln$optimization_transform$linear_term,
                             offset_term = lasso_soln$optimization_transform$offset_term+lasso_soln$observed_raw)
  }
    
  ndim = length(lasso_soln$observed_opt_state)
  
  if (sampler =="R"){
    S = sample_opt_variables(lasso_soln, jump_scale=rep(1/sqrt(n), ndim), nsample=nsample)
    opt_samples = as.matrix(S$samples[(burnin+1):nsample,,drop=FALSE])
  } else if (sampler == "A"){
    opt_samples = gaussian_sampler(noise_scale, 
                                 lasso_soln$observed_opt_state, 
                                 cur_opt_transform$linear_term,
                                 cur_opt_transform$offset_term,
                                 constraints,
                                 nsamples=nsample)
    opt_sample = opt_samples[(burnin+1):nsample,]
  }
  
  X_E = X[, active_set]
  X_minusE = X[, inactive_set]

  
  
  if (family=="gaussian"){
    lm_y = lm(y ~ X_E - 1)
    sigma_resid = sqrt(sum(resid(lm_y)^2) / lm_y$df.resid)
    observed_target = lm_y$coefficients
    W_E = diag(rep(1,n))
    observed_internal = c(observed_target, t(X_minusE) %*% (y-X_E%*% observed_target))
  } else if (family=="binomial"){
    glm_y = glm(y~X_E-1)
    sigma_resid = sqrt(sum(resid(glm_y)^2) / glm_y$df.resid)
    observed_target = as.matrix(glm_y$coefficients)
    temp = X_E%*%observed_target
    pi_vec = exp(temp)/(1+exp(temp))
    observed_internal =  c(observed_target, t(X_minusE) %*% (y-pi_vec))
    W_E=diag(as.vector(pi_vec *(1-pi_vec)))
  }
  
  # if no sigma given, use the estimate
  
  if (is.null(sigma)) {
    sigma = sigma_resid
  }        
  
  target_cov = solve(t(X_E) %*% W_E %*% X_E)*sigma^2
  cov_target_internal = rbind(target_cov, matrix(0, nrow=p-nactive, ncol=nactive))
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
    
    
    # changing dimension of density evalutaion

    if ((condition_subgrad == TRUE) & (nactive < p)) {
        target_opt_linear = cbind(target_transform$linear_term, opt_transform$linear_term)
        reduced_target_opt_linear = chol(t(target_opt_linear) %*% target_opt_linear)
        target_linear = reduced_target_opt_linear[,1]
        temp = solve(t(reduced_target_opt_linear)) %*% t(target_opt_linear)
        target_offset = temp %*% target_transform$offset_term
        target_transform = list(linear_term = as.matrix(target_linear), offset_term = target_offset)
        cur_linear = reduced_target_opt_linear[,2:ncol(reduced_target_opt_linear)]
        cur_offset = temp %*% opt_transform$offset_term
        cur_transform = list(linear_term = as.matrix(cur_linear), offset_term = cur_offset)
        raw = target_transform$linear_term * observed_target[i] + target_transform$offset_term
    } else {
        cur_transform = opt_transform
        raw = observed_raw
    }   

    target_sample = rnorm(nrow(as.matrix(opt_samples))) * sqrt(target_cov[i,i])

    pivot = function(candidate){

      weights = importance_weight(noise_scale,
                                  t(as.matrix(target_sample) + candidate),
                                  t(opt_samples),
                                  cur_transform,
                                  target_transform,
                                  raw)
      return(mean((target_sample + candidate < observed_target[i]) * weights)/mean(weights))

    }

    rootU = function(candidate){
      return (pivot(observed_target[i]+candidate)-(1-level)/2)
    }

    rootL = function(candidate){
      return(pivot(observed_target[i]+candidate)-(1+level)/2)
    }

    pvalues[i] = pivot(0)
    line_min = -10*sd(target_sample)
    line_max = 10*sd(target_sample)

    if (rootU(line_min)*rootU(line_max)<0){
      ci[i,2] = uniroot(rootU, c(line_min, line_max))$root+observed_target[i]
    } else{
      ci[i,2]=line_max
    }
    if (rootL(line_min)*rootL(line_max)<0){
      ci[i,1] = uniroot(rootL, c(line_min, line_max))$root+observed_target[i]
    } else{
      ci[i,1] = line_min
    }
  }
  return(list(active_set=active_set, pvalues=pvalues, ci=ci))
}

   
    
    
    



