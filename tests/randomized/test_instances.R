#library(devtools)
#devtools::install_github('jonathan-taylor/R-selective/selectiveInference')
library(selectiveInference, lib.loc='/Users/Jelena/anaconda/lib/R/library')


gaussian_instance = function(n, p, s, sigma=1, rho=0, signal=6, X=NA,
                             random_signs=TRUE, scale=TRUE, center=TRUE, seed=NA){
  if (!is.na(seed)){
    set.seed(seed)
  }
  
  if (is.na(X)){
    X = sqrt(1-rho)*matrix(rnorm(n*p),n) + sqrt(rho)*matrix(rep(rnorm(n), p), nrow = n)
    X = scale(X)/sqrt(n)
  }
  beta = rep(0, p)
  if (s>0){
    beta[1:s] = seq(3, 6, length.out=s)
  }
  beta = sample(beta)
  if (random_signs==TRUE & s>0){
    signs = sample(c(-1,1), s, replace = TRUE)
    beta = beta * signs
  }
  y = X %*% beta + rnorm(n)*sigma
  result = list(X=X,y=y,beta=beta)
  return(result)
}

conditional_density = function(noise_scale, lasso_soln){
  
  active_set = lasso_soln$active_set
  inactive_set = lasso_soln$inactive_set
  observed_raw = lasso_soln$observed_raw
  opt_linear = lasso_soln$optimization_transform$linear_term
  opt_offset =  lasso_soln$optimization_transform$offset_term
  observed_opt_state = lasso_soln$observed_opt_state
  
  nactive = length(active_set)
  B = opt_linear[,1:nactive]
  beta_offset = observed_raw+opt_offset
  p=length(observed_opt_state)
  if (nactive<p){
    U=opt_linear[,(nactive+1):p]
    beta_offset =+U %*% observed_opt_state[(nactive+1):p]
  }
  opt_transform = list(linear_term=B, offset_term = beta_offset)
  reduced_B = chol(t(B) %*% B)
  reduced_beta_offset = solve(t(reduced_B)) %*% (t(B) %*% beta_offset)
  
  log_condl_optimization_density = function(opt_state) {
    
    if  (sum(opt_state < 0) > 0) {
      return(-Inf)
    }
    D = selectiveInference:::log_density_gaussian_conditional_(noise_scale,
                                                               reduced_B,
                                                               as.matrix(observed_opt_state[1:nactive]),
                                                               reduced_beta_offset)
    return(D)
  }
  lasso_soln$log_optimization_density = log_condl_optimization_density
  lasso_soln$observed_opt_state = observed_opt_state[1:nactive]
  lasso_soln$optimization_transform = opt_transform
  return(lasso_soln)
}



run_instance = function(n, p, s){
  rho=0.3
  lam=1.
  sigma=1
  data = gaussian_instance(n=n,p=p,s=s, rho=rho, sigma=sigma)
  X=data$X
  print(dim(X))
  y=data$y
  ridge_term=sd(y)/sqrt(n)
  noise_scale = sd(y)/2
  lasso_soln=selectiveInference:::randomizedLASSO(X, y, lam, noise_scale, ridge_term)
  active_set = lasso_soln$active_set
  inactive_set = lasso_soln$inactive_set
  nactive = length(active_set)
  print(paste("nactive", nactive))
  
  lasso_soln = conditional_density(noise_scale, lasso_soln)
  
  S = selectiveInference:::sample_opt_variables(lasso_soln, jump_scale=rep(1/sqrt(n), nactive), nsample=10000)
  opt_samples = S$samples[2001:10000,]
  print(paste("dim opt samples", toString(dim(opt_samples))))
  
  observed_raw = lasso_soln$observed_raw
  opt_linear = lasso_soln$optimization_transform$linear_term
  opt_offset =  lasso_soln$optimization_transform$offset_term
  observed_opt_state = lasso_soln$observed_opt_state
  opt_transform = lasso_soln$optimization_transform
  
  X_E=X[, active_set]
  X_minusE=X[, inactive_set]
  target_cov = solve(t(X_E)%*%X_E)*sigma^2
  cov_target_internal = rbind(target_cov, matrix(0, nrow=p-nactive, ncol=nactive))
  observed_target = solve(t(X_E) %*% X_E) %*% t(X_E) %*% y
  observed_internal = c(observed_target, t(X_minusE) %*% (y-X_E%*% observed_target))
  internal_transform = lasso_soln$internal_transform
  
  pivots = rep(0, nactive)
  for (i in 1:nactive){
    target_transform = selectiveInference:::linear_decomposition(observed_target[i], 
                                                  observed_internal, 
                                                  target_cov[i,i], 
                                                  cov_target_internal[,i],
                                                  internal_transform)
    
    target_sample = rnorm(nrow(opt_samples)) * sqrt(target_cov[i,i])
    
    pivot = function(candidate){
      weights = selectiveInference:::importance_weight(noise_scale,
                                                     t(as.matrix(target_sample))+candidate,
                                                     t(opt_samples),
                                                     opt_transform,
                                                     target_transform,
                                                     observed_raw)
      return(mean((target_sample<observed_target[i])*weights)/mean(weights))
    }
    
    pivots[i] = pivot(0)
    print(pivots[i])
  }
  
  return(pivots)
}

collect_instances = function(n,p,s, nsim=1){
  sample_pivots = NULL
  for (i in 1:nsim){
    result = run_instance(n,p,s)
    sample_pivots = c(sample_pivots, result)
  }
  jpeg('pivots.jpg')
  plot(ecdf(sample_pivots), xlim=c(0,1),  main="Empirical CDF of null p-values", xlab="p-values", ylab="ecdf")
  abline(0, 1, lty=2)
  dev.off()
}


collect_instances(n=100, p=20, s=0)



