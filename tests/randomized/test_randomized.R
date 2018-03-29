
smoke_test = function() {
    n = 100; p = 50
    X = matrix(rnorm(n * p), n, p)
    y = rnorm(n)
    lam = 20 / sqrt(n)
    noise_scale = 0.01 * sqrt(n)
    ridge_term = .1 / sqrt(n)
    selectiveInference:::randomizedLasso(X, 
                                         y, 
                                         lam, 
                                         noise_scale=noise_scale, 
                                         ridge_term=ridge_term)
}

A = smoke_test()

sampler_test = function() {

    n = 100; p = 50
    X = matrix(rnorm(n * p), n, p)
    y = rnorm(n)
    lam = 20 / sqrt(n)
    noise_scale = 0.01 * sqrt(n)
    ridge_term = .1 / sqrt(n)
    obj = selectiveInference:::randomizedLasso(X, 
                                               y, 
                                               lam, 
                                               noise_scale=noise_scale, 
                                               ridge_term=ridge_term,
 					       condition_subgrad=FALSE)
    S = selectiveInference:::sample_opt_variables(obj$law, jump_scale=rep(1/sqrt(n), p), nsample=10000)
    return(S$samples[2001:10000,])
}
B = sampler_test()

gaussian_density_test = function() {

    noise_scale = 10.
    random_lasso = smoke_test()
    p = nrow(random_lasso$internal_transform$linear_term)
    internal_state = matrix(rnorm(p * 20), p, 20)
    optimization_state = matrix(rnorm(p * 20), p, 20)
    offset = rnorm(p)

    V1 = selectiveInference:::log_density_gaussian_(noise_scale,
                                                    random_lasso$internal_transform$linear_term,
                                                    internal_state,
                                                    random_lasso$optimization_transform$linear_term,
                                                    optimization_state,
                                                    offset)
    A1 = random_lasso$internal_transform$linear_term
    A2 = random_lasso$optimization_transform$linear_term
    arg = A1 %*% internal_state + A2 %*% optimization_state + offset
    V2 = -apply(arg^2, 2, sum) / (2 * noise_scale^2)
    print(sqrt(sum((V1-V2)^2) / sum(V1^2)))

    U1 = selectiveInference:::log_density_gaussian_conditional_(noise_scale,
                                                                random_lasso$optimization_transform$linear_term,
                                                                optimization_state,
                                                                offset)
    arg = A2 %*% optimization_state + offset
    U2 = -apply(arg^2, 2, sum) / (2 * noise_scale^2)
    print(sqrt(sum((U1-U2)^2) / sum(U1^2)))

    # test that a single column matrix works -- numeric should not

    print(selectiveInference:::log_density_gaussian_conditional_(noise_scale,
                                                                 random_lasso$optimization_transform$linear_term,
                                                                 optimization_state[,1,drop=FALSE],
                                                                 offset))
    print(selectiveInference:::log_density_gaussian_(noise_scale,
                                                     random_lasso$internal_transform$linear_term,
                                                     internal_state[,1,drop=FALSE],
                                                     random_lasso$optimization_transform$linear_term,
                                                     optimization_state[,1,drop=FALSE],
                                                     offset))

}

gaussian_density_test()

laplace_density_test = function() {

    noise_scale = 10.
    random_lasso = smoke_test()
    p = nrow(random_lasso$internal_transform$linear_term)
    internal_state = matrix(rnorm(p * 20), p, 20)
    optimization_state = matrix(rnorm(p * 20), p, 20)
    offset = rnorm(p)

    V1 = selectiveInference:::log_density_laplace_(noise_scale,
                                                    random_lasso$internal_transform$linear_term,
                                                    internal_state,
                                                    random_lasso$optimization_transform$linear_term,
                                                    optimization_state,
                                                    offset)
    A1 = random_lasso$internal_transform$linear_term
    A2 = random_lasso$optimization_transform$linear_term
    arg = A1 %*% internal_state + A2 %*% optimization_state + offset
    V2 = -apply(abs(arg), 2, sum) / noise_scale
    print(sqrt(sum((V1-V2)^2) / sum(V1^2)))

    U1 = selectiveInference:::log_density_laplace_conditional_(noise_scale,
                                                                random_lasso$optimization_transform$linear_term,
                                                                optimization_state,
                                                                offset)
    arg = A2 %*% optimization_state + offset
    U2 = -apply(abs(arg), 2, sum) / noise_scale
    print(sqrt(sum((U1-U2)^2) / sum(U1^2)))

    # test that a single column matrix works -- numeric should not

    print(selectiveInference:::log_density_laplace_conditional_(noise_scale,
                                                                 random_lasso$optimization_transform$linear_term,
                                                                 optimization_state[,1,drop=FALSE],
                                                                 offset))
    print(selectiveInference:::log_density_laplace_(noise_scale,
                                                     random_lasso$internal_transform$linear_term,
                                                     internal_state[,1,drop=FALSE],
                                                     random_lasso$optimization_transform$linear_term,
                                                     optimization_state[,1,drop=FALSE],
                                                     offset))


test_randomized = function(seed=1, outfile, type="full", nrep=1, n=1000, p=10000, s=50, rho=0.){
  
  snr = sqrt(2*log(p)/n)
  
  set.seed(seed)
  loss="ls"
  construct_ci=TRUE
  penalty_factor = rep(1, p)
  
  pvalues = NULL
  sel_intervals=NULL
  sel_coverages=NULL
  sel_lengths=NULL
  
  FDR_sample = NULL
  power_sample=NULL
  
  for (i in 1:nrep){
    data = selectiveInference:::gaussian_instance(n=n, p=p, s=s, rho=rho, sigma=1, snr=snr)
    X=data$X
    y=data$y
    beta=data$beta
    cat("true nonzero:", which(beta!=0), "\n")
    
    #CV = cv.glmnet(X, y, standardize=FALSE, intercept=FALSE, family=selectiveInference:::family_label(loss))
    #sigma_est=selectiveInference:::estimate_sigma(X,y,coef(CV, s="lambda.min")[-1]) # sigma via Reid et al.
    #print(c("sigma est", sigma_est))
    sigma_est=1
    # lambda = CV$lambda[which.min(CV$cvm+rnorm(length(CV$cvm))/sqrt(n))]  # lambda via randomized cv 
    lambda = 0.8*selectiveInference:::theoretical.lambda(X, loss, sigma_est)  # theoretical lambda
    
    
    rand_lasso_soln = selectiveInference:::randomizedLasso(X, 
                                                           y, 
                                                           lambda*n, 
                                                           family=selectiveInference:::family_label(loss),
                                                           condition_subgrad=TRUE)
    
    full_targets=selectiveInference:::set.target(rand_lasso_soln, type=type, sigma_est=sigma_est)
    
    PVS = selectiveInference:::randomizedLassoInf(rand_lasso_soln,
                                                  full_targets=full_targets,
                                                  sampler = "norejection", #"adaptMCMC", #
                                                  level=0.9, 
                                                  burnin=1000, 
                                                  nsample=10000)
    active_vars=rand_lasso_soln$active_set
    cat("active_vars:",active_vars,"\n")
    pvalues = c(pvalues, PVS$pvalues)
    sel_intervals = cbind(sel_intervals, t(PVS$ci))  # matrix with two rows
    
    
    if (length(pvalues)>0){
      plot(ecdf(pvalues))
      #lines(ecdf(naive_pvalues), col="red")
      abline(0,1)
    }
    
    if (construct_ci && length(active_vars)>0){
      
      sel_coverages=c(sel_coverages, selectiveInference:::compute_coverage(t(PVS$ci), beta[active_vars]))
      sel_lengths=c(sel_lengths, as.vector(PVS$ci[,2]-PVS$ci[,1]))
      print(c("selective coverage:", mean(sel_coverages)))
      print(c("selective length mean:", mean(sel_lengths)))
      print(c("selective length median:", median(sel_lengths)))
      #naive_coverages=c(naive_coverages, selectiveInference:::compute_coverage(PVS$naive_intervals, beta[active_vars]))
      #naive_lengths=c(naive_lengths, as.vector(PVS$naive_intervals[2,]-PVS$naive_intervals[1,]))
      #print(c("naive coverage:", mean(naive_coverages)))
      #print(c("naive length mean:", mean(naive_lengths)))
      #print(c("naive length median:", median(naive_lengths)))
    }
    
    mc = selectiveInference:::selective.plus.BH(beta, active_vars, PVS$pvalues, q=0.2)
    FDR_sample=c(FDR_sample, mc$FDR)
    power_sample=c(power_sample, mc$power)
    
    if (length(FDR_sample)>0){
      print(c("FDR:", mean(FDR_sample)))
      print(c("power:", mean(power_sample)))
    }
  }
  
  if (is.null(outfile)){
    outfile="randomized_full.rds"
  }
  
  saveRDS(list(sel_intervals=sel_intervals, sel_coverages=sel_coverages, sel_lengths=sel_lengths,
               pvalues=pvalues,
               FDR_sample=FDR_sample, power_sample=power_sample,
               n=n,p=p, s=s, snr=snr, rho=rho), file=outfile)
  
  return(list(pvalues=pvalues))
}

#test_randomized()
