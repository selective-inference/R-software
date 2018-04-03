library(selectiveInference)

get_instance = function(n, p, s, sigma=1, rho=0, signal=6, family="gaussian",
                        X=NA, random_signs=TRUE, scale=TRUE, center=TRUE, seed=NA){
  if (!is.na(seed)){
    set.seed(seed)
  }
  
  if (is.na(X)){
    X = sqrt(1-rho)*matrix(rnorm(n*p),n, p) + sqrt(rho)*matrix(rep(rnorm(n), p), nrow = n)
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
  mu = X %*% beta
  if (family=="gaussian"){
    y = mu + rnorm(n)*sigma
  } else if (family=="binomial"){
    prob = exp(mu)/(1+exp(mu))
    y= rbinom(n,1, prob)
  }
  result = list(X=X,y=y,beta=beta)
  return(result)
}




test_randomized_lasso = function(n=100,p=200,s=0){
  set.seed(1)
  data = get_instance(n=n,p=p,s=s, rho=0.3, sigma=1, family="binomial")
  X=data$X
  y=data$y
  lam = 2.
  noise_scale = 0.5
  ridge_term = 1./sqrt(n)
  result = selectiveInference:::randomizedLasso(X, y, lam, noise_scale, ridge_term)
  print(result$soln)
  print(length(which(result$soln!=0)))
  print(result$observed_opt_state) # compared with python code
}

test_randomized_logistic = function(n=100,p=20,s=0){
  set.seed(1)
  data = get_instance(n=n,p=p,s=s, rho=0.3, sigma=1, family="binomial")
  X=data$X
  y=data$y
  lam = 0.5
  noise_scale = 0.5
  ridge_term = 1./sqrt(n)
  set.seed(1)
  perturb = rnorm(p)*noise_scale
  result = selectiveInference:::solve_logistic(X,y,lam, ridge_term, perturb)
  print(result$soln)
  print(length(which(result$soln!=0)))
}

#test_randomized_logistic()


test_KKT=function(){
  set.seed(1)
  n=200
  p=100
  data = gaussian_instance(n=n,p=p,s=0, rho=0.3, sigma=3)
  X=data$X
  y=data$y
  lam = 2.
  noise_scale = 0.5
  ridge_term = 1./sqrt(n)
  result = selectiveInference:::randomizedLasso(X,y,lam, noise_scale, ridge_term)
  print("check KKT")
  opt_linear = result$optimization_transform$linear_term
  opt_offset = result$optimization_transform$offset_term
  observed_opt_state=result$observed_opt_state
  #print(dim(opt_linear))
  #print(opt_offset)
  #print(result$perturb)
  print(opt_linear %*% observed_opt_state+opt_offset+result$observed_raw-result$perturb) ## should be zero
}


collect_results = function(n,p,s, nsim=100, level=0.9, 
                           family = "gaussian",
                           condition_subgrad=TRUE, 
                           type="full",
                           lam=1.2){

  rho=0.
  sigma=1
  sample_pvalues = c()
  sample_coverage = c()
  for (i in 1:nsim){
    data = get_instance(n=n,p=p,s=s, rho=rho, sigma=sigma, family=family)
    X=data$X
    y=data$y

    rand_lasso_soln = selectiveInference:::randomizedLasso(X, 
                                                           y, 
                                                           lam, 
                                                           family=family,
                                                           condition_subgrad=condition_subgrad)
    
    targets=selectiveInference:::compute_target(rand_lasso_soln,type=type)
    print(targets$construct_pvalues)    
    result = selectiveInference:::randomizedLassoInf(rand_lasso_soln,
                                                     targets=targets,
                                                     sampler="norejection", #"adaptMCMC", #
                                                     level=level, 
                                                     burnin=1000, 
                                                     nsample=5000)
    if (length(result$pvalues)>0){
      true_beta = data$beta[rand_lasso_soln$active_set]
      coverage = rep(0, nrow(result$ci))
      for (i in 1:nrow(result$ci)){
        if (result$ci[i,1]<true_beta[i] & result$ci[i,2]>true_beta[i]){
          coverage[i]=1
        }
        print(paste("ci", toString(result$ci[i,])))
      }
      sample_pvalues = c(sample_pvalues, result$pvalues)
      sample_coverage = c(sample_coverage, coverage)
      print(paste("coverage", mean(sample_coverage)))
    }
  }
  if (length(sample_coverage)>0){
    print(paste("coverage", mean(sample_coverage)))
    jpeg('pivots.jpg')
    plot(ecdf(sample_pvalues), xlim=c(0,1),  main="Empirical CDF of null p-values", xlab="p-values", ylab="ecdf")
    abline(0, 1, lty=2)
    dev.off()
  }
}

set.seed(1)
collect_results(n=100, p=20, s=0, lam=1.)

