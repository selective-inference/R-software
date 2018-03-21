library(gglasso)
library(MASS)
library(selectiveInference)
library(glmnet)
#.libPaths("/Users/Jelena/anaconda/lib/R/library")
#.libPaths("/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
#source("/Users/Jelena/GitHub Jelena/full/code/group_lasso_fns.R")


gaussian_instance = function(n, p, s, sigma=1, rho=0, snr=7, 
                             random_signs=TRUE, scale=FALSE, center=FALSE){
  X = sqrt(1-rho)*matrix(rnorm(n*p),n) + sqrt(rho)*matrix(rep(rnorm(n), p), nrow = n)
  if (scale==TRUE){
    X = scale(X)
    X = X/sqrt(n)
  }
  beta = rep(0, p)
  beta[1:s]=snr
  if (random_signs==TRUE && s>0){
    signs = sample(c(-1,1), s, replace = TRUE)
    beta[1:s] = beta[1:s] * signs
  }
  beta=sample(beta)
  y = X %*% beta + rnorm(n)*sigma 
  result <- list(X=X,y=y,beta=beta)
  return(result)
}


test_group_lasso = function(nrep=10){
  
  set.seed(1)
  loss="ls"
  construct_ci=TRUE
  n=500
  p=100
  s=10
  ngroups=50
  lambda=60/n
  snr=0.3
  groups = rep(1:ngroups,each=p/ngroups)
  group_sizes = rep(p/ngroups,ngroups)
  penalty_factor = rep(1, length(group_sizes))
  
  
  pvalues = NULL
  naive_pvalues = NULL
  sel_intervals=NULL
  naive_intervals=NULL
  sel_coverages=NULL
  naive_coverages=NULL
  sel_lengths=NULL
  naive_lengths=NULL
  
  
  for (i in 1:nrep){
    data = gaussian_instance(n,p,s, snr=snr)
    X=data$X
    y=data$y
    beta=data$beta
    #print("beta")
    #print(which(beta!=0))
    CV = cv.glmnet(X, y, standardize=FALSE, intercept=FALSE, family=selectiveInference:::family_label(loss))
    sigma_est=selectiveInference:::estimate_sigma(X,y, coef(CV, s="lambda.min")[-1])
    print(c("sigma est", sigma_est))
    
    soln = selectiveInference:::solve_problem_gglasso(X, y, groups, lambda, penalty_factor=penalty_factor, loss=loss)
    PVS = selectiveInference:::inference_group_lasso(X,y, soln,groups, lambda=lambda, penalty_factor=penalty_factor, 
                                        sigma_est=sigma_est, loss=loss, algo="glmnet", construct_ci = construct_ci)
    
    active_vars = PVS$active_vars
    print("active_vars")
    print(active_vars)
    null_active_vars = match(which(beta==0), active_vars)
    null_active_vars = null_active_vars[!is.na(null_active_vars)]
    print("null_active_vars")
    print(null_active_vars)
    
    if (length(null_active_vars)>0){
      print(PVS$pvalues)
      print(PVS$pvalues[null_active_vars])
      print(PVS$naive_pvalues[null_active_vars])
      pvalues = c(pvalues, PVS$pvalues[null_active_vars])
      naive_pvalues = c(naive_pvalues, PVS$naive_pvalues[null_active_vars])
      plot(ecdf(pvalues))
      lines(ecdf(naive_pvalues), col="red")
      abline(0,1)
    }
    
    if (construct_ci && length(active_vars)>0){
      
      sel_coverages=c(sel_coverages, selectiveInference:::compute_coverage(PVS$sel_intervals, beta[active_vars]))
      naive_coverages=c(naive_coverages, selectiveInference:::compute_coverage(PVS$naive_intervals, beta[active_vars]))
      sel_lengths=c(sel_lengths, mean(PVS$sel_intervals[2,]-PVS$sel_intervals[1,]))
      naive_lengths=c(naive_lengths, mean(PVS$naive_intervals[2,]-PVS$naive_intervals[1,]))
      print(c("selective covarage:", mean(sel_coverages)))
      print(c("naive coverage:", mean(naive_coverages)))
      print(c("selective length:", mean(sel_lengths)))
      print(c("naive length:", mean(naive_lengths)))
    }
    
  }
  
  return(list(pvalues=pvalues, naive_pvalues=naive_pvalues))
}


test_group_lasso()
