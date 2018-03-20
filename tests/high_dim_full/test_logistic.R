#library(gglasso)
library(MASS)
library(selectiveInference)
library(glmnet)
#source("/Users/Jelena/GitHub Jelena/full/code/group_lasso_fns.R")


logistic_instance = function(n, p, s, sigma=1, rho=0, snr=7, 
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
  mu = X %*% beta
  prob = exp(mu)/(1+exp(mu))
  y = rbinom(n, 1, prob)
  
  result <- list(X=X,y=y,beta=beta)
  return(result)
}


test_logistic = function(nrep=10){
  
  set.seed(1)
  loss="logit"
  construct_ci=TRUE
  n=300
  p=50
  s=10
 
  snr=2/sqrt(n)
  
  penalty_factor = rep(1, p)
  
  pvalues = NULL
  naive_pvalues = NULL
  sel_intervals=NULL
  naive_intervals=NULL
  sel_coverages=NULL
  naive_coverages=NULL
  sel_lengths=NULL
  naive_lengths=NULL
  
  
  for (i in 1:nrep){
    #data = logistic_instance(n, p, s, snr=snr)
    #X=data$X
    #y=data$y
    #beta=data$beta
    X = matrix(rnorm(n*p),n)
    beta=rep(0,p)
    beta[1:s]=snr
    mu = X %*% beta
    prob = exp(mu)/(1+exp(mu))
    y=rbinom(n,1, prob)
    
    print("beta")
    print(which(beta!=0))
    
    lambda=10/n
    #lambda = 0.7*selectiveInference:::theoretical.lambda(X, loss, sigma=1)
    
    soln = selectiveInference:::solve_problem_glmnet(X, y, lambda, penalty_factor=penalty_factor, loss=loss)
    
    PVS = selectiveInference:::inference_group_lasso(X,y, soln, groups=1:ncol(X), lambda=lambda, penalty_factor=penalty_factor,
                                sigma_est=1, loss=loss, algo="glmnet", construct_ci = construct_ci)
    
    active_vars = PVS$active_vars
    print("active_vars")
    print(active_vars)
    null_active_vars = match(which(beta==0), active_vars)
    null_active_vars = null_active_vars[!is.na(null_active_vars)]
    #print("null_active_vars")
    #print(null_active_vars)
    
    if (length(null_active_vars)>0){
      pvalues = c(pvalues, PVS$pvalues[null_active_vars])
      naive_pvalues = c(naive_pvalues, PVS$naive_pvalues[null_active_vars])
      plot(ecdf(pvalues))
      lines(ecdf(naive_pvalues), col="red")
      abline(0,1)
    }
    
    if (construct_ci && length(active_vars)>0){
      
      sel_coverages=c(sel_coverages, selectiveInference:::compute_coverage(PVS$sel_intervals, beta[active_vars]))
      naive_coverages=c(naive_coverages, selectiveInference:::compute_coverage(PVS$naive_intervals, beta[active_vars]))
      
      #sel_coverages=c(sel_coverages, compute_coverage(PVS$sel_intervals[,null_active_vars], rep(0, length(null_active_vars))))
      #naive_coverages=c(naive_coverages,compute_coverage(PVS$naive_intervals[,null_active_vars], rep(0, length(null_active_vars))))
      
      sel_lengths=c(sel_lengths, mean(PVS$sel_intervals[2,]-PVS$sel_intervals[1,]))
      naive_lengths=c(naive_lengths, mean(PVS$naive_intervals[2,]-PVS$naive_intervals[1,]))
      print(c("selective coverage:", mean(sel_coverages)))
      print(c("naive coverage:", mean(naive_coverages)))
      print(c("selective length:", mean(sel_lengths)))
      print(c("naive length:", mean(naive_lengths)))
    }
  }
  
  return(list(pvalues=pvalues, naive_pvalues=naive_pvalues))
}

test_logistic()

