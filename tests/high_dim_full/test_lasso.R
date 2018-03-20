library(gglasso)
library(MASS)
library(selectiveInference)
library(glmnet)
#source("/Users/Jelena/GitHub Jelena/full/code/group_lasso_fns.R")

test_high_dim_lasso = function(nrep=10){
  
  set.seed(1)
  loss="ls"
  construct_ci=TRUE
  n=100
  p=300
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
    X=matrix(rnorm(n*p), nrow=n)
    beta = rep(0,p)
    #beta[1:10]=1
    y = X%*%beta+rnorm(n)
    
    CV = cv.glmnet(X, y, standardize=FALSE, intercept=FALSE, family=selectiveInference:::family_label(loss))
    # print(CV$lambda.min)
    # lambda via randomized cv 
    # lambda = CV$lambda[which.min(CV$cvm+rnorm(length(CV$cvm))/sqrt(n))]
    # sigma via Reid et al.
    sigma_est=selectiveInference:::estimate_sigma(X,y,coef(CV, s="lambda.min")[-1])
    print(c("sigma est", sigma_est))
    # theoretical lambda
    lambda = 0.8*selectiveInference:::theoretical.lambda(X, loss, sigma_est)
    
    soln = selectiveInference:::solve_problem_glmnet(X, y, lambda, penalty_factor=penalty_factor, loss=loss)
    #soln = solve_problem_gglasso(X, y, groups=1:ncol(X), lambda, penalty_factor=penalty_factor, loss=loss)
    PVS = selectiveInference:::inference_group_lasso(X, y, soln, groups=1:ncol(X), lambda=lambda, penalty_factor=penalty_factor, 
                                sigma_est, loss=loss, algo="Q", construct_ci = construct_ci)
    active_vars=PVS$active_vars
    pvalues = c(pvalues, PVS$pvalues)
    naive_pvalues = c(naive_pvalues, PVS$naive_pvalues)
    sel_intervals = cbind(sel_intervals, PVS$sel_intervals)
    naive_intervals = cbind(naive_intervals, PVS$naive_intervals)
    
    if (length(pvalues)>0){
      plot(ecdf(pvalues))
      lines(ecdf(naive_pvalues), col="red")
      abline(0,1)
    }
    if (construct_ci && length(active_vars)>0){
      
      sel_coverages=c(sel_coverages, selectiveInference:::compute_coverage(PVS$sel_intervals, beta[active_vars]))
      naive_coverages=c(naive_coverages, selectiveInference:::compute_coverage(PVS$naive_intervals, beta[active_vars]))
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

test_high_dim_lasso()

