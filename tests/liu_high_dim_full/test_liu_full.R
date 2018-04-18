library(gglasso)
library(MASS)
library(selectiveInference)
library(glmnet)

# testing Liu et al type=full in high dimensional settings -- uses debiasing matrix

test_liu_full = function(seed=1, outfile=NULL, 
                         nrep=5, n=500, p=5000, s=30, rho=0.){
  
  snr = sqrt(2*log(p)/n)
  
  set.seed(seed)
  loss="ls"
  construct_ci=TRUE
  penalty_factor = rep(1, p)
  
  pvalues=NULL
  sel_intervals=NULL
  sel_coverages=NULL
  sel_lengths=NULL
  naive_pvalues=NULL
  naive_intervals=NULL
  naive_coverages=NULL
  naive_lengths=NULL
  
  FDR_sample = NULL
  power_sample=NULL
  
  for (i in 1:nrep){
    data = selectiveInference:::gaussian_instance(n=n, p=p, s=s, rho=rho, sigma=1, snr=snr)
    X=data$X
    y=data$y
    beta=data$beta
    cat("true nonzero:", which(beta!=0), "\n")
    
    # CV = cv.glmnet(X, y, standardize=FALSE, intercept=FALSE, family=selectiveInference:::family_label(loss))
    # sigma_est=selectiveInference:::estimate_sigma(X,y,coef(CV, s="lambda.min")[-1]) # sigma via Reid et al.
    sigma_est=1
    print(c("sigma est", sigma_est))
    
    # lambda = CV$lambda[which.min(CV$cvm+rnorm(length(CV$cvm))/sqrt(n))]  # lambda via randomized cv 
    lambda = 0.8*selectiveInference:::theoretical.lambda(X, loss, sigma_est)  # theoretical lambda
    print(c("lambda", lambda))
    
    soln = selectiveInference:::solve_problem_glmnet(X, y, lambda, penalty_factor=penalty_factor, loss=loss)
    #soln = solve_problem_gglasso(X, y, groups=1:ncol(X), lambda, penalty_factor=penalty_factor, loss=loss)
    PVS = selectiveInference:::inference_group_lasso(X, y, soln, groups=1:ncol(X), lambda=lambda, penalty_factor=penalty_factor, 
                                sigma_est, loss=loss, algo="glmnet", construct_ci = construct_ci)
    active_vars=PVS$active_vars
    cat("active_vars:",active_vars,"\n")
    pvalues = c(pvalues, PVS$pvalues)
    naive_pvalues = c(naive_pvalues, PVS$naive_pvalues)
    sel_intervals = cbind(sel_intervals, PVS$sel_intervals)  # matrix with two rows
    naive_intervals = cbind(naive_intervals, PVS$naive_intervals)
    
    if (length(pvalues)>0){
      plot(ecdf(pvalues))
      lines(ecdf(naive_pvalues), col="red")
      abline(0,1)
    }
    
    if (construct_ci && length(active_vars)>0){
      
      sel_coverages=c(sel_coverages, selectiveInference:::compute_coverage(PVS$sel_intervals, beta[active_vars]))
      naive_coverages=c(naive_coverages, selectiveInference:::compute_coverage(PVS$naive_intervals, beta[active_vars]))
      sel_lengths=c(sel_lengths, as.vector(PVS$sel_intervals[2,]-PVS$sel_intervals[1,]))
      naive_lengths=c(naive_lengths, as.vector(PVS$naive_intervals[2,]-PVS$naive_intervals[1,]))
      print(c("selective coverage:", mean(sel_coverages)))
      print(c("naive coverage:", mean(naive_coverages)))
      print(c("selective length mean:", mean(sel_lengths)))
      print(c("selective length median:", median(sel_lengths)))
      print(c("naive length mean:", mean(naive_lengths)))
      print(c("naive length median:", median(naive_lengths)))
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
    outfile="liu_full.rds"
  }
  
  saveRDS(list(sel_intervals=sel_intervals, sel_coverages=sel_coverages, sel_lengths=sel_lengths,
               naive_intervals=naive_intervals, naive_coverages=naive_coverages, naive_lengths=naive_lengths,
               pvalues=pvalues, naive_pvalues=naive_pvalues,
               FDR_sample=FDR_sample, power_sample=power_sample,
               n=n, p=p, s=s, snr=snr, rho=rho), file=outfile)
  
  return(list(pvalues=pvalues, naive_pvalues=naive_pvalues))
}

#test_liu_full()

