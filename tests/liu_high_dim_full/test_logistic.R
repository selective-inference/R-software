library(MASS)
library(selectiveInference)
library(glmnet)


test_logistic = function(seed=1, outfile=NULL, nrep=10, n=1000, p=2000, s=30, rho=0.){
  
  snr = sqrt(2*log(p)/n)
  snr = 2*snr
  #print(snr)
  #empirical=empirical = apply(abs(t(matrix(rnorm(n*p), nrow=n)) %*% matrix(rnorm(n*1000), nrow=n)), 2, max)
  #print(mean(empirical)/n)
  
  #empirical = apply(abs(t(matrix(rnorm(n*p), nrow=n)) %*% matrix(sample(c(0,1), n*1000, replace = TRUE), nrow=n)), 2, max)
  #snr=mean(empirical)/n
  #print(snr)
  
  set.seed(seed)
  loss="logit"
  construct_ci=TRUE
  
  penalty_factor = rep(1, p)
  
  pvalues = NULL
  naive_pvalues = NULL
  sel_intervals=NULL
  naive_intervals=NULL
  sel_coverages=NULL
  naive_coverages=NULL
  sel_lengths=NULL
  naive_lengths=NULL
  
  FDR_sample = NULL
  power_sample=NULL
  
  for (i in 1:nrep){
    data = selectiveInference:::logistic_instance(n=n, p=p, s=s, rho=rho, sigma=1, snr=snr)
    X=data$X
    y=data$y
    beta=data$beta
    
    cat("true nonzero:", which(beta!=0), "\n")
    
    lambda = 0.4*selectiveInference:::theoretical.lambda(X, loss, sigma=1)
    
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
    
    mc = selectiveInference:::selective.plus.BH(beta, active_vars, PVS$pvalues, q=0.2)
    FDR_sample=c(FDR_sample, mc$FDR)
    power_sample=c(power_sample, mc$power)
    
    if (length(FDR_sample)>0){
      print(c("FDR:", mean(FDR_sample)))
      print(c("power:", mean(power_sample)))
    }
  }
  
  if (is.null(outfile)){
    outfile="liu_logistic_full.rds"
  }
  
  saveRDS(list(sel_intervals=sel_intervals, sel_coverages=sel_coverages, sel_lengths=sel_lengths,
               naive_intervals=naive_intervals, naive_coverages=naive_coverages, naive_lengths=naive_lengths,
               pvalues=pvalues, naive_pvalues=naive_pvalues,
               FDR_sample=FDR_sample, power_sample=power_sample,
               n=n, p=p, s=s, snr=snr, rho=rho), file=outfile)
  
  
  return(list(pvalues=pvalues, naive_pvalues=naive_pvalues))
}

test_logistic()

