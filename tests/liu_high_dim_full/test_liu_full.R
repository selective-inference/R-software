library(gglasso)
library(MASS)
library(selectiveInference)
library(glmnet)

# testing Liu et al type=full in high dimensional settings -- uses debiasing matrix

test_liu_full = function(seed=1, outfile=NULL, family="gaussian", lambda_frac=0.7,
                         nrep=50, n=200, p=500, s=20, rho=0.){
  
  snr = sqrt(2*log(p)/n)
  
  set.seed(seed)
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
    
    if (family=="gaussian"){
      sigma=1
      data = selectiveInference:::gaussian_instance(n=n, p=p, s=s, rho=rho, sigma=sigma, snr=snr)
      loss = 'ls'
    } else if (family=='binomial'){
      sigma=1
      data = selectiveInference:::logistic_instance(n=n, p=p, s=s, rho=rho, snr=snr)
      loss = 'logit'
    }

    X=data$X
    y=data$y
    beta=data$beta
    cat("true nonzero:", which(beta!=0), "\n")
    
    # CV = cv.glmnet(X, y, standardize=FALSE, intercept=FALSE, family=selectiveInference:::family_label(loss))
    # sigma_est=selectiveInference:::estimate_sigma(X,y,coef(CV, s="lambda.min")[-1]) # sigma via Reid et al.
    sigma_est=sigma
    #sigma_est = selectiveInference:::estimate_sigma_data_spliting(X,y)
    print(c("sigma est", sigma_est))
    
    lambda = lambda_frac*selectiveInference:::theoretical.lambda(X, loss, sigma_est)  # theoretical lambda
    print(c("lambda", lambda))
    
    soln = selectiveInference:::solve_problem_glmnet(X, y, lambda, penalty_factor=penalty_factor, family=family)
    PVS = ROSI(X, 
               y, 
               soln, 
               lambda=lambda*n, 
               penalty_factor=penalty_factor, 
               dispersion=sigma_est^2, 
               family=family,
               solver="QP", 
               construct_ci=construct_ci, 
               debiasing_method="JM",
               verbose=TRUE)
    
    active_vars=PVS$active_set
    cat("active_vars:",active_vars,"\n")
    pvalues = c(pvalues, PVS$pvalues)

    naive_Z = PVS$estimate / PVS$std_err
    naive_P = pnorm(naive_Z)
    naive_P = 2 * pmin(naive_P, 1 - naive_P)
    naive_pvalues = c(naive_pvalues, naive_P)
    sel_intervals = rbind(sel_intervals, PVS$intervals)  # matrix with two rows
    naive_Q = qnorm(0.95)
    naive_int = cbind(PVS$estimate - naive_Q * PVS$std_err, PVS$estimate + naive_Q * PVS$std_err)
    naive_int[is.na(PVS$pvalues),] = NA
    naive_intervals = rbind(naive_intervals, naive_int)
    print('naive intervals')
    print(naive_intervals)
    if (length(pvalues)>0){
      plot(ecdf(pvalues))
      lines(ecdf(naive_pvalues), col="red")
      abline(0,1)
    }
    
    if (construct_ci && length(active_vars)>0){
      
      sel_coverages=c(sel_coverages, selectiveInference:::compute_coverage(PVS$intervals, beta[active_vars]))
      naive_coverages=c(naive_coverages, selectiveInference:::compute_coverage(naive_int, beta[active_vars]))
      sel_lengths=c(sel_lengths, as.vector(naive_int[,2]-naive_int[,1]))
      naive_lengths=c(naive_lengths, as.vector(PVS$naive_intervals[,2]-PVS$naive_intervals[,1]))
      #cat("sel cov", sel_coverages, "\n")
      print(c("selective coverage:", mean(sel_coverages, na.rm=TRUE)))
      NA_sel = is.na(sel_coverages)
      naive_coverages[NA_sel] = NA
      print(c("naive coverage:", mean(naive_coverages, na.rm=TRUE)))
      print(c("selective length mean:", mean(sel_lengths, na.rm=TRUE)))
      print(c("selective length median:", median(sel_lengths, na.rm=TRUE)))
      naive_lengths[NA_sel] = NA
      print(c("naive length mean:", mean(naive_lengths, na.rm=TRUE)))
      print(c("naive length median:", median(naive_lengths, na.rm=TRUE)))
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

test_liu_full()

