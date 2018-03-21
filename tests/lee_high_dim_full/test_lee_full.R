library(selectiveInference)
library(glmnet)

# testing Lee et al type=full in high dimensional settings -- uses debiasing matrix

test_lee_full = function(nrep=50, n=100, p=300, s=10, rho=0){
  
  snr = sqrt(2*log(p)/n)
  
  set.seed(1)
  loss="ls"
  construct_ci=TRUE
  penalty_factor = rep(1, p)
  
  pvalues = NULL
  sel_intervals=NULL
  sel_coverages=NULL
  sel_lengths=NULL
  
  for (i in 1:nrep){
    data = selectiveInference:::gaussian_instance(n=n, p=p, s=s, rho=rho, sigma=1, snr=snr)
    X=data$X
    y=data$y
    beta=data$beta
    cat("true nonzero:", which(beta!=0), "\n")
    
    CV = cv.glmnet(X, y, standardize=FALSE, intercept=FALSE, family=selectiveInference:::family_label(loss))
    
    sigma_est=selectiveInference:::estimate_sigma(X,y,coef(CV, s="lambda.min")[-1])  # sigma via Reid et al.
    print(c("sigma est", sigma_est))
    
    # lambda = CV$lambda[which.min(CV$cvm+rnorm(length(CV$cvm))/sqrt(n))] # lambda via randomized cv 
    lambda = 0.8*selectiveInference:::theoretical.lambda(X, loss, sigma_est) # theoretical lambda
    
    lasso = glmnet(X, y, family=selectiveInference:::family_label(loss), alpha=1, standardize=FALSE, intercept=FALSE, thresh=1e-12)
    soln = as.numeric(coef(lasso,x=X,y=y, family=selectiveInference:::family_label(loss), s=lambda, exact=TRUE))[-1]

    PVS = selectiveInference:::fixedLassoInf(X,y,soln, intercept=FALSE, lambda*n, family=selectiveInference:::family_label(loss),
                                             type="full",sigma=sigma_est)
    
    abs_soln = abs(soln)
    beta_threshold = abs_soln[order(abs_soln,decreasing=TRUE)][length(PVS$pv)]
    active_vars = which(abs_soln>=beta_threshold)
    cat("nactive:", length(active_vars), "\n")
    cat("active vars:", active_vars, "\n")
    
    pvalues = c(pvalues, PVS$pv)
    sel_intervals = cbind(sel_intervals, t(PVS$ci))

    if (length(pvalues)>0){
      plot(ecdf(pvalues))
      abline(0,1)
    }
    
    if (construct_ci && length(active_vars)>0){
      sel_coverages=c(sel_coverages, selectiveInference:::compute_coverage(t(PVS$ci), beta[active_vars]))
      sel_lengths=c(sel_lengths, as.vector(PVS$ci[,2]-PVS$ci[,1]))
      print(c("selective coverage:", mean(sel_coverages)))
      print(c("selective length mean:", mean(sel_lengths)))
      print(c("selective length median:", median(sel_lengths)))
    }
  }
  
  saveRDS(list(sel_intervals=sel_intervals, sel_coverages=sel_coverages, sel_lengths=sel_lengths, 
               n=n,p=p, s=s, snr=snr, rho=rho), file="lee_full.rds")
  
  return(list(pvalues=pvalues))
}

test_lee_full()



