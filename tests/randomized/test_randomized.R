library(MASS)
library(selectiveInference)
library(glmnet)


test_randomized = function(seed=1, outfile=NULL, type="partial", loss="ls", lambda_frac=0.7,
                           nrep=50, n=200, p=800, s=30, rho=0.){
  
  snr = sqrt(2*log(p)/n)
  
  set.seed(seed)
  construct_ci=TRUE
  penalty_factor = rep(1, p)
  
  pvalues = NULL
  sel_intervals=NULL
  sel_coverages=NULL
  sel_lengths=NULL
  
  FDR_sample = NULL
  power_sample=NULL
  
  for (i in 1:nrep){
    
    if (loss=="ls"){
      data = selectiveInference:::gaussian_instance(n=n, p=p, s=s, rho=rho, sigma=1, snr=snr)
    } else if (loss=="logit"){
      data = selectiveInference:::logistic_instance(n=n, p=p, s=s, rho=rho, snr=snr)
    }
    
    X=data$X
    y=data$y
    beta=data$beta
    cat("true nonzero:", which(beta!=0), "\n")
    
    #CV = cv.glmnet(X, y, standardize=FALSE, intercept=FALSE, family=selectiveInference:::family_label(loss))
    #sigma_est=selectiveInference:::estimate_sigma(X,y,coef(CV, s="lambda.min")[-1]) # sigma via Reid et al.
    #sigma_est=1
    sigma_est = selectiveInference:::estimate_sigma_data_spliting(X,y)
    print(c("sigma est", sigma_est))
    
    # lambda = CV$lambda[which.min(CV$cvm+rnorm(length(CV$cvm))/sqrt(n))]  # lambda via randomized cv 
    lambda = lambda_frac*selectiveInference:::theoretical.lambda(X, loss, sigma_est)  # theoretical lambda
    
    
    rand_lasso_soln = selectiveInference:::randomizedLasso(X, 
                                                           y, 
                                                           lambda*n, 
                                                           family=selectiveInference:::family_label(loss),
                                                           condition_subgrad=TRUE)
    
    targets=selectiveInference:::compute_target(rand_lasso_soln, type=type, sigma_est=sigma_est)
    
    PVS = selectiveInference:::randomizedLassoInf(rand_lasso_soln,
                                                  targets=targets,
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
    outfile=paste("randomized_", type, ".rds", sep="")
  }
  
  saveRDS(list(sel_intervals=sel_intervals, sel_coverages=sel_coverages, sel_lengths=sel_lengths,
               pvalues=pvalues,
               FDR_sample=FDR_sample, power_sample=power_sample,
               n=n,p=p, s=s, snr=snr, rho=rho, type=type), file=outfile)
  
  return(list(pvalues=pvalues))
}

test_randomized()


