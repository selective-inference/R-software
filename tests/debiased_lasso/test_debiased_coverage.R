library(selectiveInference)
library(glmnet)
library(MASS)


debiased_lasso_inference=function(X, y, soln, loss, sigma_est, lambda, debias_mat){
  
  n=nrow(X)
  p=ncol(X)
  fit = X%*%soln
  if (loss=="ls"){
    diagonal = rep(1,n)
    W_root=diag(as.vector(diagonal))
    residuals = y-fit
  }
  if (loss=="logit"){
    diagonal = exp(fit/2)/(1+exp(fit))  ## sqrt(pi*(1-pi))
    W_root = diag(as.vector(diagonal))
    residuals = y-(exp(fit)/(1+exp(fit)))
  } 
  
  if (debias_mat == "JM"){
    M = selectiveInference:::approximate_JM(W_root %*% X, 1:p)
    estimator = soln + M %*% diag(as.vector(1/diagonal)) %*% residuals
  } else if (debias_mat=="BN"){
    M = selectiveInference:::approximate_BN(X, 1:p)
    estimator = M %*% y-((M %*% X-diag(p)) %*% soln)
  }  
  covariance = sigma_est^2*M %*% t(M)
  
  
  naive_pvalues=NULL
  naive_intervals = NULL
  
  for (i in 1:p){
    naive_int = selectiveInference:::naive_CI(estimator[i], covariance[i,i])
    naive_intervals = cbind(naive_intervals, naive_int)
    naive_pval=selectiveInference:::pvalue_naive_linear(estimator[i], covariance[i,i])
    naive_pvalues = c(naive_pvalues, naive_pval)
  }
  return(list(naive_pvalues=naive_pvalues, naive_intervals=naive_intervals))  
}


compute_coverage = function(ci, beta){
  nactive=length(beta)
  coverage_vector = rep(0, nactive)
  for (i in 1:nactive){
    if (beta[i]>=ci[1,i] && beta[i]<=ci[2,i]){
      coverage_vector[i]=1
    } else if (beta[i]<ci[1,i]){
      coverage_vector[i]=-1  # parameter of the left of the iterval
    }
      
  }
  return(coverage_vector)
}


test_debiased_coverage = function(seed=1, outfile=NULL, debias_mat = "BN", loss="ls", lambda_frac=0.8,
                         nrep=20, n=200, p=500, s=30, rho=0.){
  
  snr=sqrt(2*log(p)/n)

  set.seed(seed)
  construct_ci=TRUE
  penalty_factor = rep(1, p)
  
  pvalues=NULL
  naive_pvalues=NULL
  naive_intervals=NULL
  naive_coverages=NULL
  naive_lengths=NULL
  
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
    
    # CV = cv.glmnet(X, y, standardize=FALSE, intercept=FALSE, family=selectiveInference:::family_label(loss))
    # sigma_est=selectiveInference:::estimate_sigma(X,y,coef(CV, s="lambda.min")[-1]) # sigma via Reid et al.
    sigma_est=1
    #sigma_est = selectiveInference:::estimate_sigma_data_spliting(X,y)
    print(c("sigma est", sigma_est))
    
    # lambda = CV$lambda[which.min(CV$cvm+rnorm(length(CV$cvm))/sqrt(n))]  # lambda via randomized cv 
    lambda = lambda_frac*selectiveInference:::theoretical.lambda(X, loss, sigma_est)  # theoretical lambda
    print(c("lambda", lambda))
    
    soln = selectiveInference:::solve_problem_glmnet(X, y, lambda, penalty_factor=penalty_factor, loss=loss)
    print(c("nactive", length(which(soln!=0))))
    cat("selected", which(soln!=0), "\n")
    
    active_set = which(beta!=0) #1:p 
    
    PVS = debiased_lasso_inference(X,y,soln,loss=loss, sigma_est=sigma_est, lambda=lambda, debias_mat = debias_mat)
    
    naive_pvalues = c(naive_pvalues, PVS$naive_pvalues[active_set])
    naive_intervals = cbind(naive_intervals, PVS$naive_intervals[, active_set])
    
    if (length(pvalues)>0){
      lines(ecdf(naive_pvalues), col="red")
      abline(0,1)
    }
    
    if (construct_ci){
      
      naive_coverages=c(naive_coverages, compute_coverage(PVS$naive_intervals[, active_set], beta[active_set]))
      naive_lengths=c(naive_lengths, as.vector(PVS$naive_intervals[2,active_set]-PVS$naive_intervals[1,active_set]))
      print(c("naive coverage:", length(which(naive_coverages==1))/length(naive_coverages)))
      print(c("param on the left of CI:", length(which(naive_coverages==-1))/length(naive_coverages)))
      print(c("param on the right of CI:", length(which(naive_coverages==0))/length(naive_coverages)))
      
      print(c("naive length mean:", mean(naive_lengths)))
      print(c("naive length median:", median(naive_lengths)))
    }
    
    mc = selectiveInference:::selective.plus.BH(beta, 1:p, PVS$naive_pvalues[active_set], q=0.2)
    FDR_sample=c(FDR_sample, mc$FDR)
    power_sample=c(power_sample, mc$power)
    
    if (length(FDR_sample)>0){
      print(c("FDR:", mean(FDR_sample)))
      print(c("power:", mean(power_sample)))
    }
  }
  
  if (is.null(outfile)){
    outfile="debiased_coverage.rds"
  }
  
  saveRDS(list(naive_intervals=naive_intervals, naive_coverages=naive_coverages, naive_lengths=naive_lengths,
               naive_pvalues=naive_pvalues,
               FDR_sample=FDR_sample, power_sample=power_sample,
               n=n, p=p, s=s, snr=snr, rho=rho), file=outfile)
  
  return(NULL)
}

test_debiased_coverage()

