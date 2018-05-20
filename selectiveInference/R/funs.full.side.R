
estimate_sigma_data_spliting  = function(X,y, verbose=FALSE){
  nrep = 20
  sigma_est = 0
  nest = 0
  for (i in 1:nrep){
    n=nrow(X)
    m=floor(n/2)
    subsample = sample(1:n, m, replace=FALSE)
    leftover = setdiff(1:n, subsample)
    CV = cv.glmnet(X[subsample,], y[subsample], standardize=FALSE, intercept=FALSE, family="gaussian")
    beta_hat = coef(CV, s="lambda.min")[-1]
    selected = which(beta_hat!=0)
    if (verbose){
      print(c("nselected",length(selected)))
    }
    if (length(selected)>0){
      LM = lm(y[leftover]~X[leftover,][,selected])
      sigma_est = sigma_est+sigma(LM)
      nest = nest+1
    }
  }
  return(sigma_est/nest)
}

selective.plus.BH = function(beta, selected.vars, pvalues, q, verbose=FALSE){
  
  if (is.null(selected.vars)){
    return(list(power=NA, FDR=NA, pvalues=NULL, null.pvalues=NULL, ci=NULL, nselected=0))
  }
  
  nselected = length(selected.vars)
  p.adjust.BH = p.adjust(pvalues, method = "BH", n = nselected)
  rejected = selected.vars[which(p.adjust.BH<q)]
  nrejected=length(rejected)
  
  if (verbose){
      print(paste("sel+BH rejected", nrejected, "vars:",toString(rejected)))
  }
  
  true.nonzero = which(beta!=0)
  true.nulls = which(beta==0)
  if (verbose){
    print(paste("true nonzero", length(true.nonzero), "vars:", toString(true.nonzero)))
  }
  
  TP = length(intersect(rejected, true.nonzero))
  s = length(true.nonzero)
  if (s==0){
    power = NA
  } else{
    power = TP/s
  }
  
  FDR = (nrejected-TP)/max(1, nrejected)
  
  selected.nulls = NULL
  for (i in 1:nselected){
    if (any(true.nulls==selected.vars[i])){
      selected.nulls = c(selected.nulls, i)
    }
  }
  
  null.pvalues=NA
  if (length(selected.nulls)>0){
    null.pvalues = pvalues[selected.nulls]
    if (verbose){
      print(paste("selected nulls", length(selected.nulls), "vars:",toString(selected.vars[selected.nulls])))
    }
  }
  
  return(list(power=power, 
              FDR=FDR, 
              pvalues=pvalues, 
              null.pvalues=null.pvalues, 
              nselected=nselected, 
              nrejected=nrejected))
}


AR_design = function(n, p, rho, scale=FALSE){
  times = c(1:p)
  cov_mat <- rho^abs(outer(times, times, "-"))
  chol_mat = chol(cov_mat) # t(chol_mat) %*% chol_mat = cov_mat
  X=matrix(rnorm(n*p), nrow=n) %*% t(chol_mat)
  if (scale==TRUE){
    X = scale(X)
    X = X/sqrt(n)
  }
  return(X)
}


equicorrelated_design = function(n, p, rho, scale=FALSE){
  X = sqrt(1-rho)*matrix(rnorm(n*p),n) + sqrt(rho)*matrix(rep(rnorm(n), p), nrow = n)
  if (scale==TRUE){
    X = scale(X)
    X = X/sqrt(n)
  }
  return(X)
}

gaussian_instance = function(n, p, s, rho, sigma, snr, random_signs=TRUE, scale=FALSE, design="AR"){
  
  if (design=="AR"){
    X=AR_design(n,p,rho, scale)
  } else if (design=="equicorrelated"){
    X=equicorrelated_design(n,p, rho, scale)
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

logistic_instance = function(n, p, s, rho, snr, random_signs=TRUE, scale=FALSE, design="AR"){
  
  if (design=="AR"){
    X=AR_design(n,p,rho, scale)
  } else if (design=="equicorrelated"){
    X=equicorrelated_design(n,p, rho, scale)
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


estimate_sigma = function(X, y, beta_hat_cv){
  n=nrow(X)
  p=ncol(X)
  if (n<p){
    # Reid et al. (2016) "A study of error variance estimation in lasso regression"
    residuals = y-X%*% beta_hat_cv
    nactive=length(which(beta_hat_cv!=0))
    sigma_est_sq = drop(t(residuals) %*% residuals/(nrow(X)-nactive))
    sigma_est=sqrt(sigma_est_sq)
  } else{
    m = lm(y~X-1)
    sigma_est = summary(m)$sigma
  }
  return(sigma_est)
}



# theoretical lambda for glmnet
# glmnet solves: 1/2n*\|y-X\beta\|_2^2+\lambda_1\|\beta\|_1
theoretical.lambda = function(X, loss="ls", sigma=1){
  n = nrow(X); p = ncol(X)
  nsamples= 1000
  if (loss=="ls"){
    empirical = apply(abs(t(X) %*% matrix(rnorm(n*nsamples), nrow=n)), 2, max)
  } else if (loss=="logit"){
    empirical = apply(abs(t(X) %*% matrix(sample(c(-0.5,0.5), n*nsamples, replace = TRUE), nrow=n)), 2, max)
  }
  lam = mean(empirical)*sigma/n
  return(lam)
}
