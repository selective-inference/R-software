library(MASS)
library(selectiveInference)
library(glmnet)

loss="ls"
n=100
p=50
s=10
snr=2*sqrt(2*log(p)/n)
rho=0.
ridge_term = 1/sqrt(n)
noise_scale=1
sigma_est=1
lambda=10/n

set.seed(1)
data = selectiveInference:::gaussian_instance(n=n, p=p, s=s, rho=rho, sigma=1, snr=snr)

X=data$X
y=data$y
beta=data$beta
print(which(beta!=0))
#lasso = glmnet(x=X, y=y, family=selectiveInference:::family_label(loss), 
#               alpha=1, standardize=FALSE, intercept=FALSE, thresh = 1e-20)
#beta_hat = coef(lasso, s=lambda)[-1]
#print(beta_hat)

# scale=False plus the modif below should be equiv. to scale=True without the modif
#X=X/sqrt(n)
#lambda=lambda/sqrt(n)
#noise_scale = noise_scale/sqrt(n)
#ridge_term = ridge_term/n

#print(X)
#print(y)

rand_lasso_soln = selectiveInference:::randomizedLasso(X, 
                                                       y, 
                                                       lambda*n, 
                                                       family=selectiveInference:::family_label(loss),
                                                       condition_subgrad=TRUE,
                                                       ridge_term=ridge_term,
                                                       noise_scale=noise_scale)

rand_lasso_soln$soln

active_set = which(rand_lasso_soln$soln!=0)
LM = lm(y~X[, active_set])
sigma_est = summary(LM)$sigma
print(c("sigma est", sigma_est))
sigma_est = 1

targets=selectiveInference:::compute_target(rand_lasso_soln, type="full", sigma_est=sigma_est)
targets$alternatives = rep("two-sided", length(active_set))


pvalues=rep(0, length(active_set))
niters=10
for (i in 1:niters){
  PVS = selectiveInference:::randomizedLassoInf(rand_lasso_soln,
                                              targets=targets,
                                              sampler ="norejection", #"adaptMCMC", # 
                                              level=0.9, 
                                              burnin=1000, 
                                              nsample=10000)
  
  active_vars=rand_lasso_soln$active_set
  cat("active_vars:",active_vars,"\n")
  print(c("nactive", length(active_vars)))
  pvalues = pvalues+PVS$pvalues
  sel_intervals =  t(PVS$ci) 
}

print(pvalues/niters)


