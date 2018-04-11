library(MASS)
library(selectiveInference)
library(glmnet)

loss="ls"
n=100
p=50
s=0
snr=sqrt(2*log(p)/n)
rho=0.
ridge_term = 1/sqrt(n)
noise_scale=1
sigma_est=1
lambda=1/n

set.seed(1)
data = selectiveInference:::gaussian_instance(n=n, p=p, s=s, rho=rho, sigma=1, snr=snr, scale=TRUE)

X=data$X
y=data$y
beta=data$beta

#print(X)
#print(y)

rand_lasso_soln = selectiveInference:::randomizedLasso(X, 
                                                       y, 
                                                       lambda*n, 
                                                       family=selectiveInference:::family_label(loss),
                                                       condition_subgrad=TRUE,
                                                       ridge_term=ridge_term,
                                                       noise_scale=noise_scale)

print(rand_lasso_soln$soln)
active_set = which(rand_lasso_soln$soln!=0)
LM = lm(y~X[, active_set])
sigma_est = summary(LM)$sigma
print(c("sigma est", sigma_est))
targets=selectiveInference:::compute_target(rand_lasso_soln, type="partial", sigma_est=sigma_est)

PVS = selectiveInference:::randomizedLassoInf(rand_lasso_soln,
                                              targets=targets,
                                              sampler ="norejection", #"adaptMCMC", # 
                                              level=0.9, 
                                              burnin=1000, 
                                              nsample=10000)

active_vars=rand_lasso_soln$active_set
cat("active_vars:",active_vars,"\n")
print(c("nactive", length(active_vars)))
pvalues = PVS$pvalues
sel_intervals =  t(PVS$ci) 

print(pvalues)


