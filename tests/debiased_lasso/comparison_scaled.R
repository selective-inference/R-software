source('javanmard_montanari.R')

##############################################

# Runs nsims simulations under the global null, computing p-values
# using both the old code (slow one using Adel's code) and the new
# code (faster using Jon's code), and produces qq-plots for both.
# Runing 50 sims takes about 10-15 mins because old code is slow, so
# feel free to lower nsims if you want


library(selectiveInference)
library(glmnet)

# set.seed(424)

n=100
p=200

sigma=.5

theor_lambda = sigma * sqrt(2 * log(p))
lambda=c(0.25, 0.5, 1, 0.8 * theor_lambda, theor_lambda)

for (j in c(3,4,5,1,2)) {

thresh = 1e-10

beta=rep(0,p)
type="full"

nsim = 20

scaling = sqrt(n)
pvs_old = c()
pvs_new <- c()
pvs_old_0 = c()   # don't add the offset correction
pvs_new_0 = c()   # don't add the offset correction
for (i in 1:nsim) {
  cat(i,fill=T)
  x = matrix(rnorm(n*p),n,p)
  x = scale(x,T,T) / scaling
  mu = x%*%beta 
  y=mu+sigma*rnorm(n)

  # first run  glmnet
  gfit=glmnet(x,y,intercept=F,standardize=F,thresh=thresh)

  bhat = coef(gfit, s=lambda[j]/(sqrt(n) * scaling), exact=TRUE,x=x,y=y)[-1]

  if(sum(bhat != 0) > 0) {

  # compute fixed lambda p-values and selection intervals

  aa = fixedLassoInf(x,y,bhat,lambda[j]*sqrt(n) / scaling,intercept=F,sigma=sigma,type=type)
  bb = oldFixedLassoInf(x,y,bhat,lambda[j]*sqrt(n) / scaling,intercept=F,sigma=sigma,type=type)
  cc = fixedLassoInf(x,y,bhat,lambda[j]*sqrt(n) / scaling,intercept=F,sigma=sigma,type=type, offset_correction=FALSE)
  dd = oldFixedLassoInf(x,y,bhat,lambda[j]*sqrt(n) / scaling,intercept=F,sigma=sigma,type=type, offset_correction=FALSE)
  pvs_new <- c(pvs_new, aa$pv, recursive=TRUE)
  pvs_old <- c(pvs_old, bb$pv,recursive=TRUE)
  pvs_new_0 <- c(pvs_new_0, cc$pv, recursive=TRUE)
  pvs_old_0 <- c(pvs_old_0, dd$pv, recursive=TRUE)

  cat()
  }
}

#check uniformity 

png(paste('comparison_scaled', j, '.png', sep=''))
plot(ecdf(pvs_old), pch=23, col='green', xlim=c(0,1), ylim=c(0,1), main='ECDF of p-values')
plot(ecdf(pvs_new), pch=24, col='purple', add=TRUE)
plot(ecdf(pvs_old_0), pch=23, col='red', add=TRUE)
plot(ecdf(pvs_new_0), pch=24, col='black', add=TRUE)
abline(0,1)
legend("bottomright", legend=c("Old","New", "Old 0", "New 0"), pch=c(23,24,23,24), pt.bg=c("green","purple","red","black"))
dev.off()
}