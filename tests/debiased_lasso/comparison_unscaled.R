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

lambda=c(0.25, 0.5, 1)

for (j in c(3,2,1)) {

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
  pvs_new <- c(pvs_new, aa$pv, recursive=TRUE)
  pvs_old <- c(pvs_old, bb$pv,recursive=TRUE)

  cat()
  }
}

#check uniformity 

png(paste('comparison_unscaled', j, '.png', sep=''))
plot(ecdf(pvs_old), pch=23, col='green', xlim=c(0,1), ylim=c(0,1), main='ECDF of p-values')
plot(ecdf(pvs_new), pch=24, col='purple', add=TRUE)
abline(0,1)
legend("bottomright", legend=c("Old", "New"), pch=c(23,24), pt.bg=c("green","purple"))
dev.off()
}