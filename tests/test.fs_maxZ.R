library(selectiveInference)
options(error=dump.frames)

set.seed(0)
n = 20
p = 5
s = 3
size = 5

sigma = 1.5
x = matrix(rnorm(n*p),n,p)


b = c(sample(c(-1,1),s,replace=T)*rep(size,s),rep(0,p-s))
b=rep(0,p)
mu = x%*%b
nsim=200
pv=matrix(NA,nsim,p)
for(ii in 1:nsim){
    cat(ii)
y = mu + sigma*rnorm(n)

obj = fs(x,y,verb=T,intercept=T,norm=T, maxsteps=p)


# Sequential inference
out = fsInf_maxZ(obj,sigma=sigma, ndraw=5000, burnin=1000)
pv[ii,]=out$pv
}


par(mfrow=c(3,3))
for(j in 1:p){
plot((1:nsim)/nsim,sort(pv[,j]))
abline(0,1)
}
