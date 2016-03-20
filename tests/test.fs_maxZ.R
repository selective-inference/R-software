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


 apply(pv,2,quantile)
     [,1]   [,2]   [,3]   [,4]   [,5]
0%      0 0.0000 0.0000 0.0000 0.0000
25%     0 0.0000 0.0055 0.1215 0.2123
50%     0 0.0000 0.1826 0.4028 0.4824
75%     0 0.0875 0.5638 0.7002 0.7517
100%    0 0.9860 0.9916 0.9996 0.9984


plot((1:nsim)/nsim,sort(pv[,4]))
